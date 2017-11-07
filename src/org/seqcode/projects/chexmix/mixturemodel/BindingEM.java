package org.seqcode.projects.chexmix.mixturemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;


/**
 * BindingEM: run EM training with sparse prior and positional prior(s) on binding data.
 * Supports multiple binding types. hAll function is calculated for each binding type and + & - strands.
 * Each component consists of probabilistic mixture of binding subtypes. 
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class BindingEM {

	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected ChExMixConfig config;
	protected List<List<BindingSubComponents>> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;  //Assumes the same number of active+inactive components in each condition
	protected int numConditions;
	protected int trainingRound=0; //Identifier for the overall training round, used only for file names
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds; //Background models per condition
	//	EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[][]   		hitCounts;	// Hit weights
	protected int[][]      		hitPos;		// Hit positions
	protected boolean[][]  		hitPlusStr;	// Hit positive strand boolean
	protected int[]		   		hitNum;		// Number of hits in each condition 
	protected int[][]      		repIndices; // Index of replicate for the hit
	protected double[][][][][] 	hAll;		// H function values for all positions in the current window (precomputed)
	protected double[][][][][] 	h; 			// H function (binding component probability per read)
	protected double[][]  		n; 			// N function (noise component probability per read)
	protected double[][][][][] 	rBind;		// Binding component responsibilities
	protected double[][]   		rNoise;		// Noise component responsibilities
	protected double[][]   		pi;			// pi : emission probabilities for binding components
	protected double[]     		piNoise;	// pi : emission probabilities for noise components (fixed)
	protected int[][][][]  		mu;			// mu : positions of the binding subtype components
	protected int[]		   		numBindingType;// Number of binding type in each condition
	protected double[][][][] 	tau;		// tau : binding event type probabilities
	protected double []    		alphaMax;	// Maximum alpha
	protected double [][]  		betaMax;	// Maximum beta
	protected double [][][][] 	epsilonMax; // Maximum epsilon
//	protected double[][]   		motifPrior; // Motif prior (indexed by condition & base) 
	protected double[][][]   	forMotifPrior; //Motif prior for forward strand (indexed by condition & base)
	protected double[][][]   	revMotifPrior; //Motif prior for reverse strand (indexed by condition & base)
	protected TagProbabilityDensity[][] TagProbabilityDensities; //Array of binding models for convenience
	protected double[][][][][] 	lastRBind; 	//Last responsibilities (monitor convergence)
	protected double[][]   		lastPi;		//Last Pi (monitor convergence)
	protected int[][][][]  		lastMu;		//Last positions (monitor convergence)
	protected double lastLAP, LAP; 			//log-likelihood monitoring
	protected boolean plotEM=false;			//Plot the current region components
	protected Region plotSubRegion=null; 	//Sub region to plot
	protected double initCondPosPriorVar=10, finalCondPosPriorVar=0.001, currCondPosPriorVar=initCondPosPriorVar; //positional prior final Gaussian variance
	protected double numPotentialRegions;
	protected int stateEquivCount=0;
	
	/**
	 * Constructor
	 * @param c
	 * @param eMan
	 */
	public BindingEM(ChExMixConfig c, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg){
		config=c;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
	}
	
	//Accessor
	public List<List<BindingSubComponents>> getComponents(){return components;}
	
	/**
     * EM training
     *
     * Returns lists of binding components indexed by condition 
     *
     * Almost purely matrix/array operations.
     */
    public List<List<BindingSubComponents>>  train(List<List<StrandedBaseCount>> signals, 
    											  Region w, 
    											  List<NoiseComponent> noise,
    											  List<List<BindingSubComponents>> comps, 
    											  int numComp,
    											  double[][][] forMotifPrior,
    											  double[][][] revMotifPrior,
    											  int trainingRound,
    											  Region plotSubRegion){
    	components = comps;
        this.noise = noise;
        numComponents = numComp;
        this.forMotifPrior = forMotifPrior;
        this.revMotifPrior = revMotifPrior;
        this.trainingRound = trainingRound;
    	this.plotSubRegion = plotSubRegion;
    	numBindingType = new int[numConditions];	// Number of binding type in each condition	
        //Matrix initializations
        hitCounts= new double[numConditions][];	// Hit weights
    	hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	repIndices= new int[numConditions][];	    // Index of replicate for the hit
    	hAll = new double[numConditions][][][][];		// H function values for all positions in the current window (precomputed). hAll needs to be computed for components residing + and - strands
    	h= new double[numConditions][][][][]; 		// H function (binding component probability per read)
    	n= new double[numConditions][]; 			// N function (noise component probability per read)
    	rBind= new double[numConditions][][][][];		// Binding component responsibilities
    	rNoise= new double[numConditions][];		// Noise component responsibilities
    	pi = new double[numConditions][numComponents];	// pi : emission probabilities for binding components
    	piNoise = new double[numConditions];		// pi : emission probabilities for noise components (fixed)
    	alphaMax = new double[numConditions];		//Maximum alpha
    	betaMax = new double[numConditions][numComponents]; // Maximum beta
    	epsilonMax = new double[numConditions][numComponents][][]; // Maximum epsilon
        mu = new int[numConditions][numComponents][][];// mu : positions of the binding components
        tau = new double[numConditions][numComponents][][];// tau : binding event subtype probabilities
        TagProbabilityDensities = new TagProbabilityDensity[manager.getReplicates().size()][]; //Array of bindingModels for convenience
        plotEM = (plotSubRegion!=null && plotSubRegion.overlaps(w));
        //Monitor state convergence using the following last variables
        lastRBind = new double[numConditions][][][][];
        lastPi = new double[numConditions][numComponents];
        lastMu = new int[numConditions][numComponents][][];
		
        //Initializing data structures
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	numBindingType[c] = bindingManager.getNumBindingType(cond);
        	
        	//Add bindingModels to array
        	for(ControlledExperiment rep : cond.getReplicates()){
        		TagProbabilityDensities[rep.getIndex()] = new TagProbabilityDensity[numBindingType[c]];
        		for (int bt=0; bt< numBindingType[c]; bt++)
        			TagProbabilityDensities[rep.getIndex()][bt] = bindingManager.getBindingModel(rep).get(bt);
        	}
        	
        	//Set maximum alphas
        	alphaMax[c] =  config.getFixedAlpha()>0 ? config.getFixedAlpha() : 
        			config.getAlphaScalingFactor() * (double)conditionBackgrounds.get(cond).getMaxThreshold('.');       
        	        	
        	//Load Reads (merge from all replicates)
        	List<StrandedBaseCount> bases = new ArrayList<StrandedBaseCount>();
        	for(ControlledExperiment rep : cond.getReplicates())
        		bases.addAll(signals.get(rep.getIndex()));
        	int numBases = bases.size();
        	hitNum[c]=numBases;
        	
        	//Load replicate index for each read
        	repIndices[c] = new int[numBases]; 
        	int y=0, z=0;
        	for(ControlledExperiment rep : cond.getReplicates()){
        		z=0;
        		while(z<signals.get(rep.getIndex()).size()){
        			repIndices[c][y] = rep.getIndex();
        			z++; y++;
        		}
        	}
        	
            //Load read info
            double[] countc= new double[numBases];
            int[] posc= new int[numBases];
            boolean[] plusc= new boolean[numBases];
            for(int i=0;i<numBases;i++){
            	posc[i] = bases.get(i).getCoordinate();
            	plusc[i] = bases.get(i).getStrand() == '+';
                countc[i]=bases.get(i).getCount();
            }
            hitPos[c] = posc;
            hitCounts[c]=countc;
            hitPlusStr[c] = plusc;       	
            
            //Load pi for binding components
            for(int j=0;j<numComp;j++){
                BindingSubComponents comp = components.get(c).get(j);
                pi[c][j]= comp.getPi(); 
            }
            //Load pi for noise components
            piNoise[c]=noise.get(c).getPi();
            
            //Load binding component positions
            for(int j=0;j<numComp;j++){
            	mu[c][j] = new int[numBindingType[c]][2];
            	lastMu[c][j] = new int[numBindingType[c]][2];
            	for (int bt=0; bt< numBindingType[c]; bt++)
            		for (int s=0; s< 2; s++)
            			mu[c][j][bt][s] = components.get(c).get(j).getPosition();
            }
            
            //Initialize epsilonMax 
            for(int j=0;j<numComp;j++)
            	epsilonMax[c][j] = new double[numBindingType[c]][2];
            
            //Initialize tau
            for(int j=0;j<numComp;j++){
            	tau[c][j] = new double[numBindingType[c]][2];
            	for (int bt=0; bt< numBindingType[c]; bt++)
            		for (int s=0; s< 2; s++)
            			tau[c][j][bt][s] = 1.0/((double) numBindingType[c]*2); 
            }
            
            //Initialize H function for all positions in the current region
            //H function are calculated for components at positive and negative strands, and for each binding type
            double[][][][] hAllc = new double[w.getWidth()][numBases][numBindingType[c]][2];
            for(int i=0;i<numBases;i++){           	
            	for(int b=0;b<w.getWidth();b++){
            		int pos = b+w.getStart();
            		int dist = hitPos[c][i]-pos;
            		if (pos >= w.getStart() && pos < w.getEnd()){
            			for (int bt=0; bt<numBindingType[c] ; bt++){
            				if (hitPlusStr[c][i]){
            					// Case 1 : component location is at positive strand
            					hAllc[b][i][bt][0] = TagProbabilityDensities[repIndices[c][i]][bt].probability(dist, true); //Watson
            					// Case 2 : component location is at negative strand
            					hAllc[b][i][bt][1] = TagProbabilityDensities[repIndices[c][i]][bt].probability(-dist, false); //Crick
            				}else{
            					// Case 1 : component location is at positive strand
            					hAllc[b][i][bt][0] = TagProbabilityDensities[repIndices[c][i]][bt].probability(dist, false); //Crick
            					// Case 2 : component location is at negative strand
            					hAllc[b][i][bt][1] = TagProbabilityDensities[repIndices[c][i]][bt].probability(-dist, true); //Watson 
            				}}}}}   
    		hAll[c] = hAllc;
    		
            //Initialize responsibility functions
            double[][][][] hc= new double[numComp][numBases][numBindingType[c]][2];
            double[] nc = new double[numBases];
            for(int i=0;i<numBases;i++){
            	for(int j=0;j<numComp;j++){
            		for (int bt=0; bt<numBindingType[c];bt++){
            			for (int s=0; s<2 ;s++){
            				int index = mu[c][j][bt][s]-w.getStart();
            				hc[j][i][bt][s] = hAll[c][index][i][bt][s];
            			}
            		}
                }
            	nc[i] = noise.get(c).scorePosition(hitPos[c][i], repIndices[c][i]);
            }
            h[c] = hc;
            n[c] = nc;
            
            rBind[c] = new double[numComp][numBases][numBindingType[c]][2];
    		rNoise[c]= new double[numBases];
    		lastRBind[c] = new double[numComp][numBases][numBindingType[c]][2];
        }
        //End of data structure initialization
        
        
        //////////
        // Run EM steps
        //////////
        EM_MAP(w);
	
        //////////
        // re-assign EM result back to component objects
        //////////
        List<List<BindingSubComponents>> activeComponents = new ArrayList<List<BindingSubComponents>>();
        for(ExperimentCondition cond : manager.getConditions()){
        	//Binding Components
        	List<BindingSubComponents> currActiveComps = new ArrayList<BindingSubComponents>();
        	int c = cond.getIndex();     	
        	for(int j=0;j<numComp;j++){ 
        		BindingSubComponents comp = components.get(c).get(j);
        		comp.setPi(pi[c][j]);	            
        		comp.setPositions(mu[c][j]);
	            comp.setTau(tau[c][j]);
	            
	            double weightedMu=0;
	            for (int bt=0; bt< numBindingType[c]; bt++){
	            	for (int s=0; s< 2; s++){ if(tau[c][j][bt][s]>0){
	            		weightedMu += ((double) mu[c][j][bt][s])*tau[c][j][bt][s];
	            	}}}
	            int newMu = (int) Math.round(weightedMu);
	            comp.setPosition(newMu);

        		double watson_resp = 0.0; double crick_resp = 0.0;
        		for(int i=0;i<hitNum[c];i++){     
        			for (int bt=0; bt< numBindingType[c]; bt++){
        				if (hitPlusStr[c][i]){
        					watson_resp += hitCounts[c][i]*rBind[c][j][i][bt][0];
        					crick_resp += hitCounts[c][i]*rBind[c][j][i][bt][1];
        				}else{
        					crick_resp += hitCounts[c][i]*rBind[c][j][i][bt][0];
        					watson_resp += hitCounts[c][i]*rBind[c][j][i][bt][1];
        				}				
        			}
        		}
        		comp.setSumResponsibilities(watson_resp, crick_resp);
        		if(pi[c][j]>0.0){
        			currActiveComps.add(comp);        		}
        	}       	
	    	activeComponents.add(currActiveComps);
	    	//Noise Components
	    	double noise_resp = 0.0;	
            for(int i=0;i<hitNum[c];i++)
                noise_resp += hitCounts[c][i]*rNoise[c][i];
	    	noise.get(c).setSumResponsibility(noise_resp);
        }        
        setComponentResponsibilityProfiles(activeComponents, signals, rBind);

        /** test **/
    	for(ExperimentCondition cond : manager.getConditions()){
    		int c = cond.getIndex(); 
    		for(int j=0;j<activeComponents.get(c).size();j++){ 
        		BindingSubComponents comp = activeComponents.get(c).get(j);
        		boolean printRegion=false;
        		Region Reb1_1 = new Region(comp.getCoord().getGenome(), "15", 631892, 631992);
        		Region Reb1_2 = new Region(comp.getCoord().getGenome(), "15", 316518, 316618);
        		Region Reb1_3 = new Region(comp.getCoord().getGenome(), "13", 306801, 306901);
        		if (Reb1_1.contains(comp.getCoord()) || Reb1_2.contains(comp.getCoord()) || Reb1_3.contains(comp.getCoord()))
        			printRegion = true;
        		// Make arrays of responsibilities
        		double[][][] w_resp = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		double[][][] c_resp = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		for (int bt=0; bt < numBindingType[c]; bt++){
        			for (int s=0; s< 2; s++){
        				for (int i=0; i< config.MAX_BINDINGMODEL_WIDTH; i++){
        					w_resp[bt][s][i] = 0;
        					c_resp[bt][s][i] = 0;
        				}
        			}
        		}
        		for (int bt=0; bt< numBindingType[c]; bt++){
					for(int i=0;i<hitNum[c];i++){
						int offset = hitPos[c][i]-mu[c][j][bt][0]+(config.MAX_BINDINGMODEL_WIDTH/2);
						if (offset >=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							if (hitPlusStr[c][i])
								w_resp[bt][0][offset]=rBind[c][j][i][bt][0]*hitCounts[c][i];
							else
								c_resp[bt][0][offset]=rBind[c][j][i][bt][0]*hitCounts[c][i];
						}
						offset = mu[c][j][bt][1]-hitPos[c][i]+(config.MAX_BINDINGMODEL_WIDTH/2);
						if (offset >=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							if (hitPlusStr[c][i])
								c_resp[bt][1][offset]=rBind[c][j][i][bt][1]*hitCounts[c][i];
							else
								w_resp[bt][1][offset]=rBind[c][j][i][bt][1]*hitCounts[c][i];
						}
					}
				}
				
        		// Make arrays of responsibilities * model read distribution
        		double[][][] w_resp_distrib = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		double[][][] c_resp_distrib = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		double[][][] w_motif = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		double[][][] c_motif = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
        		for (int bt=0; bt < numBindingType[c]; bt++){
        			for (int s=0; s< 2; s++){
        				for (int i=0; i< config.MAX_BINDINGMODEL_WIDTH; i++){
        					w_resp_distrib[bt][s][i] = 0;
        					c_resp_distrib[bt][s][i] = 0;
        					w_motif[bt][s][i]=0;
        					c_motif[bt][s][i]=0;
        				}
        			}
        		}
        		for (int bt=0; bt< numBindingType[c]; bt++){
					for(int i=0;i<hitNum[c];i++){
						int dist = hitPos[c][i]-mu[c][j][bt][0];
						int offset = hitPos[c][i]-mu[c][j][bt][0]+(config.MAX_BINDINGMODEL_WIDTH/2);
						if (offset >=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							if (hitPlusStr[c][i])
								w_resp_distrib[bt][0][offset]=rBind[c][j][i][bt][0]*hitCounts[c][i]*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist,true);	//Watson
							else
								c_resp_distrib[bt][0][offset]=rBind[c][j][i][bt][0]*hitCounts[c][i]*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist,false);	//Crick;
						}
						dist = mu[c][j][bt][1]-hitPos[c][i];
						offset = mu[c][j][bt][1]-hitPos[c][i]+(config.MAX_BINDINGMODEL_WIDTH/2);
						if (offset >=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							if (hitPlusStr[c][i])
								c_resp_distrib[bt][1][offset]=rBind[c][j][i][bt][1]*hitCounts[c][i]*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist, false);	//Crick
							else
								w_resp_distrib[bt][1][offset]=rBind[c][j][i][bt][1]*hitCounts[c][i]*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist, true);	//Watson
						}
					}

					ExperimentCondition ec = manager.getConditions().get(c); 
					// add motif score
					if (bindingManager.getMotifs(ec)!=null){
						for (int index=0; index< 10; index++){
							int offset = mu[c][j][bt][0]-w.getStart()-5+index;
							if (offset >=0 && offset<forMotifPrior[c][bt].length){
								w_motif[bt][0][index]=forMotifPrior[c][bt][offset];
								c_motif[bt][0][index]=revMotifPrior[c][bt][offset];
							}
							
							offset = mu[c][j][bt][1]-w.getStart()-5+index;
							if (offset >=0 && offset<forMotifPrior[c][bt].length){
								w_motif[bt][1][index]=revMotifPrior[c][bt][offset];
								c_motif[bt][1][index]=forMotifPrior[c][bt][offset];
							}
						}
					}
        		}
        		printRegion= false;
        		if (printRegion){
        			System.out.println("print responsibilities for point "+comp.toString());
        			for (int bt=0; bt<numBindingType[c]; bt++){
        				for (int s=0; s< 2; s++){
        					System.out.println("binding type "+bt+" strand "+s+" tau "+comp.getTau()[bt][s]);
        					System.out.println("#watson");
        					for (int i=0; i< w_resp[bt][s].length; i++)
        						System.out.print(w_resp[bt][s][i]+",");
        					System.out.println("\n"+"#crick");
        					for (int i=0; i< c_resp[bt][s].length; i++)
        						System.out.print(c_resp[bt][s][i]+",");
        					System.out.println();
        				}
        			}
      
        			System.out.println("print responsibilities times tag probability densities "+comp.toString());
        			for (int bt=0; bt<numBindingType[c]; bt++){
        				for (int s=0; s< 2; s++){
        					System.out.println("binding type "+bt+" strand "+s+" tau "+comp.getTau()[bt][s]);
        					System.out.println("#watson");
        					for (int i=0; i< w_resp_distrib[bt][s].length; i++)
        						System.out.print(w_resp_distrib[bt][s][i]+",");
        					System.out.println("\n"+"#crick");
        					for (int i=0; i< c_resp_distrib[bt][s].length; i++)
        						System.out.print(c_resp_distrib[bt][s][i]+",");
        					System.out.println();
        				}
        			}
        			System.out.println("motif score");
        			for (int bt=0; bt<numBindingType[c]; bt++){
        				for (int s=0; s< 2; s++){
        					System.out.println("binding type "+bt+" strand "+s+" tau "+comp.getTau()[bt][s]);
        					System.out.println("#watson");
        					for (int i=0; i< w_motif[bt][s].length; i++)
        						System.out.print(w_motif[bt][s][i]+",");
        					System.out.println("\n"+"#crick");
        					for (int i=0; i< c_motif[bt][s].length; i++)
        						System.out.print(c_motif[bt][s][i]+",");
        					System.out.println();
        				}
        			}
        		}
    		}	  		
    	} 	
		
        
        return activeComponents;
    }//end of EMTrain method


    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void EM_MAP (Region currRegion) {
        int numComp = numComponents;
        double [][] totalResp = new double[numConditions][];
        int regStart = currRegion.getStart();
        int regEnd = currRegion.getEnd();
        
        //Variables for tracking mu maximization. Defined early to avoid memory assignment during main EM loop. 
        double[][][] muSums = new double[numConditions][numComp][]; //Results of mu maximization summations for individual components across genome
        int[][] muSumStarts = new int[numConditions][numComp]; //Start positions of muSum arrays (start of maximization window).
        int[][] muSumWidths = new int[numConditions][numComp]; //Effective widths of muSum arrays (width of maximization window).
        int[][][][] newMu = new int[numConditions][numComponents][][];// mu update  
        
        //Initialize arrays
        for(int c=0; c<numConditions; c++)
        	for(int j=0;j<numComp;j++)
        		newMu[c][j] = new int[numBindingType[c]][2];       
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		int numBases = hitNum[c];
    		totalResp[c] = new double[numBases];
    		for(int i=0;i<numBases;i++)
    			totalResp[c][i]=0;
    	}
        //Alpha is annealed in. Alpha=0 during ML steps
        double[] currAlpha = new double[numConditions];
        for(int c=0; c<numConditions; c++)
        	currAlpha[c] = 0;
        
        //Beta is annealed in. Beta=0 during ML steps
        double[][] currBeta = new double[numConditions][numComponents];
        for(int c=0; c<numConditions; c++)
        	for(int j=0;j<numComp;j++)
        		currBeta[c][j] = 0;
        
        //Epsilon is annealed in. Epsilon=0 during ML steps
        double[][][][] currEpsilon = new double[numConditions][numComponents][][];
        for(int c=0; c<numConditions; c++){
        	for(int j=0;j<numComp;j++){
        		currEpsilon[c][j] = new double[numBindingType[c]][2];
        		for (int bt=0; bt< numBindingType[c]; bt++)
        			for (int s=0; s< 2; s++)
        				currEpsilon[c][j][bt][s]=0;
        	}
        }       
        
        
    	////////////
        //Plot the initial pi & priors if plotting
    	////////////
        if(plotEM){
        	String regStr = currRegion.getLocationString().replaceAll(":", "-");
        	String outName = "EM_"+regStr+"_r"+trainingRound+"_t0";
        	int trimLeft = Math.max(0, plotSubRegion.getStart()-regStart);
        	int trimRight = Math.max(0, currRegion.getEnd()-plotSubRegion.getEnd());
        	
//        	if(config.useMotifPrior())
//       			EMStepPlotter.execute(outName, currRegion, mu, pi, forMotifPrior, revMotifPrior, numConditions, numComponents, 0, trimLeft, trimRight);
//       		else
//       			EMStepPlotter.execute(outName, currRegion, mu, pi, null,null, numConditions, numComponents, 0, trimLeft, trimRight);
        }
        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run EM while not converged
        // Note: iterations during which we eliminate a binding component don't count towards "t"
    	//////////////////////////////////////////////////////////////////////////////////////////
        int t=0, iter=0;
        while(t<config.MAX_EM_ITER){ //System.out.println(t); 

    		////////
    		//E-step
    		////////
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];
        		//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
        		for(int i=0;i<numBases;i++)
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				for (int bt=0; bt<numBindingType[c];bt++){
        					for (int s=0; s<2; s++){ if (tau[c][j][bt][s]>0){
        						int index = mu[c][j][bt][s]-regStart;
        						h[c][j][i][bt][s] = hAll[c][index][i][bt][s];
        		}}}}}

        		//Compute responsibilities
        		for(int i=0;i<numBases;i++)
        			totalResp[c][i] = 0;
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				for (int bt=0; bt<numBindingType[c]; bt++){
        					for (int s=0; s<2; s++){ if (tau[c][j][bt][s]>0){
        						// rBind is computed by multiplying h function by pi and pType
        						rBind[c][j][i][bt][s]= h[c][j][i][bt][s]*pi[c][j]*tau[c][j][bt][s];
        						totalResp[c][i] +=rBind[c][j][i][bt][s];
        					}}}}}
        			rNoise[c][i] = n[c][i] * piNoise[c];
        			totalResp[c][i] +=rNoise[c][i];
        		}
        		//Normalize responsibilities
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				for (int bt=0; bt< numBindingType[c];bt++){
        					for (int s=0; s<2 ; s++){ if (tau[c][j][bt][s]>0){
        						rBind[c][j][i][bt][s]/=totalResp[c][i];
        			}}}}}  
        			rNoise[c][i]/=totalResp[c][i];
        		}
    		}  		
        		
    		/////////////////////
    		//M-step: maximize mu (positions)
    		/////////////////////
    		//Set up variable arrays if necessary (assign memory only once to non-zero components)
			if(numConditions>1 && t==config.ALPHA_ANNEALING_ITER)
				for(int c=0; c<numConditions; c++)
					for(int j=0;j<numComp;j++)
						if(pi[c][j]>0)
							muSums[c][j] = new double[config.EM_MU_UPDATE_WIN*2];
			
			//Maximize mu part 1: calculate maximization sums assuming no events shared across conditions
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				double weightedMu = 0;
    				for (int bt=0; bt< numBindingType[c]; bt++)
    					for (int s=0; s<2; s++){if (tau[c][j][bt][s]>0){
    						weightedMu += ((double) mu[c][j][bt][s])*tau[c][j][bt][s];
    					}}
    				int wMu = (int) Math.round(weightedMu);    				
    				int start=Math.max(wMu-config.EM_MU_UPDATE_WIN, regStart);
        			int end = Math.min(currRegion.getEnd(), wMu+config.EM_MU_UPDATE_WIN);
        			//Assign special variables
        			if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER){
        				muSumStarts[c][j] = start; muSumWidths[c][j] = end-start;
        			}
        			//Score the current window
        			double[][] currScore = new double[numBindingType[c]][2]; // scores are indexed by binding type and strand        			
        			double[][] maxScore = new double[numBindingType[c]][2];
        			int[][] maxPos = new int[numBindingType[c]][2];
        			for (int bt=0; bt< numBindingType[c]; bt++){
    					for (int s=0; s<2; s++){
    						currScore[bt][s]=0;
    						maxScore[bt][s]=-Double.MAX_VALUE;
    						maxPos[bt][s]=0;
    					}
        			}
        			
        			for(int x=start; x<end; x++){
        				for (int bt=0; bt<numBindingType[c];bt++)
        					for (int s=0; s<2 ; s++)
        						currScore[bt][s]=0;
        				for(int i=0;i<numBases;i++){
        					int dist = hitPos[c][i]-x;
        					for (int bt=0; bt<numBindingType[c];bt++){
        						if (hitPlusStr[c][i]){
        							// Case 1 : component location is at positive strand
        							currScore[bt][0]+=(rBind[c][j][i][bt][0]*hitCounts[c][i])*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist,true); //Watson
        							// Case 2 : component location is at negative strand
        							currScore[bt][1]+=(rBind[c][j][i][bt][1]*hitCounts[c][i])*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(-dist, false); //Crick
        							
        						}else{
        							// Case 1 : component location is at positive strand
        							currScore[bt][0]+=(rBind[c][j][i][bt][0]*hitCounts[c][i])*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(dist,false);	//Crick
        							// Case 2 : component location is at negative strand
        							currScore[bt][1]+=(rBind[c][j][i][bt][1]*hitCounts[c][i])*TagProbabilityDensities[repIndices[c][i]][bt].logProbability(-dist, true); //Watson
        				}}}		
       				
        				ExperimentCondition ec = manager.getConditions().get(c); 
        				// add scores from motif prior
        				if(bindingManager.getMotifs(ec)!=null && config.useMotifPrior()){
        					for (int bt=0; bt < numBindingType[c]; bt++){
        						int motifIndex = bindingManager.getMotifIndexes(ec).get(bt);
        						if (motifIndex != -1){
        							currScore[bt][0] += (forMotifPrior[c][motifIndex][x-regStart]*config.getPosPriorScaling());
        							currScore[bt][1] += (revMotifPrior[c][motifIndex][x-regStart]*config.getPosPriorScaling());
        						}
        				}}	
        				
        				// Needed in multi-condition
//        				if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER)   //Save the score
//            				muSums[c][j][x-start] = currScore;
        				
        				// Store positions and strands for each binding type
        				for (int bt=0; bt<numBindingType[c]; bt++){
        					for (int s=0; s<2 ; s++){
        						if (currScore[bt][s] > maxScore[bt][s]){
        							maxPos[bt][s]=x;
        							maxScore[bt][s]=currScore[bt][s];
        				}}}
        			}										
     
        			for (int bt=0; bt<numBindingType[c]; bt++)
    					for (int s=0; s<2 ; s++)
    						newMu[c][j][bt][s]=maxPos[bt][s];

//        			System.out.println("printing log likelihood scores");
//        			for (int bt=0; bt<numBindingType[c]; bt++){
//    					for (int s=0; s<2 ; s++){
//    						System.out.print("binding type "+bt+" strand "+s+" : "+maxScore[bt][s]+"\t");			
 //   					}
//        			}
        		}}
    		}
    		
    		//Update mu values, component strand, and binding type probability
    		for(int c=0; c<numConditions; c++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				for (int bt=0; bt< numBindingType[c]; bt++)
    					for (int s=0; s< 2; s++)
    						mu[c][j][bt][s] = newMu[c][j][bt][s];				
    			}}
    		}
    		
    		//Maximize mu part 3: Resolve duplicate positions (combine & delete one copy)
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];	
        		HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>(); //Position to array index map 
        		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        			double weightedMu = 0;
    				for (int bt=0; bt< numBindingType[c]; bt++)
    					for (int s=0; s<2; s++){ if (tau[c][j][bt][s]>0){
    						weightedMu += ((double) mu[c][j][bt][s])*tau[c][j][bt][s];
    					}}
    				int wMu = (int) Math.round(weightedMu);
        			if(pos2index.containsKey(wMu)){ 
        				int orig = pos2index.get(wMu);
        				//Combine
        				pi[c][orig]+=pi[c][j];
                       	for(int i=0; i<numBases;i++){
                       		for (int bt=0; bt<numBindingType[c];bt++)
                       			for (int s=0; s<2; s++)
                       				rBind[c][orig][i][bt][s] += rBind[c][j][i][bt][s];
                       	}
                       	//Delete
                       	pi[c][j]=0.0;
                       	for(int i=0; i<numBases;i++){
                       		for (int bt=0; bt<numBindingType[c];bt++)
                       			for (int s=0; s<2; s++)
                       				rBind[c][j][i][bt][s] = 0;
                       	}
        			}else{
        				pos2index.put(wMu, j);
        			}
        		}}
    		}
        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
    		boolean componentEliminated=false;
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];	
        		//Maximize pi
        		double[] sumR=new double[numComponents];
        		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        			for(int i=0;i<numBases;i++){
        				for (int bt=0; bt< numBindingType[c]; bt++)
        					for (int s=0; s< 2; s++){ if (tau[c][j][bt][s]>0){
        						sumR[j] += rBind[c][j][i][bt][s]*hitCounts[c][i];	
        			}}}
                }}
        		int minIndex=0; double minVal=Double.MAX_VALUE;
        		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        			if(sumR[j]<minVal){ minVal=sumR[j]; minIndex=j;}
        		}}                
                if(minVal>currAlpha[c]){
                    // No component to be eliminated, update pi(j)
                	for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
                		pi[c][j]=Math.max(0, sumR[j]-currAlpha[c]); 
                	}}
                }else{
                    // Eliminate worst binding component
                    // Responsibilities will be redistributed in the E step
                   	pi[c][minIndex]=0.0; sumR[minIndex]=0.0;
                   	for(int i=0; i<numBases;i++){
                   		for (int bt=0; bt<numBindingType[c]; bt++)
                   			for (int s=0; s<2; s++)
                   				rBind[c][minIndex][i][bt][s] = 0;	
                   	}
                   	//I discussed this bit with Chris, and we decided that the best thing to do is
                   	//to re-estimate pi values for non-eliminated components using the current responsibility assignments
                   	for(int j=0;j<numComp;j++){ 
                   		if(j!=minIndex)
                   			pi[c][j]=Math.max(0, sumR[j]); 
                	}
                   	componentEliminated=true;
                }
                //Normalize pi (accounting for piNoise)
                double totalPi=0;
                for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
            		totalPi+=pi[c][j];
            	}}
                for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
            		if(totalPi>0)
            			pi[c][j]=pi[c][j]/(totalPi/(1-piNoise[c]));
            	}}
            	
        		
            	/////////////
            	//Anneal alpha
            	//////////////
        		if (t >config.EM_ML_ITER && t <= config.ALPHA_ANNEALING_ITER)
        			currAlpha[c] = alphaMax[c] * (t-config.EM_ML_ITER)/(config.ALPHA_ANNEALING_ITER-config.EM_ML_ITER);
        		else if(t > config.ALPHA_ANNEALING_ITER)
        			currAlpha[c] = alphaMax[c];
        	}
    		
        	//Non-zero components count
        	int nonZeroComps=0;
        	for(int c=0; c<numConditions; c++)
        		for(int j=0;j<numComp;j++)
        			if(pi[c][j]>0)
        				nonZeroComps++;
    		
    		/////////////////////
    		//M-step: maximize tau
    		/////////////////////
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];	   		
    			ExperimentCondition ec = manager.getIndexedCondition(c);
    			
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				double[][] compSumR = new double[numBindingType[c]][2];
    				for (int bt=0; bt< numBindingType[c]; bt++)
						for (int s=0; s< 2; s++)
							if (tau[c][j][bt][s]>0){ compSumR[bt][s] = 0;}
    				for(int i=0;i<numBases;i++)
    					for (int bt=0; bt< numBindingType[c]; bt++)
    						for (int s=0; s< 2; s++)
    							if (tau[c][j][bt][s]>0){ compSumR[bt][s] += rBind[c][j][i][bt][s]*hitCounts[c][i];;}
    				int minBtIndex=0; int minSIndex=0; double minVal=Double.MAX_VALUE;
    				for (int bt=0; bt< numBindingType[c]; bt++){
						for (int s=0; s< 2; s++){ if (tau[c][j][bt][s]>0){
							if((compSumR[bt][s]+currEpsilon[c][j][bt][s])<minVal){ minVal=(compSumR[bt][s]+currEpsilon[c][j][bt][s]) ;minBtIndex=bt;minSIndex=s;}
					}}}
    				if(minVal>currBeta[c][j]){
    					// No component to be eliminated, update tau
    					for (int bt=0; bt< numBindingType[c]; bt++){
    						for (int s=0; s< 2; s++){ if (tau[c][j][bt][s]>0){ 
    								tau[c][j][bt][s]=Math.max(0, compSumR[bt][s]-currBeta[c][j]+currEpsilon[c][j][bt][s]);
    						}}}
    				}else{
    					// Eliminate worst binding component subtype
                        // Responsibilities will be redistributed in the E step
    					tau[c][j][minBtIndex][minSIndex]=0.0; compSumR[minBtIndex][minSIndex]=0.0;
    					for(int i=0; i<numBases;i++)
    						rBind[c][j][i][minBtIndex][minSIndex]=0;
    					//to re-estimate tau values for non-eliminated components using the current responsibility assignments
    					for (int bt=0; bt< numBindingType[c]; bt++)
    						for (int s=0; s< 2; s++)
    							if ( bt!=minBtIndex || s!= minSIndex)
    								tau[c][j][bt][s]=Math.max(0, compSumR[bt][s]);
    				}    				
    				// Normalize tau
    				double totalTau = 0;
    				for (int bt=0; bt< numBindingType[c]; bt++)
						for (int s=0; s< 2; s++)
							if (tau[c][j][bt][s]>0){totalTau += tau[c][j][bt][s];}
    				for (int bt=0; bt< numBindingType[c]; bt++){
						for (int s=0; s< 2; s++){ if (tau[c][j][bt][s]>0){
							if (totalTau>0)
								tau[c][j][bt][s]/=totalTau;
					}}}		
    				
    				//Update betaMax
    				double sumR=0.0;
    				for (int bt=0; bt< numBindingType[c]; bt++)
						for (int s=0; s< 2; s++)
							sumR+=compSumR[bt][s];
    				betaMax[c][j]=sumR*config.getBetaScalingFactor();
    				   				
    				/////////////
                	//Anneal beta
                	//////////////
            		if (t >config.EM_ML_ITER && t <= config.ALPHA_ANNEALING_ITER)
            			currBeta[c][j] = betaMax[c][j] * (t-config.EM_ML_ITER)/(config.ALPHA_ANNEALING_ITER-config.EM_ML_ITER);
            		else if(t > config.ALPHA_ANNEALING_ITER)
            			currBeta[c][j] = betaMax[c][j];
    				
    				//Update epsilonMax
    				if(bindingManager.getMotifs(ec)!=null && config.useMotifPrior()){
    					double[][] motif_w = new double[numBindingType[c]][2];
    					for (int bt=0; bt< numBindingType[c]; bt++){
    						int motifIndex = bindingManager.getMotifIndexes(ec).get(bt);
    						if (motifIndex != -1){
    							int start = Math.max(mu[c][j][bt][0]-config.MOTIF_FINDING_SEQWINDOW/2 , regStart);
    							int end = Math.min(mu[c][j][bt][0]+config.MOTIF_FINDING_SEQWINDOW/2, regEnd);
    							double maxMotifScore=0.0;
    							for (int i=start; i<=end; i++)
    								if (forMotifPrior[c][motifIndex][i-regStart] > maxMotifScore)
    									maxMotifScore=forMotifPrior[c][motifIndex][i-regStart];   						
    							motif_w[bt][0]= maxMotifScore;
						
    							// Evaluate negative strands
    							start = Math.max(mu[c][j][bt][1]-config.MOTIF_FINDING_SEQWINDOW/2 , regStart);
    							end = Math.min(mu[c][j][bt][1]+config.MOTIF_FINDING_SEQWINDOW/2, regEnd);
    							maxMotifScore=0.0;
    							for (int i=start; i<=end; i++)
    								if (revMotifPrior[c][motifIndex][i-regStart] > maxMotifScore)
    									maxMotifScore=revMotifPrior[c][motifIndex][i-regStart];						
    							motif_w[bt][1]= maxMotifScore;
    						}else{
    							motif_w[bt][0]= 0;
    							motif_w[bt][0]= 0;
    						}
    					}   				
    					//Scale according to sum of responsibilities and scaling factor, and normalize epsilon 
    					for (int bt=0; bt< numBindingType[c]; bt++){
    						int motifIndex = bindingManager.getMotifIndexes(ec).get(bt);
    						if (motifIndex != -1){
    							double motifMaxScore = bindingManager.getMotifs(ec).get(motifIndex).getMaxScore();
    							for (int s=0; s< 2; s++)
    								epsilonMax[c][j][bt][s] = motif_w[bt][s]*sumR*config.getEpsilonScalingFactor()/motifMaxScore;
    						}
    					}					
    					
    					/////////////
    					//Anneal Epsilon
    					//////////////   					
    					if (t <= config.ALPHA_ANNEALING_ITER){ // There is checking condition before that t >config.EM_EPSILON_ITER
    						for (int bt=0; bt< numBindingType[c]; bt++)
    							for (int s=0; s< 2; s++){
    								currEpsilon[c][j][bt][s] = epsilonMax[c][j][bt][s] * t/config.ALPHA_ANNEALING_ITER; 
    							}
    					}else if(t > config.ALPHA_ANNEALING_ITER){
    						for (int bt=0; bt< numBindingType[c]; bt++)
    							for (int s=0; s< 2; s++)
    								currEpsilon[c][j][bt][s] = epsilonMax[c][j][bt][s];
    					}
    				}
    			}}
            }
        	
        	////////////
        	//Compute LL
        	////////////
        	LAP=0;
        	if(config.CALC_LL){
	        	//Log-likelihood calculation
	            double LL =0;
	            for(int c=0; c<numConditions; c++){
	        		int numBases = hitNum[c];
	        		for(int i=0;i<numBases;i++){
	        			// for each read, each event will give a conditional prob or bg prob
	                    double j_sum=0;
	        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0.0){
	        				for (int bt=0; bt< numBindingType[c]; bt++)
	        					for (int s=0; s< 2; s++)
	        						j_sum += Math.log(rBind[c][j][i][bt][s])/config.LOG2;
	                    }}
	        			j_sum += Math.log(rNoise[c][i])/config.LOG2;
	                    
	        			LL += j_sum*hitCounts[c][i];                        
	                }
	            }
	            //Log priors
	            double LP=0;
	            for(int c=0; c<numConditions; c++){
	            	ExperimentCondition cond = manager.getIndexedCondition(c);
	            	//sum of pi
	            	double sum_log_pi=0;
	            	for(int j=0;j<numComp;j++){ if(pi[c][j]>0.0){
	            		sum_log_pi+=Math.log(pi[c][j])/config.LOG2;
	            	}}
	            	//Positional priors
	            	double sum_pos_prior=0;
	            	//TODO: how do we account for multi-condition "prior" here?
	            	//for(int x=0; x<currRegion.getWidth(); x++)
	            	//	sum_pos_prior+=;
	            	
	            	//Motif prior
//TODO:	            	//I don't think I'm doing this correctly
	            	if(forMotifPrior!=null && config.useMotifPrior()){
	            		for (int i=0; i < bindingManager.getMotifs(cond).size();i++){
	            			for(int x=0; x<currRegion.getWidth(); x++)
	            				sum_pos_prior += Math.max(forMotifPrior[c][i][x], revMotifPrior[c][i][x]);
	            		}	            		
	            	}
	            	LP+=-(currAlpha[c]*sum_log_pi)+sum_pos_prior;
	            }
	            LAP = LL+LP;
	
	            System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
            }
        	
        	////////////
            //Plot the current pi & priors if plotting
        	////////////
            if(plotEM){
            	String regStr = currRegion.getLocationString().replaceAll(":", "-");
            	String outName = "EM_"+regStr+"_r"+trainingRound+"_t"+(t+1)+"_i"+(iter+1);
            	int trimLeft = Math.max(0, plotSubRegion.getStart()-regStart);
            	int trimRight = Math.max(0, currRegion.getEnd()-plotSubRegion.getEnd());

//           		if(config.useMotifPrior())
//           			EMStepPlotter.execute(outName, currRegion, mu, pi, forMotifPrior, revMotifPrior, numConditions, numComponents, t+1, trimLeft, trimRight);
//           		else
//           			EMStepPlotter.execute(outName, currRegion, mu, pi, null,null, numConditions, numComponents, t+1, trimLeft, trimRight);
            }

            //Is current state equivalent to the last?
            if(((numConditions>1 && t>config.POSPRIOR_ITER) || (numConditions==1 && t>config.ALPHA_ANNEALING_ITER)) && 
            		lastEquivToCurr())
            	stateEquivCount++;
            else
            	stateEquivCount=0;
                        

    		//Tick the clock forward
    		if(!componentEliminated)
    			t++;
    		iter++;
    		
            ////////////
          	//Check Stopping condition
          	////////////
            if (nonZeroComps>0 && ((numConditions>1 && t<=config.POSPRIOR_ITER) || (numConditions==1 && t<=config.ALPHA_ANNEALING_ITER) || (config.CALC_LL && Math.abs(LAP-lastLAP)>config.EM_CONVERGENCE) || stateEquivCount<config.EM_STATE_EQUIV_ROUNDS)){
                copyStateToLast();
                lastLAP = LAP;
                continue;
            }else{
            	copyStateToLast();
            	lastLAP = LAP;
            	//if(config.isVerbose())
            		//System.err.println("\tRegTrain:"+trainingRound+"\t"+currRegion.getLocationString()+"\t"+currRegion.getWidth()+"\t"+t+"\t"+iter+"\t"+nonZeroComps);
            	break;
            }
        } //LOOP: Run EM while not converged
    }//end of EM_MAP method
      
  
    /**
     * Set responsibility profile for each component (for kernel update)
     * @param bindComponents
     * @param signals
     * @param responsibilities
     * @param c2b
     */
    private void setComponentResponsibilityProfiles(List<List<BindingSubComponents>> bindComponents, List<List<StrandedBaseCount>> signals, 
            									double[][][][][] responsibilities) {
		for(ExperimentCondition cond : manager.getConditions()){
			int c = cond.getIndex();			
			for(int j=0;j<bindComponents.get(c).size();j++){
				
				boolean printRegion=false;				
        		Region top1Abf1 = new Region(bindComponents.get(c).get(j).getCoord().getGenome(), "7", 596312, 596372);
        		Region top2Abf1 = new Region(bindComponents.get(c).get(j).getCoord().getGenome(), "12", 332253, 332313); 
        		Region top3Abf1 = new Region(bindComponents.get(c).get(j).getCoord().getGenome(), "6", 5235, 5295); 
//        		if (top1Abf1.contains(bindComponents.get(c).get(j).getCoord()) || top2Abf1.contains(bindComponents.get(c).get(j).getCoord()) || top3Abf1.contains(bindComponents.get(c).get(j).getCoord()))
//        			printRegion = true;	 
				
				BindingSubComponents comp = bindComponents.get(c).get(j);
				int jr = comp.getIndex();		
		    	for(ControlledExperiment rep : cond.getReplicates()){
		    		List<StrandedBaseCount> bases = signals.get(rep.getIndex());		
			    	double[][][][] rc = responsibilities[c];		   
			    	int center = config.MAX_BINDINGMODEL_WIDTH/2;
			   		// store binding profile (read responsibilities in c condition) of this component
					double[][][] sub_profile_plus = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
					double[][][] sub_profile_minus = new double[numBindingType[c]][2][config.MAX_BINDINGMODEL_WIDTH];
					for(int i=0;i<bases.size();i++){
						StrandedBaseCount base = bases.get(i);	
						int offset = base.getCoordinate()-comp.getPosition()+center;
						if(offset>=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							// Do not flip arrays for now to assign responsibilities for each sub component
							for (int bt=0; bt< numBindingType[c]; bt++){
								for (int s=0; s< 2; s++){							
									if (base.getStrand()=='+')
										sub_profile_plus[bt][s][offset]=rc[jr][i][bt][s]*base.getCount();
									else // if components at negative strand
										sub_profile_minus[bt][s][offset]=rc[jr][i][bt][s]*base.getCount();
								}
							}
						}
					}
					// Set subcomponents profiles
					comp.setSubReadProfile(rep.getIndex(), sub_profile_plus, '+');
					comp.setSubReadProfile(rep.getIndex(), sub_profile_minus, '-');
					
					double[] profile_plus = new double[config.MAX_BINDINGMODEL_WIDTH];
					double[] profile_minus = new double[config.MAX_BINDINGMODEL_WIDTH];
					for (int bt=0; bt< numBindingType[c]; bt++){
						for (int s=0; s< 2; s++){
							for (int i=0; i< config.MAX_BINDINGMODEL_WIDTH; i++){
								profile_plus[i]+=sub_profile_plus[bt][s][i];
								profile_minus[i]+=sub_profile_minus[bt][s][i];
							}	
						}
					}
					// Set overall read profiles
					comp.setReadProfile(rep.getIndex(), profile_plus, '+');
					comp.setReadProfile(rep.getIndex(), profile_minus, '-');
					
					if (printRegion){
						System.out.println("====Printing responsibilities ====");
						StrandedPoint spoint = new StrandedPoint(bindComponents.get(c).get(j).getCoord().getGenome(), bindComponents.get(c).get(j).getCoord().getChrom(), bindComponents.get(c).get(j).getPosition(), '+'); 
						System.out.println(bindComponents.get(c).get(j).getCoord().toString());
						System.out.println("#total resp for watson");
						for (int i=0; i < profile_plus.length; i++)
							System.out.print(","+profile_plus[i]);
						System.out.println("#total resp for crick");
						for (int i=0; i < profile_plus.length; i++)
							System.out.print(","+profile_minus[i]);
						for (int bt=0; bt< numBindingType[c]; bt++){
							System.out.println("binding type "+bt);
							System.out.println("#watson");
							for (int i=0; i < sub_profile_plus[bt][0].length; i++)
								System.out.print(","+sub_profile_plus[bt][0][i]);
							System.out.println("#crick");
							for (int i=0; i < sub_profile_minus[bt][1].length; i++)
								System.out.print(","+sub_profile_minus[bt][0][i]);
						}
						System.out.println("printing actual tag landscape");
						List<StrandedPoint> plist = new ArrayList<StrandedPoint>();
		        		plist.add(spoint);
		        		CompositeTagDistribution maker = new CompositeTagDistribution(plist, manager, config.MAX_BINDINGMODEL_WIDTH, true);
		        		double[] watsonTags =maker.getPointWatson(spoint, cond);
		        		double[] crickTags = maker.getPointCrick(spoint, cond);
		        		// reverse watson and crick tags
		        		double[] rWatsonTags = new double[config.MAX_BINDINGMODEL_WIDTH];
		        		double[] rCrickTags = new double[config.MAX_BINDINGMODEL_WIDTH];
		        		for (int i=0; i< config.MAX_BINDINGMODEL_WIDTH; i++){
		        			rWatsonTags[i] = crickTags[config.MAX_BINDINGMODEL_WIDTH-i-1];
		        			rCrickTags[i] = watsonTags[config.MAX_BINDINGMODEL_WIDTH-i-1];
		        		}
		        		if (printRegion){
		        			System.out.println("actual tag landscape");
		        			System.out.println("#watson");
		        			for (int i=0; i< watsonTags.length; i++)
		        				System.out.print(watsonTags[i]+",");
		        			System.out.println("\n"+"#crick");
		        			for (int i=0; i< crickTags.length; i++)
		        				System.out.print(crickTags[i]+",");
		        		}					
					}
		    	}
			}
		}
	}//end of setComponentResponsibilities method
    
    /**
     * Copy current variables to last variables (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     */
    private void copyStateToLast(){
    	int numC = manager.getNumConditions();
    	for(int c=0; c<numC; c++){
    		for(int j=0; j<numComponents; j++){
    			lastPi[c][j] = pi[c][j];
    			for (int bt=0; bt< numBindingType[c];bt++)
    				for (int s=0; s< 2; s++)
    					lastMu[c][j][bt][s] = mu[c][j][bt][s];
    			for(int x=0; x<rBind[c][j].length; x++){
    				lastRBind[c][j][x] = rBind[c][j][x];
    			}
    		}
    	}
    }
    
    /**
     * Compare last variables to current (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     * @return
     */
    private boolean lastEquivToCurr(){
    	int numC = manager.getNumConditions();
    	int currNZ=0, lastNZ=0;
    	for(int c=0; c<numConditions; c++)
    		for(int j=0;j<numComponents;j++){
    			if(pi[c][j]>0)
    				currNZ++;
    			if(lastPi[c][j]>0)
    				lastNZ++;
    		}
    	boolean numCompEqual = currNZ==lastNZ;
    	boolean compPosEqual=true;
    	if(numCompEqual){
    		for(int c=0; c<numC; c++)
    			for(int j=0; j<mu[c].length; j++){if(pi[c][j]>0){
    				for (int bt=0; bt< numBindingType[c]; bt++)
    					for (int s=0; s< 2; s++)
    						compPosEqual = compPosEqual && mu[c][j][bt][s]==lastMu[c][j][bt][s];
    			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
    	for(int c=0; c<numC; c++)
			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<config.EM_STATE_EQUIV_THRES);
			}}
    	boolean rBindEquivalent=true;
    	for(int c=0; c<numC; c++)
    		for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    			for(int x=0; x<rBind[c][j].length; x++){
    				for (int bt=0; bt<numBindingType[c]; bt++)
    					for (int s=0; s<2; s++)
    						rBindEquivalent = rBindEquivalent && (Math.abs(rBind[c][j][x][bt][s]-lastRBind[c][j][x][bt][s])<config.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
}
