package org.seqcode.projects.chexmix.mixturemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingModel;
import org.seqcode.projects.chexmix.framework.XOGPSConfig;



/**
 * BindingEM: run EM training with sparse prior and positional prior(s) on binding data.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class MultiGPSBindingEM {

	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected XOGPSConfig config;
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
	protected double[][]   hitCounts;	// Hit weights
	protected int[][]      hitPos;		// Hit positions
	protected boolean[][]  hitPlusStr;	// Hit positive strand boolean
	protected int[]		   hitNum;		// Number of hits in each condition 
	protected int[][]      repIndices;  // Index of replicate for the hit
	protected double[][][] hAll;		// H function values for all positions in the current window (precomputed)
	protected double[][][] h; 			// H function (binding component probability per read)
	protected double[][]   n; 			// N function (noise component probability per read)
	protected double[][][] rBind;		// Binding component responsibilities
	protected double[][]   rNoise;		// Noise component responsibilities
	protected double[][]   pi;			// pi : emission probabilities for binding components
	protected double[]     piNoise;		// pi : emission probabilities for noise components (fixed)
	protected int[][]      mu;			// mu : positions of the binding components
	protected double []    alphaMax;	// Maximum alpha
	protected double[][]   motifPrior;  // Motif prior (indexed by condition & base) 
	protected BindingModel[] bindingModels; //Array of binding models for convenience
	protected double[][][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[][]   lastPi;		//Last Pi (monitor convergence)
	protected int[][]      lastMu;		//Last positions (monitor convergence)
	protected double lastLAP, LAP; 		//log-likelihood monitoring
	protected boolean plotEM=false;		//Plot the current region components
	protected Region plotSubRegion=null; //Sub region to plot
	protected double initCondPosPriorVar=10, finalCondPosPriorVar=0.001, currCondPosPriorVar=initCondPosPriorVar; //positional prior final Gaussian variance
	protected double numPotentialRegions;
	protected double probAgivenB, probAgivenNOTB;
	protected int stateEquivCount=0;
	
	/**
	 * Constructor
	 * @param c
	 * @param eMan
	 */
	public MultiGPSBindingEM(XOGPSConfig c, ExperimentManager eMan, BindingManager bMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg){
		config=c;
		manager = eMan;
		bindingManager = bMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
		
		//Positional prior constants & initialization for ML steps
        //TODO: The prior could be worked out explicitly for the multicondition case, but that's not a priority since it shouldn't make a huge difference. 
        //As it stands, we are assuming that we know how many total binding sites there are (# potential regions), and that each condition has an equal number of binding sites, 
        //and that those binding sites are shared randomly between conditions with the fixed sharing rate between two conditions.
        double N = numPotentialRegions;
		double S = N*config.getProbSharedBinding();
		double L = (double)config.getGenome().getGenomeLength();
		probAgivenB = Math.log(config.getProbSharedBinding())/Math.log(2);
        probAgivenNOTB =  Math.log((N-S)/(L-N))/Math.log(2);
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
    											  int trainingRound,
    											  Region plotSubRegion){
    	
    	components = comps;
        this.noise = noise;
        numComponents = numComp;
        this.motifPrior = motifPrior;
        this.trainingRound = trainingRound;
    	this.plotSubRegion = plotSubRegion;
        //Matrix initializations
        hitCounts= new double[numConditions][];	// Hit weights
    	hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	repIndices= new int[numConditions][];	    // Index of replicate for the hit
    	hAll = new double[numConditions][][];		// H function values for all positions in the current window (precomputed)
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	n= new double[numConditions][]; 			// N function (noise component probability per read)
    	rBind= new double[numConditions][][];		// Binding component responsibilities
    	rNoise= new double[numConditions][];		// Noise component responsibilities
    	pi = new double[numConditions][numComponents];	// pi : emission probabilities for binding components
    	piNoise = new double[numConditions];		// pi : emission probabilities for noise components (fixed)
    	alphaMax = new double[numConditions];		//Maximum alpha
        mu = new int[numConditions][numComponents];// mu : positions of the binding components
        bindingModels = new BindingModel[manager.getReplicates().size()]; //Array of bindingModels for convenience
        plotEM = (plotSubRegion!=null && plotSubRegion.overlaps(w));
        //Monitor state convergence using the following last variables
        lastRBind = new double[numConditions][][];
        lastPi = new double[numConditions][numComponents];
        lastMu = new int[numConditions][numComponents];
		
        //Initializing data structures
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	
        	//Add bindingModels to array
        	for(ControlledExperiment rep : cond.getReplicates())
        		bindingModels[rep.getIndex()] = bindingManager.getUnstrandedBindingModel(rep);
        	
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
            for(int j=0;j<numComp;j++)
            	mu[c][j] = components.get(c).get(j).getPosition();
            
            //Initialize H function for all positions in the current region
            double[][] hAllc = new double[w.getWidth()][numBases];
            for(int i=0;i<numBases;i++)
            	for(int b=0;b<w.getWidth();b++){
            		int pos = b+w.getStart();
                    int dist = hitPlusStr[c][i] ? hitPos[c][i]-pos: pos-hitPos[c][i];
                    //Standard ChIP-seq / ChIP-exo
                    hAllc[b][i] = bindingModels[repIndices[c][i]].probability(dist);
            	}
    		hAll[c] = hAllc;
    		
            //Initialize responsibility functions
            double[][] hc= new double[numComp][numBases];
            double[] nc = new double[numBases];
            for(int i=0;i<numBases;i++){
            	for(int j=0;j<numComp;j++){
            		int index = mu[c][j]-w.getStart();
                    hc[j][i] = hAll[c][index][i];
                }
            	nc[i] = noise.get(c).scorePosition(hitPos[c][i], repIndices[c][i]);
            }
            h[c] = hc;
            n[c] = nc;
            
            rBind[c] = new double[numComp][numBases];
    		rNoise[c]= new double[numBases];
    		lastRBind[c] = new double[numComp][numBases];
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
	            comp.setPosition(mu[c][j]);
	            double sum_resp = 0.0;	
                for(int i=0;i<hitNum[c];i++){
                    sum_resp += hitCounts[c][i]*rBind[c][j][i];
                }
	            comp.setUnstrandedSumResponsibilities(sum_resp);
	            if(pi[c][j]>0.0){
	            	currActiveComps.add(comp);
	            }
	    	}
	    	activeComponents.add(currActiveComps);
	    	//Noise Components
	    	double noise_resp = 0.0;	
            for(int i=0;i<hitNum[c];i++)
                noise_resp += hitCounts[c][i]*rNoise[c][i];
	    	noise.get(c).setSumResponsibility(noise_resp);
        }        
        setComponentResponsibilityProfiles(activeComponents, signals, rBind);
        
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
        
        //Variables for tracking mu maximization. Defined early to avoid memory assignment during main EM loop. 
        double[][][] muSums = new double[numConditions][numComp][]; //Results of mu maximization summations for individual components across genome
        int[][] muSumStarts = new int[numConditions][numComp]; //Start positions of muSum arrays (start of maximization window).
        int[][] muSumWidths = new int[numConditions][numComp]; //Effective widths of muSum arrays (width of maximization window).
        int[][] muSumMaxPos = new int[numConditions][numComp]; //Positions of maxima in mu maximization summations
        int[] muJoinClosestComps = new int[numConditions]; //Indices of nearest components in other conditions
        boolean[] muJoinSharedBetter = new boolean[numConditions]; //Indicator that sharing components across conditions is better than not
        int[][] newMu = new int[numConditions][numComponents];// mu update
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		int numBases = hitNum[c];
    		totalResp[c] = new double[numBases];
            for(int i=0;i<numBases;i++)
                totalResp[c][i] = 0;
    	}
        //Alpha is annealed in. Alpha=0 during ML steps
        double[] currAlpha = new double[numConditions];
        for(int c=0; c<numConditions; c++)
        	currAlpha[c] = 0;
        
        /**
    	////////////
        //Plot the initial pi & priors if plotting
    	////////////
        if(plotEM){
        	String regStr = currRegion.getLocationString().replaceAll(":", "-");
        	String outName = "EM_"+regStr+"_r"+trainingRound+"_t0";
        	int trimLeft = Math.max(0, plotSubRegion.getStart()-regStart);
        	int trimRight = Math.max(0, currRegion.getEnd()-plotSubRegion.getEnd());
        	
        	if(config.useMotifPrior())
       			EMStepPlotter.execute(outName, currRegion, mu, pi, motifPrior, numConditions, numComponents, 0, trimLeft, trimRight);
       		else
       			EMStepPlotter.execute(outName, currRegion, mu, pi, null, numConditions, numComponents, 0, trimLeft, trimRight);
        }
        **/
        
    	
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
                		int index = mu[c][j]-regStart;
                        h[c][j][i] = hAll[c][index][i];
                	}}
        		//Compute responsibilities
    			for(int i=0;i<numBases;i++)
                    totalResp[c][i] = 0;
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				rBind[c][j][i] = h[c][j][i]*pi[c][j];
        				totalResp[c][i] +=rBind[c][j][i]; 
        			}}
        			rNoise[c][i] = n[c][i] * piNoise[c];
        			totalResp[c][i] +=rNoise[c][i];
        		}
        		//Normalize responsibilities
        		for(int i=0;i<numBases;i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				rBind[c][j][i]/=totalResp[c][i];
        			}}
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
    				int start=Math.max(mu[c][j]-config.EM_MU_UPDATE_WIN, regStart);
        			int end = Math.min(currRegion.getEnd(), mu[c][j]+config.EM_MU_UPDATE_WIN);
        			//Assign special variables
        			if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER){
        				muSumStarts[c][j] = start; muSumWidths[c][j] = end-start;
        			}
        			//Score the current window
        			double currScore=0, maxScore=-Double.MAX_VALUE;
        			int maxPos = 0;
        			for(int x=start; x<end; x++){
        				currScore=0;
        				for(int i=0;i<numBases;i++){
        					int dist = hitPlusStr[c][i] ? hitPos[c][i]-x: x-hitPos[c][i];
                            //Standard ChIP-seq / ChIP-exo
                            currScore+=(rBind[c][j][i]*hitCounts[c][i]) * bindingModels[repIndices[c][i]].logProbability(dist);
        				}
        				if(motifPrior!=null && config.useMotifPrior())
        					currScore += motifPrior[c][x-regStart];
        				
        				if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER)   //Save the score
            				muSums[c][j][x-start] = currScore;
        				 
        				if(currScore>maxScore){
        					maxPos=x;
        					maxScore=currScore;
        				}
        			}
        			muSumMaxPos[c][j] = maxPos; 
        		}}
    		}
    		//Maximize mu part 2: evaluate whether joining nearby components across conditions is more favorable 
    		for(int c=0; c<numConditions; c++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER && config.useMultiConditionPosPrior()){
    					//mu2.a: find the closest components to j in each condition
    					int closestComp=-1; int closestDist = Integer.MAX_VALUE;
    					for(int d=0; d<numConditions; d++){ if(d!=c){
    						closestComp=-1; closestDist = Integer.MAX_VALUE;
    		    			for(int k=0;k<numComp;k++){ if(pi[d][k]>0){
    		    				int dist = Math.abs(mu[c][j]-mu[d][k]);
    		    				if(dist<closestDist && dist<config.EM_MU_UPDATE_WIN){
    		    					closestDist = dist; closestComp=k;
    		    				}
    		    			}
    		    			muJoinClosestComps[d]=closestComp;
    		    		}}}
    					//mu2.b: evaluate each pair of conditions, asking if a shared event involving j and its closest component would be better than independent events
    					int maxMuStart=muSumStarts[c][j];
    					int minMuEnd=muSumStarts[c][j]+muSumWidths[c][j];
    					int numSharedBetter=0;
    					for(int d=0; d<numConditions; d++){ if(d!=c){
    						int k = muJoinClosestComps[d];
    						if(k==-1)
    							muJoinSharedBetter[d]=false;
    						else{
    							//Case 1: two independent components
    							double indepScore = muSums[c][j][muSumMaxPos[c][j]-muSumStarts[c][j]] + probAgivenNOTB +
    												muSums[d][k][muSumMaxPos[d][k]-muSumStarts[d][k]] + probAgivenNOTB;
    							//Case 2: single shared components
    							double maxSharedScore=-Double.MAX_VALUE; int maxSharedPos = 0; double currScore=0;
    		        			for(int y=muSumStarts[c][j]; y<muSumStarts[c][j]+muSumWidths[c][j]; y++){
    		        				if(y>=muSumStarts[d][k] && y<muSumStarts[d][k]+muSumWidths[d][k]){
    		        					currScore = muSums[c][j][y-muSumStarts[c][j]] + probAgivenB +
    		        								muSums[d][k][y-muSumStarts[d][k]] + probAgivenB;
    		        					if(currScore > maxSharedScore){
    		        						maxSharedScore = currScore; maxSharedPos = y;
    		        					}
    		        			}}
    							muJoinSharedBetter[d] = maxSharedScore>indepScore ? true : false;
    							maxMuStart = muJoinSharedBetter[d] ? Math.max(maxMuStart, muSumStarts[d][k]) : maxMuStart; 
    							minMuEnd = muJoinSharedBetter[d] ? Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]) : minMuEnd;
    							if(muJoinSharedBetter[d])
    								numSharedBetter++;
    							
    							//mu2.d: update mu (Shortcut for numConditions==2)
    							if(numConditions == 2){
    								if(muJoinSharedBetter[d]) 
    									newMu[c][j] = maxSharedPos;
    								else
    									newMu[c][j] = muSumMaxPos[c][j];
    							}
    						}
    					}}
    					
    					//mu2.c: for all conditions that passed the pairwise test, evaluate if a single shared event is better than all independent
    					if(numConditions>2){  //Shortcut for 2 conditions above
	    					//Case 1: sum of all independent components
	    					double allIndepScore = muSums[c][j][muSumMaxPos[c][j]-muSumStarts[c][j]] + probAgivenNOTB;
	    					for(int d=0; d<numConditions; d++){ if(d!=c){
	    						int k = muJoinClosestComps[d];
	    						if(k!=-1){
	    							allIndepScore+=muSums[d][k][muSumMaxPos[d][k]-muSumStarts[d][k]] + probAgivenNOTB;
	    						}
	    					}}
	    					//Case 2: sum of shared component and non-shared
	    					double maxSomeSharedScore=-Double.MAX_VALUE; int maxSomeSharedPos = 0; double currScore=0;
	    					for(int y = maxMuStart; y<minMuEnd; y++){ 
	    						currScore=muSums[c][j][y-muSumStarts[c][j]] + probAgivenB;
	    						for(int d=0; d<numConditions; d++){ if(d!=c){
		    						int k = muJoinClosestComps[d];
		    						if(k!=-1){
		    							if(muJoinSharedBetter[d])
		    								currScore += muSums[d][k][y-muSumStarts[d][k]] + probAgivenB;
		    							else
		    								currScore +=muSums[d][k][muSumMaxPos[d][k]-muSumStarts[d][k]] + probAgivenNOTB;
		    						}
	    						}}
	    						if(currScore > maxSomeSharedScore){
	        						maxSomeSharedScore = currScore; maxSomeSharedPos = y;
	        					}
	    					}
	    					
	    					double maxAllSharedScore=-Double.MAX_VALUE; int maxAllSharedPos = 0; currScore=0;
	    					if(numSharedBetter==numConditions-1){
	    						maxAllSharedScore=maxSomeSharedScore;
	    						maxAllSharedPos = maxSomeSharedPos;
	    					}else{
		    					//mu2.d: Case 3: single shared position, regardless of what happened in the pairwise tests
		    					for(int d=0; d<numConditions; d++){ if(c!=d){//Update window
			    						int k = muJoinClosestComps[d];
			    						if(k!=-1){
			    							maxMuStart = Math.max(maxMuStart, muSumStarts[d][k]); 
			    							minMuEnd = Math.min(minMuEnd, muSumStarts[d][k]+muSumWidths[d][k]);
			    						}
		    					}}
		    					for(int y = maxMuStart; y<minMuEnd; y++){ 
		    						currScore=muSums[c][j][y-muSumStarts[c][j]] + probAgivenB;
		    						for(int d=0; d<numConditions; d++){ if(d!=c){
			    						int k = muJoinClosestComps[d];
			    						if(k!=-1){
			    							currScore += muSums[d][k][y-muSumStarts[d][k]] + probAgivenB;
			    						}
		    						}}
		    						if(currScore > maxAllSharedScore){
		        						maxAllSharedScore = currScore; maxAllSharedPos = y;
		        					}
		    					}
	    					}
	    					
	    					//mu2.e: update mu
	    					if(maxAllSharedScore >=allIndepScore && maxAllSharedScore >=maxSomeSharedScore)
	    						newMu[c][j] = maxAllSharedPos;
	    					else if(maxSomeSharedScore >=allIndepScore)
	    						newMu[c][j] = maxSomeSharedPos;
	    					else
	    						newMu[c][j] = muSumMaxPos[c][j];
	    					
	    				}

    				}else{
    					//Ignore other conditions in first phases of training (until many components are eliminated)
    					newMu[c][j] = muSumMaxPos[c][j];
    				}
    			}}
    		}//Update mu values
    		for(int c=0; c<numConditions; c++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				mu[c][j] = newMu[c][j];
    			}}
    		}
    		//Maximize mu part 3: Resolve duplicate positions (combine & delete one copy)
    		for(int c=0; c<numConditions; c++){ int numBases = hitNum[c];	
        		HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>(); //Position to array index map 
        		for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        			if(pos2index.containsKey(mu[c][j])){ 
        				int orig = pos2index.get(mu[c][j]);
        				//Combine
        				pi[c][orig]+=pi[c][j];
                       	for(int i=0; i<numBases;i++)
                       		rBind[c][orig][i] += rBind[c][j][i];
                       	//Delete
                       	pi[c][j]=0.0;
                       	for(int i=0; i<numBases;i++)
                       		rBind[c][j][i] = 0;
        			}else{
        				pos2index.put(mu[c][j], j);
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
        			for(int i=0;i<numBases;i++)
        				sumR[j] += rBind[c][j][i]*hitCounts[c][i];
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
                   	for(int i=0; i<numBases;i++)
                   		rBind[c][minIndex][i] = 0;
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
	        				j_sum += Math.log(rBind[c][j][i])/config.LOG2;
	                    }}
	        			j_sum += Math.log(rNoise[c][i])/config.LOG2;
	                    
	        			LL += j_sum*hitCounts[c][i];                        
	                }
	            }
	            //Log priors
	            double LP=0;
	            for(int c=0; c<numConditions; c++){
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
	            	if(motifPrior!=null && config.useMotifPrior())
	            		for(int x=0; x<currRegion.getWidth(); x++)
	            			sum_pos_prior += motifPrior[c][x];
	            	
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
//           			EMStepPlotter.execute(outName, currRegion, mu, pi, motifPrior, numConditions, numComponents, t+1, trimLeft, trimRight);
//           		else
//           			EMStepPlotter.execute(outName, currRegion, mu, pi, null, numConditions, numComponents, t+1, trimLeft, trimRight);
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
            									double[][][] responsibilities) {
		for(ExperimentCondition cond : manager.getConditions()){
			int c = cond.getIndex();
			
			for(int j=0;j<bindComponents.get(c).size();j++){
				BindingSubComponents comp = bindComponents.get(c).get(j);
				int jr = comp.getIndex();
			
		    	for(ControlledExperiment rep : cond.getReplicates()){
		    		List<StrandedBaseCount> bases = signals.get(rep.getIndex());
		
			    	double[][] rc = responsibilities[c];
			   
			    	int center = config.MAX_BINDINGMODEL_WIDTH/2;
			   		// store binding profile (read responsibilities in c condition) of this component
					double[] profile_plus = new double[config.MAX_BINDINGMODEL_WIDTH];
					double[] profile_minus = new double[config.MAX_BINDINGMODEL_WIDTH];
					for(int i=0;i<bases.size();i++){
						StrandedBaseCount base = bases.get(i);
						/**
						if (base.getStrand()=='+'){
							int offset = base.getCoordinate()-comp.getPosition()+center;
							if(offset>=0 && offset<config.MAX_BINDINGMODEL_WIDTH)
								profile_plus[offset]=rc[jr][i]*base.getCount();
						}else{
							int offset = comp.getPosition()-base.getCoordinate()+center;
							if(offset>=0 && offset<config.MAX_BINDINGMODEL_WIDTH)
								profile_minus[offset]=rc[jr][i]*base.getCount();
						}
						**/
						int offset = base.getCoordinate()-comp.getPosition()+center;
						if(offset>=0 && offset<config.MAX_BINDINGMODEL_WIDTH){
							if (base.getStrand()=='+')
								profile_plus[offset]=rc[jr][i]*base.getCount();
							else
								profile_minus[offset]=rc[jr][i]*base.getCount();
						}
					}
					comp.setReadProfile(rep.getIndex(), profile_plus,  '+');
					comp.setReadProfile(rep.getIndex(), profile_minus, '-');
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
    			lastMu[c][j] = mu[c][j];
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
    				compPosEqual = compPosEqual && (mu[c][j] == lastMu[c][j]);
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
    				rBindEquivalent = rBindEquivalent && (Math.abs(rBind[c][j][x]-lastRBind[c][j][x])<config.EM_STATE_EQUIV_THRES);
    			}
			}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
}
