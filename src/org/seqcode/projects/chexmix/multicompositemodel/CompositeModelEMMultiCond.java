package org.seqcode.projects.chexmix.multicompositemodel;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.imageio.ImageIO;
import javax.imageio.stream.FileImageOutputStream;
import javax.imageio.stream.ImageOutputStream;

import org.apache.commons.math3.fitting.GaussianFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.projects.chexmix.composite.CompositeModelComponent;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.utilities.TrainingStepPlotter;
import org.seqcode.viz.utils.GifSequenceWriter;


/**
 * CompositeModelEM: run EM training with sparse prior on a composite tag distribution.
 * 
 * This training method do not assume that the same protein-DNA interaction model is valid (i.e. constant) across all examined conditions.
 * Supports joint inference of positions across-conditions.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class CompositeModelEMMultiCond {

	protected ExperimentManager manager;
	protected ChExMixConfig config;
	protected CompositeTagDistribution composite;
	protected ProteinDNAInteractionModelMultiCond model;
	protected int numComponents;  //The count of all components (active +inactive) in the model
	protected int numConditions;
	protected int trainingRound=0; //Identifier for the overall training round, used only for file names
	//	EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[][]   hitCounts;	// Hit weights
	protected int[][]      hitPos;		// Hit positions
	protected boolean[][]  hitPlusStr;	// Hit positive strand boolean
	protected int[]		   hitNum;		// Number of hits in each condition 
	protected double[][][] hAll;		// H function values for all positions in the current window (precomputed)
	protected double[][][] h; 			// H function (binding component probability per read)
	protected double[][][] r;		// Binding component responsibilities
	protected double[][]   pi;			// pi : emission probabilities for binding components
	protected int[][]      mu;			// mu : positions of the binding components
	protected double[]    alphaMax;	// Maximum alpha
	protected double[][][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[][]   lastPi;		//Last Pi (monitor convergence)
	protected int[][]      lastMu;		//Last positions (monitor convergence)
	protected double lastLAP, LAP; 		//log-likelihood monitoring
	protected boolean plotEM=false;		//Plot the current region components
	protected TrainingStepPlotter plotter = null;
	protected List<BufferedImage> fullImages, zoomImages;
	protected double probAgivenB, probAgivenNOTB;
	protected int stateEquivCount=0;
	
	/**
	 * Constructor
	 * 
	 * @param c: XLAnalysisConfig
	 * @param eMan: ExperimentManager
	 * @throws Exception 
	 */
	public CompositeModelEMMultiCond(CompositeTagDistribution composite, ChExMixConfig c, ExperimentManager eMan){
		this.composite=composite;
		config=c;
		manager = eMan;
		numConditions = manager.getNumConditions();		
		//Plotter?
    	plotEM = config.getPlotEM();
        if(plotEM){
        	plotter = new TrainingStepPlotter();
        	fullImages = new ArrayList<BufferedImage>();
        	zoomImages = new ArrayList<BufferedImage>();
        }
        
        //Positional prior constants & initialization for ML steps
        //TODO: The prior could be worked out explicitly for the multicondition case, but that's not a priority since it shouldn't make a huge difference. 
        //As it stands, we are assuming that we know how many total binding sites there are (# potential regions), and that each condition has an equal number of binding sites, 
        //and that those binding sites are shared randomly between conditions with the fixed sharing rate between two conditions.
        double N = config.getN();
		double S = N*config.getProbSharedBinding();
		double L = (double)config.getCompositeWinSize();
		probAgivenB = Math.log(config.getProbSharedBinding())/Math.log(2);
		probAgivenNOTB =  Math.log((N-S)/(L-N))/Math.log(2);
        
        System.out.println("probAgivenB: "+probAgivenB+ "\tprobAgivenNOTB: "+probAgivenNOTB);
        
	}
	
	
	/**
     * EM training
     *
     * Almost purely matrix/array operations.
     * 
     * Returns an updated protein-DNA interaction model
     *
     * @param model: initial model
	 * @param trainingRound: training round number
	 * @param runEM: if false, just assign responsibilities instead of EM (used to get composite-level responsibilities using a trained model) 
	 * @return
	 * @throws Exception
	 */
    public ProteinDNAInteractionModelMultiCond  train(ProteinDNAInteractionModelMultiCond model,
    											  int trainingRound,
    											  boolean runEM
    											  ) throws Exception{
    	this.model=model;
    	//Need to ensure that the composite and model have the same center offsets, or the coordinate system will not be consistent
    	if(composite.getWinSize()!=model.getWidth() || composite.getCenterOffset()!=model.getCenterOffset())
    		throw new Exception("CompositeModelEM: Composite distribution coordinate system not consistent with protein-DNA interaction model");
    			
    	numComponents = model.getNumComponents();
        this.trainingRound = trainingRound;
        //Matrix initializations
        hitCounts= new double[numConditions][];		// Hit weights
    	hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	r= new double[numConditions][][];			// Binding component responsibilities
    	pi = new double[numConditions][numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numConditions][numComponents]; // mu : positions of the binding components (fixed across conditions)
    	alphaMax = new double[numConditions];		// Maximum alpha
        //Monitor state convergence using the following last variables
        lastRBind = new double[numConditions][][];
        lastPi = new double[numConditions][numComponents];
        lastMu = new int[numConditions][numComponents];

    	//Initializing data structures
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	
        	for(int j=0;j<numComponents;j++){
            	//Load pi for binding components
            	CompositeModelComponent comp = model.getAllComponents(c).get(j);
                pi[c][j]= comp.getPi();
                //Load binding component positions
            	mu[c][j] = model.getAllComponents(c).get(j).getPosition();
            }
        	
        	//Set maximum alpha
        	alphaMax[c] =  config.getFixedAlpha()>0 ? config.getFixedAlpha() :
        			Math.max(config.MIN_ALPHA, (config.getAlphaScalingFactor() * model.getBackgroundComponent(c).getPi())/composite.getWinSize()); 
        	System.out.println("\n\tT="+trainingRound+", Condition="+cond.getName()+", Alpha= "+alphaMax[c]);
        	
        	//Number of unique positions in the composite distribution (watson + crick)
        	hitNum[c]=composite.getWinSize()*2;

            //Load read info
            double[] countc= new double[hitNum[c]];
            int[] posc= new int[hitNum[c]];
            boolean[] plusc= new boolean[hitNum[c]];
            for(int i=0;i<composite.getWinSize();i++){ //Watson
            	posc[i] = i;
            	plusc[i] = true;
                countc[i]=composite.getCompositeWatson(cond)[i];
            }
            for(int i=0;i<composite.getWinSize();i++){ //Crick
            	posc[i+composite.getWinSize()] = i;
            	plusc[i+composite.getWinSize()] = false;
                countc[i+composite.getWinSize()]=composite.getCompositeCrick(cond)[i];
            }
            hitPos[c] = posc;
            hitCounts[c]=countc;
            hitPlusStr[c] = plusc;
        	
            //Initialize responsibility functions
            double[][] hc= new double[numComponents][hitNum[c]];
            for(int i=0;i<hitNum[c];i++){
            	for(int j=0;j<numComponents;j++){
            		int dist = hitPos[c][i]-mu[c][j];
                    hc[j][i] = model.getAllComponents(c).get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
                }
            }
            h[c] = hc;
            
            r[c] = new double[numComponents][hitNum[c]];
    		lastRBind[c] = new double[numComponents][hitNum[c]];
    		        	
        }
        //End of data structure initialization
        
        
        //////////
        // Run EM steps
        //////////
        if(runEM)
        	EM_MAP();
        else
        	responsibilityAssignment();
	
        //////////
        // re-assign EM result back to component objects
        //////////

        for(int c=0; c<numConditions;c++){
        	//Binding Components
        	for(int j=0;j<numComponents;j++){ 
        		CompositeModelComponent comp = model.getAllComponents(c).get(j);
        		comp.setPi(pi[c][j]);
        		comp.setPosition(mu[c][j]);
        		double sumRespW=0.0, sumRespC=0.0;	
        		for(int i=0;i<hitNum[c];i++){
        			if(hitPlusStr[c][i])
        				sumRespW += hitCounts[c][i]*r[c][j][i];
        			else
        				sumRespC += hitCounts[c][i]*r[c][j][i];
        		}
        		comp.setSumResponsibilities(sumRespW, sumRespC);
        	}
        	//Responsibility profiles
        	setComponentResponsibilityProfiles(r);
        
        	//Print the responsibilities to files
        	if(config.getPrintCompositeResponsibilities()){
        		for(ExperimentCondition cond : manager.getConditions()){
        			String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()
        					+"_responsibilities."+cond.getName()+"T"+trainingRound+".txt";
        			printResponsibilitiesToFile(cond, filename);
        		}
        	}
        }
        	
        return model;
    }//end of EMTrain method


    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void EM_MAP () {
        int numComp = numComponents;
        double [][] totalResp = new double[numConditions][];
        int[] csCompIndex = new int[numConditions];
        for (int c=0; c< numConditions; c++)
        	csCompIndex[c] = model.getCSComponent(c).getIndex(); 
        int[] bgCompIndex = new int[numConditions];
        for (int c=0; c < numConditions; c++)
        	bgCompIndex[c] = model.getBackgroundComponent(c).getIndex();
        
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
        

    	////////////
        //Plot the initial pi & priors if plotting
    	////////////
        if(plotEM && plotter!=null){
        	for(ExperimentCondition cond : manager.getConditions()){
            	int c = cond.getIndex();
            	String outName = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EM_r"+trainingRound+"_t0_"+cond.getName() : null;
            	String outNameZoom = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EMzoom_r"+trainingRound+"_t0_"+cond.getName() : null;
        		int trimLeft = (model.getWidth()/2)-50;
        		int trimRight = (model.getWidth()/2)-50;
 //       		fullImages.add(plotter.plotCompositeEM(outName, composite, model, mu[c], pi[c], trainingRound, 0, 0,0));
 //       		zoomImages.add(plotter.plotCompositeEM(outNameZoom, composite, model, mu[c], pi[c], trainingRound, 0, trimLeft, trimRight));
        	}
        }
        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run EM while not converged
        // Note: iterations during which we eliminate a binding component don't count towards "t"
    	//////////////////////////////////////////////////////////////////////////////////////////
        int t=0, iter=0;
        while(t<config.MAX_EM_ITER){  
        	
    		////////
    		//E-step
    		////////
    		for(int c=0; c<numConditions; c++){ 
        		//Recompute h function, given binding component positions
        		for(int i=0;i<hitNum[c];i++)
                	for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
                		int dist = hitPos[c][i]-mu[c][j];
                		h[c][j][i] = model.getAllComponents(c).get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
                	}}
        		//Compute responsibilities
    			for(int i=0;i<hitNum[c];i++)
                    totalResp[c][i] = 0;
        		for(int i=0;i<hitNum[c];i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				r[c][j][i] = h[c][j][i]*pi[c][j];
        				totalResp[c][i] +=r[c][j][i]; 
        			}}
        		}
        		//Normalize responsibilities
        		for(int i=0;i<hitNum[c];i++){
        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
        				r[c][j][i]/=totalResp[c][i];
        			}}
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
			//Maximize mu: calculate maximization sums assuming no events shared across conditions
			for(int c=0; c<numConditions; c++){
				for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePosition()){
					//Define window to look for new component position
					int start=Math.max(mu[c][j]-config.EM_MU_UPDATE_WIN, 0);
					int end = Math.min(model.getWidth(), mu[c][j]+config.EM_MU_UPDATE_WIN);
					//Assign special variables
        			if(numConditions>1 && t>config.ALPHA_ANNEALING_ITER){
        				muSumStarts[c][j] = start; muSumWidths[c][j] = end-start;
        			}
					//Score the current window
					double currScore=0, maxScore=-Double.MAX_VALUE;
					int maxPos = 0;
    			
					for(int x=start; x<end; x++){
						currScore=0;
						for(int i=0;i<hitNum[c];i++){
							int dist = hitPos[c][i]-x;
							currScore+=(r[c][j][i]*hitCounts[c][i]) * model.getAllComponents(c).get(j).getTagDistribution().logProbability(dist, hitPlusStr[c][i]);
						}
						
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
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePosition()){
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
    		}		
    		//Update mu values
    		for(int c=0; c<numConditions; c++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePosition()){
    				mu[c][j] = newMu[c][j];
    			}}
    		}

    		//Maximize mu follow-up: Resolve duplicate positions (combine & delete one copy)
    		for(int c=0; c<numConditions; c++){
    			HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>(); //Position to array index map 
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				if(model.getAllComponents(c).get(j).hasUpdatablePosition()){
    					if(pos2index.containsKey(mu[c][j])){ 
    						int orig = pos2index.get(mu[c][j]);
    						//Combine
    						pi[c][orig]+=pi[c][j];
    						for(int i=0; i<hitNum[c];i++)
    							r[c][orig][i] += r[c][j][i];
    						//Delete
    						pi[c][j]=0.0;
    						for(int i=0; i<hitNum[c];i++)
    							r[c][j][i] = 0;
    					}else{
    						pos2index.put(mu[c][j], j);
    					}
    				}
    			}}
    		}    		
        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
    		boolean componentEliminated=false;
    		for(int c=0; c<numConditions; c++){
    			double[] sumR=new double[numComponents];
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				for(int i=0;i<hitNum[c];i++)
    					sumR[j] += r[c][j][i]*hitCounts[c][i];
    			}}
    			int minIndex=0; double minVal=Double.MAX_VALUE;
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				if(sumR[j]<minVal){ minVal=sumR[j]; minIndex=j;}
    			}}                
    			if(minVal>currAlpha[c]){
    				// No component to be eliminated, update pi(j)
    				for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePi()){
    					pi[c][j]=Math.max(0, sumR[j]-currAlpha[c]); 
    				}}
    			}else{
    				// Eliminate worst binding component
    				// Responsibilities will be redistributed in the E step
    				pi[c][minIndex]=0.0; sumR[minIndex]=0.0;
    				for(int i=0; i<hitNum[c];i++)
    					r[c][minIndex][i] = 0;
    				//I discussed this bit with Chris, and we decided that the best thing to do is
    				//to re-estimate pi values for non-eliminated components using the current responsibility assignments
    				for(int j=0;j<numComp;j++){ 
    					if(j!=minIndex)
    						pi[c][j]=Math.max(0, sumR[j]); 
    				}
    				componentEliminated=true;
    			}
    			//Normalize pi
    			double totalRPi=0, totalNonUpdateablePi=0;
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && !model.getAllComponents(c).get(j).hasUpdatablePi()){
    				totalNonUpdateablePi+=pi[c][j];
    			}}
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePi()){
    				totalRPi+=pi[c][j];
    			}}
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePi()){
    				if(totalRPi>0){
    					pi[c][j]/=totalRPi;
    					pi[c][j]*=1-totalNonUpdateablePi;
    				}
    			}}
    			
    			//Special case for CS component:
    			if(pi[c][csCompIndex[c]]<config.MIN_CS_PI){
    				pi[c][csCompIndex[c]]=config.MIN_CS_PI;
    				//Renormalize
    				double totPi=0;
    				for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && j!=csCompIndex[c] && model.getAllComponents(c).get(j).hasUpdatablePi()){ 
    					totPi +=  pi[c][j];
    				}}
    				for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && j!=csCompIndex[c] && model.getAllComponents(c).get(j).hasUpdatablePi()){
	        			pi[c][j]/=totPi;
	        			pi[c][j]*=1-totalNonUpdateablePi-config.MIN_CS_PI;
    				}}
    			}
        	
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
        			if (pi[c][j]>0 && j!=csCompIndex[c] && j!=bgCompIndex[c])
        				nonZeroComps++;
        	
        	////////////
        	//Compute LL
        	////////////
        	LAP=0;
        	if(config.CALC_LL){
	        	//Log-likelihood calculation
	            double LL =0;
	            for(int c=0; c<numConditions; c++){
	        		for(int i=0;i<hitNum[c];i++){
	        			// for each read, each event will give a conditional prob or bg prob
	                    double j_sum=0;
	        			for(int j=0;j<numComp;j++){ if(pi[c][j]>0.0){
	        				j_sum += Math.log(r[c][j][i])/config.LOG2;
	                    }}
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
	            	
	            	LP+=-(currAlpha[c]*sum_log_pi);
	            }
	            LAP = LL+LP;
	
	            System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
            }
        	
        	////////////
        	//Update tag distributions
        	////////////
        	setComponentResponsibilityProfiles(r);
        	
        	for(ExperimentCondition cond : manager.getConditions()){
        		int c= cond.getIndex();
        		//CS component (empirical or gaussian fit)
        		TagProbabilityDensity CSdistrib = model.getCSTagDistribution(c);
        		double[] csW = model.getCSComponent(c).getTagProfile(true);
        		double[] csC = model.getCSComponent(c).getTagProfile(false);
        		//smooth with arbitrary gaussian
        		csW = StatUtil.gaussianSmoother(csW, 5);
        		csC = StatUtil.gaussianSmoother(csC, 5);

        		GaussianFitter fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
        		for(int i=0; i<CSdistrib.getWinSize(); i++){
        			fitter.addObservedPoint((double)(i+CSdistrib.getLeft()), csW[i]+csC[CSdistrib.getWinSize()-i-1]);
        		}
        		double[] parameters = fitter.fit();;
        		double newOffset = parameters[1];
        		double newSigma = parameters[2];
        		CSdistrib.loadGaussianDistrib(newOffset, newSigma, -newOffset, newSigma); //Symmetric CS component
        		
        		//XL components (fit gaussian)
        		for(CompositeModelComponent xlComp : model.getXLComponents(c)){
        			if(xlComp.isNonZero()){
        				TagProbabilityDensity XLdistrib = model.getXLTagDistribution(c, xlComp);

        				double[] xlW = xlComp.getTagProfile(true);
        				double[] xlC = xlComp.getTagProfile(false);
    	    		
        				if(config.fixedXLOffset()){
        					//Mean stays constant - calculate sigmas 
        					double newSigmaW=config.getXLDistribSigma(), newSigmaC=config.getXLDistribSigma(),
        							offsetW = XLdistrib.getGaussOffsetW(), offsetC = XLdistrib.getGaussOffsetC(),
        							sumDiffW=0, sumDiffC=0, totW=0, totC=0;
        					int l=XLdistrib.getLeft();

        					// Trying none symmetric sigma
//        					if(config.XL_DISTRIB_SYMMETRIC){
//        						for(int i=0; i<XLdistrib.getWinSize(); i++){
//        							//Note that the mean (i.e. position) is at position zero w.r.t XL distribution 
//        							sumDiffW += (xlW[i]*(l+i-offsetW)*(l+i-offsetW)) + 
//    		    						    (xlC[XLdistrib.getWinSize()-i-1]*(l+i-offsetW)*(l+i-offsetW)); //note... -oldOffsetW is correct here
//        							totW+=xlW[i]+xlC[XLdistrib.getWinSize()-i-1];
//        						}
//        						sumDiffC = sumDiffW; totC = totW;
//        					}else{
        						for(int i=0; i<XLdistrib.getWinSize(); i++){
        							//Note that the mean (i.e. position) is at position zero w.r.t XL distribution 
        							sumDiffW += (xlW[i]*(l+i-offsetW)*(l+i-offsetW));
        							totW+=xlW[i];
        							sumDiffC += (xlC[i]*(l+i-offsetC)*(l+i-offsetC));
        							totC+=xlC[i];
        						}
//        					}
        					if(sumDiffW>0)
        						newSigmaW = Math.sqrt(sumDiffW/totW);
        					if(sumDiffC>0)
        						newSigmaC = Math.sqrt(sumDiffC/totC);
        					XLdistrib.loadGaussianDistrib(-config.getXLDistribOffset(), newSigmaW, config.getXLDistribOffset(), newSigmaC);
        				}else{
        					/**
        					//Calculate new mean and sigma (symmetric)
        					fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
        					for(int i=0; i<XLdistrib.getWinSize(); i++){
        						fitter.addObservedPoint((double)(i+XLdistrib.getLeft()), xlW[i]+xlC[XLdistrib.getWinSize()-i-1]);
        					}
        					parameters = fitter.fit();;
        					newOffset = parameters[1];
        					newSigma = parameters[2];
        					XLdistrib.loadGaussianDistrib(newOffset, newSigma, -newOffset, newSigma); //symmetric
							**/
        					
        					//Calculate new mean and sigma (non-symmetric)
        					double newSigmaW=config.getXLDistribSigma(), newSigmaC=config.getXLDistribSigma(),
        							newOffsetW = XLdistrib.getGaussOffsetW(), newOffsetC = XLdistrib.getGaussOffsetC();
        					// Watson
        					fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
        					for (int i=0; i < XLdistrib.getWinSize(); i++){
        						fitter.addObservedPoint((double)(i+XLdistrib.getLeft()), xlW[i]);
        					}
        					parameters = fitter.fit();;
        					newOffsetW = parameters[1];
        					newSigmaW = parameters[2];
        					//Crick
        					fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
        					for (int i=0; i < XLdistrib.getWinSize(); i++){
        						fitter.addObservedPoint((double)(i+XLdistrib.getLeft()), xlC[i]);
        					}
        					parameters = fitter.fit();;
        					newOffsetC = parameters[1];
        					newSigmaC = parameters[2];
        					
        					XLdistrib.loadGaussianDistrib(newOffsetW, newSigmaW, newOffsetC, newSigmaC); //non-symmetric
        				}
        			}
        		}
        	}
        	
        	////////////
            //Plot the current pi & priors if plotting
        	////////////
        	if(plotEM && plotter!=null){
        		for(ExperimentCondition cond : manager.getConditions()){
        			int c = cond.getIndex();
        			String outName = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EM_r"+trainingRound+"_t"+(t+1)+"_"+cond.getName() : null;
        			String outNameZoom = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EMzoom_r"+trainingRound+"_t"+(t+1)+"_"+cond.getName() : null;
        			int trimLeft = (model.getWidth()/2)-50;
        			int trimRight = (model.getWidth()/2)-50;
//        			fullImages.add(plotter.plotCompositeEM(outName, composite, model, mu[c], pi[c], trainingRound, t, 0,0));
//        			zoomImages.add(plotter.plotCompositeEM(outNameZoom, composite, model, mu[c], pi[c], trainingRound, t, trimLeft, trimRight));
        		}
            }

            //Is current state equivalent to the last?
            if( t>config.ALPHA_ANNEALING_ITER && lastEquivToCurr())
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
            if (nonZeroComps>0 && ((numConditions==1 && t<=config.ALPHA_ANNEALING_ITER) || (config.CALC_LL && Math.abs(LAP-lastLAP)>config.EM_CONVERGENCE) || stateEquivCount<config.EM_STATE_EQUIV_ROUNDS)){
                copyStateToLast();
                lastLAP = LAP;
                continue;
            }else{
            	copyStateToLast();
            	lastLAP = LAP;
            	break;
            }
        } //LOOP: Run EM while not converged
    }//end of EM_MAP method
    
    
    /**
     * Responsibility assignment in the composite using trained model. No maximization. 
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void responsibilityAssignment () {
        int numComp = numComponents;
        double [][] totalResp = new double[numConditions][];
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		int numBases = hitNum[c];
    		totalResp[c] = new double[numBases];
            for(int i=0;i<numBases;i++)
                totalResp[c][i] = 0;
    	}        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run E step once
    	//////////////////////////////////////////////////////////////////////////////////////////
        for(int c=0; c<numConditions; c++){ 
    		//Recompute h function, given binding component positions
    		for(int i=0;i<hitNum[c];i++)
            	for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
            		int dist = hitPos[c][i]-mu[c][j];
            		h[c][j][i] = model.getAllComponents(c).get(j).getTagDistribution().probability(dist, hitPlusStr[c][i]);
            	}}
    		//Compute responsibilities
			for(int i=0;i<hitNum[c];i++)
                totalResp[c][i] = 0;
    		for(int i=0;i<hitNum[c];i++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				r[c][j][i] = h[c][j][i]*pi[c][j];
    				totalResp[c][i] +=r[c][j][i]; 
    			}}
    		}
    		//Normalize responsibilities
    		for(int i=0;i<hitNum[c];i++){
    			for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
    				r[c][j][i]/=totalResp[c][i];
    			}}
    		}
		}
    }//end of responsibilityAssignment method
    
    /**
     * Set responsibility profile for each component (for kernel update)
     * @param bindComponents
     * @param signals
     * @param responsibilities
     */
    private void setComponentResponsibilityProfiles(double[][][] responsibilities) {
    	for(ExperimentCondition cond : manager.getConditions()){
			int c = cond.getIndex();
			for(int j=0;j<numComponents;j++){
				CompositeModelComponent comp = model.getAllComponents(c).get(j);
				int width = comp.getTagDistribution().getWinSize();
				int center = -1*comp.getTagDistribution().getLeft();
        	
		    	double[][] rc = responsibilities[c];
		   
		   		// store binding profile (read responsibilities in c condition) of this component
				double[] profile_plus = new double[width];
				double[] profile_minus = new double[width];
				
				for(int i=0;i<hitNum[c];i++){
					int offset = hitPos[c][i]-comp.getPosition()+center;
					if(offset>=0 && offset<width)
						if (hitPlusStr[c][i])
							profile_plus[offset]=rc[j][i]*hitCounts[c][i];
						else
							profile_minus[offset]=rc[j][i]*hitCounts[c][i];
				}
				comp.setTagProfile(profile_plus,  true);
				comp.setTagProfile(profile_minus, false);
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
	    		for(int x=0; x<r[c][j].length; x++){
	    			lastRBind[c][j][x] = r[c][j][x];
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
    	for(int c=0; c<numC; c++){
    		for(int j=0;j<numComponents;j++){
    			if(pi[c][j]>0)
    				currNZ++;
    			if(lastPi[c][j]>0)
    				lastNZ++;
    		}
    	}
    	boolean numCompEqual = currNZ==lastNZ;
    	boolean compPosEqual=true;
    	if(numCompEqual){
    		for(int c=0; c<numC; c++){
    			for(int j=0; j<mu[c].length; j++){if(pi[c][j]>0){
    				compPosEqual = compPosEqual && (mu[c][j] == lastMu[c][j]);
    			}}}
    		}else{
    			compPosEqual=false;
    		}
    		boolean piBindEquivalent=true;
    		for(int c=0; c<numC; c++){
    			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    				piBindEquivalent = piBindEquivalent && (Math.abs(pi[c][j]-lastPi[c][j])<config.EM_STATE_EQUIV_THRES);
    		}}}
    		boolean rBindEquivalent=true;
    		for(int c=0; c<numC; c++){
    			for(int j=0; j<pi[c].length; j++){if(pi[c][j]>0){
    				for(int x=0; x<r[c][j].length; x++){
    					rBindEquivalent = rBindEquivalent && (Math.abs(r[c][j][x]-lastRBind[c][j][x])<config.EM_STATE_EQUIV_THRES);
    				}
    			}
    		}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
    
	/**
	 * Print responsibilities for each base in the composite to the file
	 * @param cond
	 * @param filename
	 */
	private void printResponsibilitiesToFile(ExperimentCondition cond, String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			int c = cond.getIndex();
			fout.write("#Pos\tTotal");
			for(int j=0;j<numComponents;j++){ if(pi[c][j]>0){
				fout.write("\t"+j);
			}}fout.write("\n");
			
			for(int i=0;i<hitNum[c];i++){
				fout.write(hitPos[c][i]+"\t"+hitCounts[c][i]);
    			for(int j=0;j<numComponents;j++){ if(pi[c][j]>0){
    				double resp = hitPlusStr[c][i] ? r[c][j][i]*hitCounts[c][i] : -r[c][j][i]*hitCounts[c][i]; 
    				fout.write("\t"+resp);
    			}}
    			fout.write("\n");
    		}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	/**
	 * Make Gifs from stored images
	 */
	public void makeGifs(){
		GifSequenceWriter gifWriter, gifZoomWriter;
		try {
			if(fullImages.size()>0 && zoomImages.size()>0){
				//GIFs
		    	String gifName = config.getOutputImagesDir()+File.separator+"EM.gif";
		    	String gifNameZoom = config.getOutputImagesDir()+File.separator+"EMzoom.gif";
		    	ImageOutputStream output1 = new FileImageOutputStream(new File(gifName));
				ImageOutputStream output2 = new FileImageOutputStream(new File(gifNameZoom));
		    	gifWriter = new GifSequenceWriter(output1, fullImages.get(0).getType(), 5, false);
		    	gifZoomWriter = new GifSequenceWriter(output2, zoomImages.get(0).getType(), 5, false);
		    	
		    	//Write the images
		    	for(int i=0; i<20; i++){ //Pause on the first image for 2 seconds
		    		gifWriter.writeToSequence(fullImages.get(0));
		    		gifZoomWriter.writeToSequence(zoomImages.get(0));
		    	}
		    	for(BufferedImage im : fullImages)
		    		gifWriter.writeToSequence(im);
		    	for(BufferedImage im : zoomImages)
		    		gifZoomWriter.writeToSequence(im);
		    	for(int i=0; i<300; i++){ //Pause on the last image for 15 seconds
		    		gifWriter.writeToSequence(fullImages.get(fullImages.size()-1));
		    		gifZoomWriter.writeToSequence(zoomImages.get(zoomImages.size()-1));
		    	}
		    	
		    	gifWriter.close();
		    	gifZoomWriter.close();
		    	output1.close();
		    	output2.close();
		    	
		    	
				//Final frames only
		    	String finalFrameName = config.getOutputImagesDir()+File.separator+"EM_final.png";
		    	String finalFrameNameZoom = config.getOutputImagesDir()+File.separator+"EMzoom_final.png";
		    	ImageIO.write(fullImages.get(fullImages.size()-1),"png",new File(finalFrameName));
		    	ImageIO.write(zoomImages.get(zoomImages.size()-1),"png",new File(finalFrameNameZoom));
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
