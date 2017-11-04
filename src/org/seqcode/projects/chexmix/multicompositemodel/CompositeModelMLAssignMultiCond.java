package org.seqcode.projects.chexmix.multicompositemodel;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.chexmix.composite.CompositeModelComponent;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;


/**
 * CompositeModelMLAssign: run ML assignment on a composite tag distribution's collection of individual points.
 * 
 * This method assumes that the same protein-DNA interaction model is valid (i.e. constant) across all examined conditions. 
 * 
 * TODO: this class should be made threadable. Be careful to instantiate one of these objects in/as each thread (i.e. don't do multithreaded calls to assign()
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CompositeModelMLAssignMultiCond {

	protected ExperimentManager manager;
	protected ChExMixConfig config;
	protected CompositeTagDistribution composite;
	protected ProteinDNAInteractionModelMultiCond model;
	protected int numComponents;  //The count of all components (active +inactive) in the model
	protected int numConditions;
	//	ML VARIABLES
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
	protected double[][][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[][]   lastPi;		//Last Pi (monitor convergence)
	protected int[][]      lastMu;		//Last positions (monitor convergence)
	protected double lastLAP, LAP; 		//log-likelihood monitoring
	protected boolean plotRegions=false;		//Plot the region components
	protected List<Region> regionsToPlot=null; //Plot these points
	protected int stateEquivCount=0;
	
	/**
	 * Constructor
	 * 
	 * @param c: XLAnalysisConfig
	 * @param eMan: ExperimentManager
	 * @throws Exception 
	 */
	public CompositeModelMLAssignMultiCond(CompositeTagDistribution composite, ProteinDNAInteractionModelMultiCond model, ChExMixConfig config, ExperimentManager eMan){
		this.composite=composite;
		this.model=model;
		this.config=config;
		manager = eMan;
		numConditions = manager.getNumConditions();
		
		//Initialize data structures
		numComponents = model.getNumComponents();
        hitCounts= new double[numConditions][];	// Hit weights
        hitPos= new int[numConditions][];			// Hit positions
    	hitPlusStr= new boolean[numConditions][];	// Hit positive strand boolean
    	hitNum = new int[numConditions];			// Number of hits in each condition
    	h= new double[numConditions][][]; 			// H function (binding component probability per read)
    	r= new double[numConditions][][];		// Binding component responsibilities
    	pi = new double[numConditions][numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numConditions][numComponents]; // mu : positions of the binding components (fixed across conditions)
    	plotRegions = config.getRegionsToPlot().size()>0;
    	regionsToPlot = config.getRegionsToPlot();
        lastRBind = new double[numConditions][][];
        lastPi = new double[numConditions][numComponents];
        lastMu = new int[numConditions][numComponents];
        
        //All calls to this ML assigner will necessarily have:
        // - same number of positions (hitNum)
        // - same relative positions (hitPos), due to model size staying constant
        // - same split of plus/minus coordinates (hitPlusStr)
        // - same responsibility function relationships to positions (h)
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
        	hitNum[c]=composite.getWinSize()*2;
        	int[] posc= new int[hitNum[c]];
            boolean[] plusc= new boolean[hitNum[c]];
            for(int i=0;i<composite.getWinSize();i++){ //Watson
            	posc[i] = i;
            	plusc[i] = true;
            }
            for(int i=0;i<composite.getWinSize();i++){ //Crick
            	posc[i+composite.getWinSize()] = i;
            	plusc[i+composite.getWinSize()] = false;
            }
            hitPos[c] = posc;
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
        
	}
	
	
	/**
     * ML assignment
     *
     * Almost purely matrix/array operations.
     * 
     * Input are a pair of watson & crick tag profiles corresponding to a genomic locus. 
     * Dimensions of tag profiles are: number of conditions x width of profile.
     * Width of profiles needs to be equal to the model width.
     * 
     * UpdatePi: flag to update the component probabilities (true ML) versus straight assignment using trained weights. 
     *  
     * Returns an updated protein-DNA interaction model
     */
    public CompositeModelSiteAssignmentMultiCond  assignToPointFromComposite(int pointIndex, boolean updatePi
    											  ) throws Exception{
    	
    	//Need to ensure that the composite and model have the same center offsets, or the coordinate system will not be consistent
    	if(composite.getWinSize()!=model.getWidth() || composite.getCenterOffset()!=model.getCenterOffset())
    		throw new Exception("CompositeModelMLAssign: Composite distribution coordinate system not consistent with protein-DNA interaction model");
    	
    	List<List<CompositeModelComponent>> nonZeroComponents = new ArrayList<List<CompositeModelComponent>>();
    	for (int c=0; c < numConditions ; c++)
    		nonZeroComponents.add(model.getNonZeroComponents(c));
    	
    	//Reset data structures
    	for (int c=0; c < numConditions ; c++){
    		for(int j=0;j<numComponents;j++){
    			//Load pi for binding components
    			CompositeModelComponent comp = model.getAllComponents(c).get(j);
    			pi[c][j]= comp.getPi();
    			//Load binding component positions
    			mu[c][j] = model.getAllComponents(c).get(j).getPosition();
    		}
        }
    	//Condition-specific tag info
        double[] hitCountTotals= new double[numConditions];	// Hit weight totals
        for(ExperimentCondition cond : manager.getConditions()){
        	int c = cond.getIndex();
            double[] countc= new double[hitNum[c]];
            double[][] pointWatson = composite.getPointWatsons(pointIndex);
            double[][] pointCrick = composite.getPointCricks(pointIndex);
            double countctot = 0.0;
            for(int i=0;i<composite.getWinSize();i++){ //Watson
                countc[i]=pointWatson[c][i];
                countctot+=countc[i];
            }
            for(int i=0;i<composite.getWinSize();i++){ //Crick
                countc[i+composite.getWinSize()]=pointCrick[c][i];
                countctot+=countc[i+composite.getWinSize()];
            }
            hitCounts[c]=countc;
            hitCountTotals[c] = countctot;
            
        }
        //End of data structure initialization
        
        
        //////////
        // Run EM steps
        //////////
        ML(composite.getPoint(pointIndex), updatePi);
	
        //////////
        // Define the results of ML for this point
        //////////
        int [][] compIndices = new int[numConditions][];
        double[][] compResp = new double[numConditions][];
        for (int c=0; c < numConditions ; c++){
    		compIndices[c] = new int[nonZeroComponents.get(c).size()];
    		compResp = new double[numConditions][nonZeroComponents.get(c).size()];
        }
    	for(int c=0; c<numConditions;c++){
    		int x=0;
    		for(CompositeModelComponent comp : nonZeroComponents.get(c)){ 
    			int j = comp.getIndex();
    			compIndices[c][x]=j;
    		
            	double sumRespW=0.0, sumRespC=0.0;
            	for(int i=0;i<hitNum[c];i++){
            		if(hitPlusStr[c][i])
            			sumRespW += hitCounts[c][i]*r[c][j][i];
            		else
            			sumRespC += hitCounts[c][i]*r[c][j][i];
            	}
            	compResp[c][x] = sumRespW+sumRespC;
            	x++;
            }
    	}
    	CompositeModelSiteAssignmentMultiCond currSiteAssign = new CompositeModelSiteAssignmentMultiCond(composite.getPoint(pointIndex), hitCountTotals, compIndices, compResp);
    	
        return currSiteAssign;
    }//end of EMTrain method


    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void ML (StrandedPoint currPoint, boolean updatePi) {
        int numComp = numComponents;
        double [][] totalResp = new double[numConditions][];
        boolean plotThisPoint=false;
        
        //Initialize responsibilities
        for(int c=0; c<numConditions; c++){
    		int numBases = hitNum[c];
    		totalResp[c] = new double[numBases];
            for(int i=0;i<numBases;i++)
                totalResp[c][i] = 0;
    	}
        

    	////////////
        //Plot the initial pi & priors if plotting
    	////////////
        if(plotRegions && currPoint!=null){
        	//Check if current Point is in the list of regions to print
        	for(Region preg : regionsToPlot)
        		if(preg.contains(currPoint))
        			plotThisPoint =true;
        }
        if(plotThisPoint){
        	/*String regStr = currRegion.getLocationString().replaceAll(":", "-");
        	String outName = "EM_"+regStr+"_r"+trainingRound+"_t0";
        	int trimLeft = Math.max(0, plotSubRegion.getStart()-regStart);
        	int trimRight = Math.max(0, currRegion.getEnd()-plotSubRegion.getEnd());
        	
        	if(config.useMotifPrior())
       			EMStepPlotter.execute(outName, currRegion, mu, pi, motifPrior, numConditions, numComponents, 0, trimLeft, trimRight);
       		else
       			EMStepPlotter.execute(outName, currRegion, mu, pi, null, numConditions, numComponents, 0, trimLeft, trimRight);
       		*/
        }
        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run ML while not converged
        // Note: iterations during which we eliminate a binding component don't count towards "t"
    	//////////////////////////////////////////////////////////////////////////////////////////
        int t=0;
        while(t<config.ML_ML_ITER){  
        	
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
    		//M-step: maximize pi
    		/////////////////////
    		if(updatePi){
    			for(int c=0; c<numConditions; c++){
    				double[] sumR=new double[numComponents];
    				for(int j=0;j<numComp;j++){ if(pi[c][j]>0){
	    				for(int i=0;i<hitNum[c];i++)
	    					sumR[j] += r[c][j][i]*hitCounts[c][i];
	            }}
	
	    		//Update pi(j)
	            for(int j=0;j<numComp;j++){ if(pi[c][j]>0 && model.getAllComponents(c).get(j).hasUpdatablePi()){
	            	pi[c][j]=Math.max(0, sumR[j]); 
	            }}
	            
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
    		}}
            
            
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
	            	
	            	LP+=-(sum_log_pi);
	            }
	            LAP = LL+LP;
	
	            System.out.println("ML: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
            }
        	
        	////////////
            //Plot the current pi & priors if plotting
        	////////////
            if(plotThisPoint){
            	/*String regStr = currRegion.getLocationString().replaceAll(":", "-");
            	String outName = "EM_"+regStr+"_r"+trainingRound+"_t"+(t+1)+"_i"+(iter+1);
            	int trimLeft = Math.max(0, plotSubRegion.getStart()-regStart);
            	int trimRight = Math.max(0, currRegion.getEnd()-plotSubRegion.getEnd());

           		if(config.useMotifPrior())
           			EMStepPlotter.execute(outName, currRegion, mu, pi, motifPrior, numConditions, numComponents, t+1, trimLeft, trimRight);
           		else
           			EMStepPlotter.execute(outName, currRegion, mu, pi, null, numConditions, numComponents, t+1, trimLeft, trimRight);
           		*/
            }

            //Is current state equivalent to the last?
            if(lastEquivToCurr())
            	stateEquivCount++;
            else
            	stateEquivCount=0;
                        

    		//Tick the clock forward
    		t++;
    		
            ////////////
          	//Check Stopping condition
          	////////////
    		if (updatePi || (nonZeroComps>0 && (t==0 || !lastEquivToCurr()))){
                copyStateToLast();
                lastLAP = LAP;
                continue;
            }else{
            	copyStateToLast();
            	lastLAP = LAP;
            	break;
            }
        } //LOOP: Run ML while not converged
    }//end of ML method
    
    
    
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
    
	
}
