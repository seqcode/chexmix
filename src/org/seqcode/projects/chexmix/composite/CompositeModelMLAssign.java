package org.seqcode.projects.chexmix.composite;

import java.util.List;

import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
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
public class CompositeModelMLAssign {

	protected ExperimentCondition condition;
	protected XLAnalysisConfig config;
	protected CompositeTagDistribution composite;
	protected ProteinDNAInteractionModel model;
	protected int numComponents;  //The count of all components (active +inactive) in the model
	//	ML VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[]     hitCounts;	// Hit weights
	protected int[]        hitPos;		// Hit positions
	protected boolean[]    hitPlusStr;	// Hit positive strand boolean
	protected int		   hitNum;		// Number of hits in each condition 
	protected double[][]   hAll;		// H function values for all positions in the current window (precomputed)
	protected double[][]   h; 			// H function (binding component probability per read)
	protected double[][]   r;		// Binding component responsibilities
	protected double[]     pi;			// pi : emission probabilities for binding components
	protected int[]        mu;			// mu : positions of the binding components
	protected double[][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[]   lastPi;		//Last Pi (monitor convergence)
	protected int[]      lastMu;		//Last positions (monitor convergence)
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
	public CompositeModelMLAssign(CompositeTagDistribution composite, ProteinDNAInteractionModel model, XLAnalysisConfig config, ExperimentCondition cond){
		this.composite=composite;
		this.model=model;
		this.config=config;
		condition = cond;
		
		//Initialize data structures
		numComponents = model.getNumComponents();
    	pi = new double[numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numComponents]; // mu : positions of the binding components (fixed across conditions)
    	plotRegions = config.getRegionsToPlot().size()>0;
    	regionsToPlot = config.getRegionsToPlot();
        lastPi = new double[numComponents];
        lastMu = new int[numComponents];
        
        //All calls to this ML assigner will necessarily have:
        // - same number of positions (hitNum)
        // - same relative positions (hitPos), due to model size staying constant
        // - same split of plus/minus coordinates (hitPlusStr)
        // - same responsibility function relationships to positions (h)
        hitNum=composite.getWinSize()*2;
        hitPos= new int[hitNum];
        hitPlusStr= new boolean[hitNum];
        for(int i=0;i<composite.getWinSize();i++){ //Watson
        	hitPos[i] = i;
        	hitPlusStr[i] = true;
        }
        for(int i=0;i<composite.getWinSize();i++){ //Crick
        	hitPos[i+composite.getWinSize()] = i;
        	hitPlusStr[i+composite.getWinSize()] = false;
        }
        	
        //Initialize responsibility functions
        h= new double[numComponents][hitNum];
        for(int i=0;i<hitNum;i++){
        	for(int j=0;j<numComponents;j++){
            	int dist = hitPos[i]-mu[j];
                h[j][i] = model.getAllComponents().get(j).getTagDistribution().probability(dist, hitPlusStr[i]);
            }
        }
        r = new double[numComponents][hitNum];
        lastRBind = new double[numComponents][hitNum];
        
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
    public CompositeModelSiteAssignment  assignToPointFromComposite(int pointIndex, boolean updatePi
    											  ) throws Exception{
    	
    	//Need to ensure that the composite and model have the same center offsets, or the coordinate system will not be consistent
    	if(composite.getWinSize()!=model.getWidth() || composite.getCenterOffset()!=model.getCenterOffset())
   			throw new Exception("CompositeModelMLAssign: Composite distribution coordinate system not consistent with protein-DNA interaction model");
    			
    	List<CompositeModelComponent> nonZeroComponents = model.getNonZeroComponents();
    	
    	//Reset data structures
        for(int j=0;j<numComponents;j++){
        	//Load pi for binding components
        	CompositeModelComponent comp = model.getAllComponents().get(j);
            pi[j]= comp.getPi();
            //Load binding component positions
        	mu[j] = model.getAllComponents().get(j).getPosition();
        }
    	//Condition-specific tag info
        double hitCountTotals = 0.0;	// Hit weight totals
        hitCounts = new double[hitNum];
        double[] pointWatson = composite.getPointWatsons(pointIndex);
        double[] pointCrick = composite.getPointCricks(pointIndex);
        for(int i=0;i<composite.getWinSize();i++){ //Watson
        	hitCounts[i]=pointWatson[i];
        	hitCountTotals+=hitCounts[i];
        }
        for(int i=0;i<composite.getWinSize();i++){ //Crick
        	hitCounts[i+composite.getWinSize()]=pointCrick[i];
        	hitCountTotals+=hitCounts[i+composite.getWinSize()];
        }

        //End of data structure initialization
        
        
        //////////
        // Run EM steps
        //////////
        ML(composite.getPoint(pointIndex), updatePi);
	
        //////////
        // Define the results of ML for this point
        //////////
    	int[] compIndices = new int[nonZeroComponents.size()];
    	double[] compResp = new double[nonZeroComponents.size()];
    	int x=0;
    	for(CompositeModelComponent comp : nonZeroComponents){ 
    		int j = comp.getIndex();
    		compIndices[x]=j;
    		
            double sumRespW=0.0, sumRespC=0.0;
            for(int i=0;i<hitNum;i++){
            	if(hitPlusStr[i])
            		sumRespW += hitCounts[i]*r[j][i];
            	else
            		sumRespC += hitCounts[i]*r[j][i];
            }
            compResp[x] = sumRespW+sumRespC;

            x++;
    	}
    	CompositeModelSiteAssignment currSiteAssign = new CompositeModelSiteAssignment(composite.getPoint(pointIndex), hitCountTotals, compIndices, compResp);
    	
        return currSiteAssign;
    }//end of EMTrain method


    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void ML (StrandedPoint currPoint, boolean updatePi) {
        int numComp = numComponents;
        boolean plotThisPoint=false;
        
        //Initialize responsibilities
    	double[] totalResp = new double[hitNum];
        for(int i=0;i<hitNum;i++)
            totalResp[i] = 0;
        

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
        	//Recompute h function, given binding component positions
        	for(int i=0;i<hitNum;i++)
                for(int j=0;j<numComp;j++){ if(pi[j]>0){
                	int dist = hitPos[i]-mu[j];
                	h[j][i] = model.getAllComponents().get(j).getTagDistribution().probability(dist, hitPlusStr[i]);
                }}
        	//Compute responsibilities
    		for(int i=0;i<hitNum;i++)
                   totalResp[i] = 0;
        	for(int i=0;i<hitNum;i++){
        		for(int j=0;j<numComp;j++){ if(pi[j]>0){
        			r[j][i] = h[j][i]*pi[j];
        			totalResp[i] +=r[j][i]; 
        		}}
        	}
        	//Normalize responsibilities
        	for(int i=0;i<hitNum;i++){
        		for(int j=0;j<numComp;j++){ if(pi[j]>0){
        			r[j][i]/=totalResp[i];
        		}}
        	}
    		
        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
    		if(updatePi){
	    		double[] sumR=new double[numComponents];
	    		for(int j=0;j<numComp;j++){ if(pi[j]>0){
	    			for(int i=0;i<hitNum;i++)
	    				sumR[j] += r[j][i]*hitCounts[i];
	            }}
	
	    		//Update pi(j)
	            for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePi()){
	            	pi[j]=Math.max(0, sumR[j]); 
	            }}
	            
	            //Normalize pi
	            double totalRPi=0, totalNonUpdateablePi=0;
	            for(int j=0;j<numComp;j++){ if(pi[j]>0 && !model.getAllComponents().get(j).hasUpdatablePi()){
	        		totalNonUpdateablePi+=pi[j];
	        	}}
	            for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePi()){
	        		totalRPi+=pi[j];
	        	}}
	            for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePi()){
	        		if(totalRPi>0){
	        			pi[j]/=totalRPi;
	        			pi[j]*=1-totalNonUpdateablePi;
	        		}
	        	}}
    		}
            
            
        	//Non-zero components count
        	int nonZeroComps=0;
        	for(int j=0;j<numComp;j++)
        		if(pi[j]>0)
        			nonZeroComps++;
        	
        	////////////
        	//Compute LL
        	////////////
        	LAP=0;
        	if(config.CALC_LL){
	        	//Log-likelihood calculation
	            double LL =0;
	        	for(int i=0;i<hitNum;i++){
	        		// for each read, each event will give a conditional prob or bg prob
	                double j_sum=0;
	                for(int j=0;j<numComp;j++){ if(pi[j]>0.0){
	        			j_sum += Math.log(r[j][i])/config.LOG2;
	                }}
	                LL += j_sum*hitCounts[i];                        
	            }

	            //Log priors
	            double LP=0;
	            //sum of pi
	            double sum_log_pi=0;
	            for(int j=0;j<numComp;j++){ if(pi[j]>0.0){
	            	sum_log_pi+=Math.log(pi[j])/config.LOG2;
	            }}
	            	
	            LP+=-(sum_log_pi);

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
		for(int j=0; j<numComponents; j++){
			lastPi[j] = pi[j];
			lastMu[j] = mu[j];
	    	for(int x=0; x<r[j].length; x++){
	    		lastRBind[j][x] = r[j][x];
	    	}
		}
    }
    
    /**
     * Compare last variables to current (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     * @return
     */
    private boolean lastEquivToCurr(){
    	int currNZ=0, lastNZ=0;
		for(int j=0;j<numComponents;j++){
			if(pi[j]>0)
				currNZ++;
			if(lastPi[j]>0)
				lastNZ++;
		}
    	boolean numCompEqual = currNZ==lastNZ;
    	boolean compPosEqual=true;
    	if(numCompEqual){
			for(int j=0; j<mu.length; j++){if(pi[j]>0){
				compPosEqual = compPosEqual && (mu[j] == lastMu[j]);
			}}
    	}else{
    		compPosEqual=false;
    	}
    	boolean piBindEquivalent=true;
		for(int j=0; j<pi.length; j++){if(pi[j]>0){
			piBindEquivalent = piBindEquivalent && (Math.abs(pi[j]-lastPi[j])<config.EM_STATE_EQUIV_THRES);
		}}
    	boolean rBindEquivalent=true;
		for(int j=0; j<pi.length; j++){if(pi[j]>0){
			for(int x=0; x<r[j].length; x++){
				rBindEquivalent = rBindEquivalent && (Math.abs(r[j][x]-lastRBind[j][x])<config.EM_STATE_EQUIV_THRES);
			}
		}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
    
	
}
