package org.seqcode.projects.chexmix.mixturemodel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.chexmix.events.BindingEvent;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingModel;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;



/**
 * BindingMLAssignment: Maximum likelihood assignment of reads to a configuration of binding components.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class MultiGPSMLAssignment {

	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected ChExMixConfig config;
	protected List<BindingSubComponents> components;
	protected List<NoiseComponent> noise;
	protected int numComponents;  //Assumes the same number of active+inactive components in each condition
	protected int numConditions;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds; //Background models per condition
	//	EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[]   sigHitCounts;	// Hit weights
	protected int[]      sigHitPos;		// Hit positions
	protected boolean[]  sigHitPlusStr;	// Hit positive strand boolean
	protected int		 sigHitNum;		// Number of hits in each condition 
	protected int[]      sigRepIndices;  // Index of replicate for the hit
	protected double[]   ctrlHitCounts;	// Hit weights
	protected int[]      ctrlHitPos;		// Hit positions
	protected boolean[]  ctrlHitPlusStr;	// Hit positive strand boolean
	protected int		 ctrlHitNum;		// Number of hits in each condition 
	protected int[]      ctrlRepIndices;  // Index of replicate for the hit
	protected double[][] h; 			// H function (binding component probability per read)
	protected double[]   n; 			// N function (noise component probability per read)
	protected double[][] rBindSig;		// Binding component responsibilities (signal reads)
	protected double[]   rNoiseSig;		// Noise component responsibilities (signal reads)
	protected double[][] rBindCtrl;		// Binding component responsibilities (control reads)
	protected double[]   rNoiseCtrl;		// Noise component responsibilities (control reads)
	protected double[]   pi;			// pi : emission probabilities for binding components
	protected double     piNoise;		// pi : emission probabilities for noise components (fixed)
	protected int[]      mu;			// mu : positions of the binding components
	protected double[][] compLL;		//Log-likelihood for each component in each condition
	protected double[][] lastRBind;	//Last responsibilities (monitor convergence)
	protected double[]   lastPi;		//Last Pi (monitor convergence)
	protected int[]      lastMu;		//Last positions (monitor convergence)
	protected double[]   tmp_pi;			// pi used in ML calc
	protected double[][] tmp_h;			// h used in ML calc
	protected double[][] tmp_rBindSig;	// rBindSig used in ML calc
	protected double[]   tmp_rNoiseSig;	// rNoiseSig used in ML calc
	protected double	 tmp_piNoise; 	// piNoise used in ML calc
	protected BindingModel[] bindingModels; //Array of binding models for convenience
	protected double[]	   sigRepHitCountTotals; //Hit count totals counted by replicate (for convenience)
	protected double[]	uniformRepHitCountTotals; //Hit count totals by replicate if signal read counts were distributed uniformly (used only if there is no control) 
	protected double numPotentialRegions;
	
	/**
	 * Constructor
	 * @param c
	 * @param eMan
	 */
	public MultiGPSMLAssignment(ExptConfig econ, EventsConfig evcon, ChExMixConfig c, ExperimentManager eMan, BindingManager bindMan, HashMap<ExperimentCondition, BackgroundCollection> condBacks, int numPotReg){
		config=c;
		evconfig = evcon;
		econfig = econ;
		manager = eMan;
		bindingManager = bindMan;
		conditionBackgrounds = condBacks;
		numConditions = manager.getNumConditions();
		numPotentialRegions = (double)numPotReg;
	}
	
	/**
     * ML assignment
     *
     * Takes as input a SINGLE list of binding components (i.e. a single configuraton)
     * Returns lists of binding events indexed by condition. 
     * Pi values are calculated from the signal hits and then those same components are directly applied to the control hits. 
     *
     * Almost purely matrix/array operations.
     */
    public List<BindingEvent>  assign(List<List<StrandedBaseCount>> signals,
    								  List<List<StrandedBaseCount>> controls,
    								  Region w, 
    								  List<NoiseComponent> noise,
    								  List<BindingSubComponents> comps, 
    								  int numComp,
    								  ExperimentCondition cond){ 
    	components = comps;
        this.noise = noise;
        numComponents = numComp;
        //Matrix initializations
    	pi = new double[numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numComponents];// mu : positions of the binding components
    	compLL = new double [numConditions][numComponents];		//Log-likelihood for each component in each condition
        bindingModels = new BindingModel[manager.getReplicates().size()]; //Array of bindingModels for convenience
        sigRepHitCountTotals = new double[manager.getReplicates().size()]; //Hit count totals counted by replicate (for convenience)
        uniformRepHitCountTotals = new double[manager.getReplicates().size()]; //Hit count totals by replicate if reads were distributed uniformly 
        //Monitor state convergence using the following last variables
        lastPi = new double[numComponents];
        lastMu = new int[numComponents];
        //Temporary variables
        tmp_pi = new double[numComponents];			// pi used in ML calc
        
        //Initializing data structures
        int c = cond.getIndex();
        	
        //Add bindingModels to array
        for(ControlledExperiment rep : cond.getReplicates())
        	bindingModels[rep.getIndex()] = bindingManager.getUnstrandedBindingModel(rep);
        	
        //Load Reads (merge from all replicates)
        List<StrandedBaseCount> sigBases = new ArrayList<StrandedBaseCount>();
        List<StrandedBaseCount> ctrlBases = new ArrayList<StrandedBaseCount>();
        for(ControlledExperiment rep : cond.getReplicates()){
        	sigBases.addAll(signals.get(rep.getIndex()));
        	if(controls.get(rep.getIndex())!=null)
        		ctrlBases.addAll(controls.get(rep.getIndex()));
        }
        sigHitNum = sigBases.size();
        ctrlHitNum = ctrlBases.size();
        	
        //Count total weights for convenience
        for(ControlledExperiment rep : cond.getReplicates()){
        	sigRepHitCountTotals[rep.getIndex()]=0;
        	for(StrandedBaseCount s : signals.get(rep.getIndex()))
        		sigRepHitCountTotals[rep.getIndex()]+=s.getCount();
        	uniformRepHitCountTotals[rep.getIndex()] = (((rep.getSignal().getHitCount()*(1-rep.getSignalVsNoiseFraction()))/econfig.getMappableGenomeLength())*(double)w.getWidth())/rep.getControlScaling(); //Normalizing by control scaling is a hack - usually control scaling will be 1 when the replicate has no control... however, it is not 1 for SES. 
        }
        	
        //Load replicate index for each read
        sigRepIndices = new int[sigHitNum];
        ctrlRepIndices = ctrlHitNum==0 ? null : new int[ctrlHitNum];
        int ys=0, yc=0, z=0;
        for(ControlledExperiment rep : cond.getReplicates()){
        	z=0;
        	while(z<signals.get(rep.getIndex()).size()){
        		sigRepIndices[ys] = rep.getIndex();
        		z++; ys++;
        	}
        	z=0;
        	while(z<controls.get(rep.getIndex()).size()){
        		ctrlRepIndices[yc] = rep.getIndex();
        		z++; yc++;
        	}
        }
            
        //Load signal read info
        sigHitCounts= new double[sigHitNum];
        sigHitPos= new int[sigHitNum];
        sigHitPlusStr= new boolean[sigHitNum];
        for(int i=0;i<sigHitNum;i++){
        	sigHitPos[i] = sigBases.get(i).getCoordinate();
            sigHitPlusStr[i] = sigBases.get(i).getStrand() == '+';
            sigHitCounts[i]=sigBases.get(i).getCount();
        }
            
        //Load control read info
        if(ctrlHitNum>0){
	        ctrlHitCounts= new double[ctrlHitNum];
	        ctrlHitPos= new int[ctrlHitNum];
	        ctrlHitPlusStr= new boolean[ctrlHitNum];
	        for(int i=0;i<ctrlHitNum;i++){
	            ctrlHitPos[i] = ctrlBases.get(i).getCoordinate();
	            ctrlHitPlusStr[i] = ctrlBases.get(i).getStrand() == '+';
	            ctrlHitCounts[i]=ctrlBases.get(i).getCount();
	        }
        }

        //Load pi for binding components
        for(int j=0;j<numComp;j++){
        	BindingSubComponents comp = components.get(j);
            pi[j]= comp.getPi(); 
        }
        //Load pi for noise components
        piNoise=noise.get(c).getPi();
            
        //Load binding component positions
        for(int j=0;j<numComp;j++)
            mu[j] = components.get(j).getPosition();
    		
        //Initialize responsibility functions
        double[][] hc= new double[numComp][sigHitNum];
        double[][] thc= new double[numComp][sigHitNum];
        double[] nc = new double[sigHitNum];
        for(int i=0;i<sigHitNum;i++){
        	for(int j=0;j<numComp;j++){
            	int dist = sigHitPlusStr[i] ? sigHitPos[i]-mu[j]: mu[j]-sigHitPos[i];
            		
                hc[j][i] = bindingModels[sigRepIndices[i]].probability(dist);
                thc[j][i] = bindingModels[sigRepIndices[i]].probability(dist);

            }
            nc[i] = noise.get(c).scorePosition(sigHitPos[i],sigRepIndices[i]);
        }
        h = hc;
        n = nc;
        tmp_h = thc;
            
        rBindSig  = new double[numComp][sigHitNum];
        rNoiseSig = new double[sigHitNum];
    	lastRBind = new double[numComp][sigHitNum];
    	tmp_rBindSig  = new double[numComp][sigHitNum];
    	tmp_rNoiseSig = new double[sigHitNum];
        
    	//End of data structure initialization
        
        
        //////////
        // Run ML steps
        //////////
        ML(w, cond);
        
        
        //////////
        // Assign ML result to BindingEvents
        //////////
        List<BindingEvent> events = new ArrayList<BindingEvent>();
	   	for(int j=0;j<numComp;j++){ 
	   		BindingEvent event = new BindingEvent(components.get(j).getCoord(), w);
	   		
    		ArrayList<Sample> controlsSeen = new ArrayList<Sample>();
    		boolean uniformBackAdded=false;
    		double condSigResp = 0.0, condCtrlResp = 0.0;
    		c = cond.getIndex();
    		for(ControlledExperiment rep : cond.getReplicates()){
    			int r = rep.getIndex();
	    		double repSigResp = 0.0, repCtrlResp = 0.0; 
	    		if(pi[j]>0){
	    			double scount=0;
			           for(int i=0;i<sigHitNum;i++)
			            if(sigRepIndices[i]==r)
			            	scount += sigHitCounts[i]*rBindSig[j][i];
			        repSigResp+=scount;
			        condSigResp+=scount;
			            
			        double ccount=0;
			        if(rep.hasControl()){
			        	for(int i=0;i<ctrlHitNum;i++)
			            	if(ctrlRepIndices[i]==r)
			            		ccount += ctrlHitCounts[i]*rBindCtrl[j][i];
			            repCtrlResp+=ccount;
				        if(!controlsSeen.contains(rep.getControl()))
				        	condCtrlResp+=ccount;
				        controlsSeen.add(rep.getControl());
			        }else{  //If there is no control channel, assign pseudo-control counts as if the noise reads in the IP channel were distributed perfectly uniformly
			        	repCtrlResp = uniformRepHitCountTotals[rep.getIndex()]*pi[j];
			            if(!uniformBackAdded)
			            	condCtrlResp+=repCtrlResp;
			            uniformBackAdded=true;
			        }
	    		}
		        event.setRepSigHits(rep, repSigResp);
		        event.setRepCtrlHits(rep, repCtrlResp);
    		}
    		event.setCondSigHits(cond, condSigResp);
	        event.setCondCtrlHits(cond, condCtrlResp);
	        if(evconfig.CALC_EVENTS_LL)
	        	event.setLLd(cond, compLL[c][j]);
	            
    		if(nonZeroInAny(event, cond))
    			events.add(event);
        }        
        return events;
    }//end of EMTrain method

    private boolean nonZeroInAny(BindingEvent ev, ExperimentCondition cond){
    	boolean nonZero=false;
    	for(ControlledExperiment rep : cond.getReplicates())
    		nonZero = nonZero || ev.getRepSigHits(rep)>=1;
    	return nonZero;
    }

    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void ML (Region currRegion, ExperimentCondition cond) {
    	int c = cond.getIndex();
        int numComp = numComponents;
        double [] totalRespCtrl = null;
        
        //Initialize responsibilities
        double [] totalRespSig = new double[sigHitNum];
    	for(int i=0;i<sigHitNum;i++)
            totalRespSig[i] = 0;
    	if(ctrlHitNum>0){
	    	rBindCtrl = new double[numComp][ctrlHitNum];
	    	rNoiseCtrl= new double[ctrlHitNum];
	    	totalRespCtrl = new double[ctrlHitNum];
	        for(int i=0;i<ctrlHitNum;i++)
	            totalRespCtrl[i] = 0;
    	}        
    	
    	////////////////////////////
        //Run ML -- this should only need one or two rounds
    	////////////////////////////
        for(int t=0; t<config.EM_ML_ITER ; t++){ //System.out.println(t); 
        	
    		////////
    		//E-step
    		////////
        	//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
        	for(int i=0;i<sigHitNum;i++)
                for(int j=0;j<numComp;j++){ if(pi[j]>0){
                    int dist = sigHitPlusStr[i] ? sigHitPos[i]-mu[j]: mu[j]-sigHitPos[i];
                    h[j][i] = bindingModels[sigRepIndices[i]].probability(dist);
                }}
        	//Compute responsibilities
    		for(int i=0;i<sigHitNum;i++)
                totalRespSig[i] = 0;
        	for(int i=0;i<sigHitNum;i++){
        		for(int j=0;j<numComp;j++){ if(pi[j]>0){
        			rBindSig[j][i] = h[j][i]*pi[j];
        			totalRespSig[i] +=rBindSig[j][i]; 
        		}}
        		rNoiseSig[i] = n[i] * piNoise;
        		totalRespSig[i] +=rNoiseSig[i];
        	}
        	//Normalize responsibilities
        	for(int i=0;i<sigHitNum;i++){
        		for(int j=0;j<numComp;j++){ if(pi[j]>0){
        			rBindSig[j][i]/=totalRespSig[i];
        		}}
        		rNoiseSig[i]/=totalRespSig[i];
        	}
    		        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
        	//Maximize pi
        	double[] sumR=new double[numComponents];
        	for(int j=0;j<numComp;j++){ if(pi[j]>0){
        		for(int i=0;i<sigHitNum;i++)
        			sumR[j] += rBindSig[j][i]*sigHitCounts[i];
            }}
                
        	// No components to be eliminated in ML, update pi(j)
        	for(int j=0;j<numComp;j++){ 
        		pi[j]=Math.max(0, sumR[j]); 
        	}
                
            //Normalize pi (accounting for piNoise)
            double totalPi=0;
            for(int j=0;j<numComp;j++){ if(pi[j]>0){
            	totalPi+=pi[j];
            }}
            for(int j=0;j<numComp;j++){ if(pi[j]>0){
            	if(totalPi>0)
            		pi[j]=pi[j]/(totalPi/(1-piNoise));
            }}

        	
        	//Non-zero components count
        	int nonZeroComps=0;
        	for(int j=0;j<numComp;j++)
        		if(pi[j]>0.0)
        			nonZeroComps++;
        	
        	////////////
        	//Check Stopping condition
        	////////////	
            if (nonZeroComps>0 && (t==0 || !lastEquivToCurr())){
            	copyStateToLast();
                continue;
            }else{
            	copyStateToLast();
            	break;
            }
        } //LOOP: Run ML while not converged
        //Base log-likelihood calculation
        double baseLL=0;
    	for(int i=0;i<sigHitNum;i++){
    		// for each read, each event will give a conditional prob or bg prob
            double j_sum=0;
    		for(int j=0;j<numComp;j++){ if(pi[j]>0.0){
    			j_sum += Math.log(rBindSig[j][i])/config.LOG2;
            }}
    		j_sum += Math.log(rNoiseSig[i])/config.LOG2;
            baseLL += j_sum*sigHitCounts[i];                        
        }

        
        //ML assignment of signal reads to components is finished
        //Assign control reads with converged pi values here
        double[][] hCtrl= new double[numComp][sigHitNum];
        double[] nCtrl = new double[sigHitNum];

		//Recompute h & n functions for control reads, given binding component positions 
		for(int i=0;i<sigHitNum;i++){
	        for(int j=0;j<numComp;j++){ if(pi[j]>0){
	            int dist = ctrlHitPlusStr[i] ? ctrlHitPos[i]-mu[j]: mu[j]-ctrlHitPos[i];
                hCtrl[j][i] = bindingModels[ctrlRepIndices[i]].probability(dist);
	        }}
	        nCtrl[i] = noise.get(c).scorePosition(ctrlHitPos[i], ctrlRepIndices[i]);
		}
		//Compute responsibilities
		for(int i=0;i<sigHitNum;i++)
	        totalRespCtrl[i] = 0;
		for(int i=0;i<sigHitNum;i++){
			for(int j=0;j<numComp;j++){ if(pi[j]>0){
				rBindCtrl[j][i] = hCtrl[j][i]*pi[j];
				totalRespCtrl[i] +=rBindCtrl[j][i]; 
			}}
			rNoiseCtrl[i] = nCtrl[i] * piNoise;
			totalRespCtrl[i] +=rNoiseCtrl[i];
		}
		//Normalize responsibilities
		for(int i=0;i<sigHitNum;i++){
			for(int j=0;j<numComp;j++){ if(pi[j]>0){
				rBindCtrl[j][i]/=totalRespCtrl[i];
			}}
			rNoiseCtrl[i]/=totalRespCtrl[i];
		}

        
        if(evconfig.CALC_EVENTS_LL){
	        /////////////////////////////////////////////////////////////////////////////////
	        //Finally, run ML on alternate models, where each component is eliminated in turn
	        /////////////////////////////////////////////////////////////////////////////////
	        for(int elim=0;elim<numComp;elim++){ 
	        	//If this component is active in any condition
	        	boolean active=false;
	        	if(pi[elim]>0){ active=true;}
	        	if(active){
	        		for(int t=0; t<=2 ; t++){  
	        			//Initialize temporary variables
	        			for(int j=0;j<numComp;j++)
	        				if(j==elim)
	        					tmp_pi[j]=0;
	        				else
	        					tmp_pi[j]= pi[j];
	        			if(piNoise+pi[elim]==1.0)//Single component case
	        				tmp_piNoise=1.0;
	        			else{
	        				tmp_piNoise = piNoise;
	        				double totalPi=0;
	                        for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                        	totalPi+=tmp_pi[j];
	                        }}
	                        for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                        	if(totalPi>0)
	                        		tmp_pi[j]=tmp_pi[j]/(totalPi/(1-tmp_piNoise));
	                        }}	
	        			}

	        			
	                	////////
	            		//E-step
	            		////////
	                	//Recompute h function, given binding component positions (n function is constant because noise model doesn't move)
	                	for(int i=0;i<sigHitNum;i++)
	                        for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                            int dist = sigHitPlusStr[i] ? sigHitPos[i]-mu[j]: mu[j]-sigHitPos[i];
	                            tmp_h[j][i] = bindingModels[sigRepIndices[i]].probability(dist);
	                        }}
	                	//Compute responsibilities
	            		for(int i=0;i<sigHitNum;i++)
	                        totalRespSig[i] = 0;
	                	for(int i=0;i<sigHitNum;i++){
	                		for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                			tmp_rBindSig[j][i] = tmp_h[j][i]*tmp_pi[j];
	                			totalRespSig[i] +=tmp_rBindSig[j][i]; 
	                		}}
	                		tmp_rNoiseSig[i] = n[i] * tmp_piNoise;
	                		totalRespSig[i] +=rNoiseSig[i];
	                	}
	                	//Normalize responsibilities
	                	for(int i=0;i<sigHitNum;i++){
	                		for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                			tmp_rBindSig[j][i]/=totalRespSig[i];
	                		}}
	                		tmp_rNoiseSig[i]/=totalRespSig[i];
	                	}
	            		        
	            		/////////////////////
	            		//M-step: maximize pi
	            		/////////////////////
	            		if(tmp_piNoise<1.0){//No need for maximization in single component cases
	                	//Maximize pi
	                	double[] sumR=new double[numComponents];
	                	for(int j=0;j<numComp;j++){ if(pi[j]>0){
	                		for(int i=0;i<sigHitNum;i++)
	                			sumR[j] += tmp_rBindSig[j][i]*sigHitCounts[i];
	                    }}
	                        
	                	// No components to be eliminated in ML, update pi(j)
	                	for(int j=0;j<numComp;j++){ 
	                		tmp_pi[j]=Math.max(0, sumR[j]); 
	                	}
	                        
	                    //Normalize pi (accounting for piNoise)
	                    double totalPi=0;
	                    for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                    	totalPi+=tmp_pi[j];
	                    }}
	                    for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0){
	                    	if(totalPi>0)
	                    		tmp_pi[j]=tmp_pi[j]/(totalPi/(1-tmp_piNoise));
	                    }}
	            		}
	        		}
	        		// log-likelihood calculation
	                compLL[c][elim] =-2*baseLL;
	            	int numBases = sigHitNum;
	            	for(int i=0;i<numBases;i++){
	            		// for each read, each event will give a conditional prob or bg prob
	                    double j_sum=0;
	            		for(int j=0;j<numComp;j++){ if(tmp_pi[j]>0.0){
	            			j_sum += Math.log(tmp_rBindSig[j][i])/config.LOG2;
	                    }}
	            		j_sum += Math.log(tmp_rNoiseSig[i])/config.LOG2;
	                    compLL[c][elim] += 2*j_sum*sigHitCounts[i];                        
	                }
	        	}
	        }
        }        
    }//end of ML method
 
    /**
     * Copy current variables to last variables (lastRBind, lastPi, lastMu).
     * Assumes visibility of both.
     */
    private void copyStateToLast(){
    	for(int j=0; j<numComponents; j++){
    		lastPi[j] = pi[j];
    		lastMu[j] = mu[j];
    		for(int x=0; x<rBindSig[j].length; x++){
    			lastRBind[j][x] = rBindSig[j][x];
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
    		for(int x=0; x<rBindSig[j].length; x++){
    			rBindEquivalent = rBindEquivalent && (Math.abs(rBindSig[j][x]-lastRBind[j][x])<config.EM_STATE_EQUIV_THRES);
    		}
		}}
		return numCompEqual && compPosEqual && piBindEquivalent && rBindEquivalent;
    }
}
