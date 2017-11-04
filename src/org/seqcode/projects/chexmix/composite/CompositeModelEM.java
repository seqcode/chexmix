package org.seqcode.projects.chexmix.composite;


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

import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.utilities.TrainingStepPlotter;
import org.seqcode.viz.utils.GifSequenceWriter;


/**
 * CompositeModelEM: run EM training with sparse prior on a composite tag distribution.
 * 
 * This training method assumes that the same protein-DNA interaction model is valid (i.e. constant) across all examined conditions. 
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CompositeModelEM {

	protected ExperimentCondition condition;
	protected ChExMixConfig config;
	protected CompositeTagDistribution composite;
	protected ProteinDNAInteractionModel model;
	protected int numComponents;  //The count of all components (active +inactive) in the model
	protected int trainingRound=0; //Identifier for the overall training round, used only for file names
	//	EM VARIABLES
	// H function and responsibility have to account for all reads in region now, as they will be updated 
    // once the component positions change (i.e. we can't do the trick where we restrict to reads within 
    // range of the components).
	protected double[]   hitCounts;		// Hit weights
	protected int[]      hitPos;		// Hit positions
	protected boolean[]  hitPlusStr;	// Hit positive strand boolean
	protected int		   hitNum;		// Number of hits in each condition 
	protected double[][] hAll;			// H function values for all positions in the current window (precomputed)
	protected double[][] h; 			// H function (binding component probability per read)
	protected double[][] r;				// Binding component responsibilities
	protected double[]   pi;			// pi : emission probabilities for binding components
	protected int[]      mu;			// mu : positions of the binding components
	protected double    alphaMax;		// Maximum alpha
	protected double[][] lastRBind;		//Last responsibilities (monitor convergence)
	protected double[]   lastPi;		//Last Pi (monitor convergence)
	protected int[]      lastMu;		//Last positions (monitor convergence)
	protected double lastLAP, LAP; 		//log-likelihood monitoring
	protected boolean plotEM=false;		//Plot the current region components
	protected TrainingStepPlotter plotter = null;
	protected List<BufferedImage> fullImages, zoomImages;
	protected int stateEquivCount=0;
	
	/**
	 * Constructor
	 * 
	 * @param c: XLAnalysisConfig
	 * @param eMan: ExperimentManager
	 * @throws Exception 
	 */
	public CompositeModelEM(CompositeTagDistribution composite, ChExMixConfig c, ExperimentCondition cond){
		this.composite=composite;
		config=c;
		condition = cond;
		//Plotter?
    	plotEM = config.getPlotEM();
        if(plotEM){
        	plotter = new TrainingStepPlotter();
        	fullImages = new ArrayList<BufferedImage>();
        	zoomImages = new ArrayList<BufferedImage>();
        }
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
    public ProteinDNAInteractionModel  train(ProteinDNAInteractionModel model,
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
    	pi = new double[numComponents];	// pi : emission probabilities for binding components
    	mu = new int[numComponents]; // mu : positions of the binding components (fixed across conditions)
        //Monitor state convergence using the following last variables
        lastPi = new double[numComponents];
        lastMu = new int[numComponents];
        
        //Initializing data structures
        for(int j=0;j<numComponents;j++){
        	//Load pi for binding components
        	CompositeModelComponent comp = model.getAllComponents().get(j);
            pi[j]= comp.getPi();
            //Load binding component positions
        	mu[j] = model.getAllComponents().get(j).getPosition();
        }
        
		//Set maximum alpha
    	alphaMax =  config.getFixedAlpha()>0 ? config.getFixedAlpha() :
    			Math.max(config.MIN_ALPHA, (config.getAlphaScalingFactor() * model.getBackgroundComponent().getPi())/composite.getWinSize()); 
    	System.out.println("\n\tT="+trainingRound+", Alpha= "+alphaMax);

    	//Condition-specific stuff
        int c = condition.getIndex();
        	
        //Number of unique positions in the composite distribution (watson + crick)
        hitNum=composite.getWinSize()*2;

        //Load read info
        double[] countc= new double[hitNum];
        int[] posc= new int[hitNum];
        boolean[] plusc= new boolean[hitNum];
        for(int i=0;i<composite.getWinSize();i++){ //Watson
            posc[i] = i;
            plusc[i] = true;
            countc[i]= composite.getCompositeWatson()[i];
        }
        for(int i=0;i<composite.getWinSize();i++){ //Crick
        	posc[i+composite.getWinSize()] = i;
            plusc[i+composite.getWinSize()] = false;
            countc[i+composite.getWinSize()]=composite.getCompositeCrick()[i];
        }
        hitPos = posc;
        hitCounts=countc;
        hitPlusStr = plusc;
        	
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

    	//Binding Components
    	for(int j=0;j<numComponents;j++){ 
            CompositeModelComponent comp = model.getAllComponents().get(j);
            comp.setPi(pi[j]);
            comp.setPosition(mu[j]);
            double sumRespW=0.0, sumRespC=0.0;	
            for(int i=0;i<hitNum;i++){
            	if(hitPlusStr[i])
            		sumRespW += hitCounts[i]*r[j][i];
            	else
            		sumRespC += hitCounts[i]*r[j][i];
            }
            comp.setSumResponsibilities(sumRespW, sumRespC);
    	}
    	//Responsibility profiles
        setComponentResponsibilityProfiles(r);
        
        //Print the responsibilities to files
        if(config.getPrintCompositeResponsibilities()){
        	String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()
        		+"_responsibilities."+condition.getName()+"T"+trainingRound+".txt";
        	printResponsibilitiesToFile(filename);
        }
        	
        return model;
    }//end of EMTrain method


    /**
     * Core EM iterations with sparse prior (component elimination) & multi-condition positional priors.
     * Assumes H function, pi, and responsibilities have all been initialized
     */
    private void EM_MAP () {
        int numComp = numComponents;
        int[] newMu = new int[numComponents];// mu update
        int csCompIndex = model.getCSComponent().getIndex(); 
        
        //Initialize responsibilities
    	int numBases = hitNum;
    	double [] totalResp = new double[hitNum];
        for(int i=0;i<hitNum;i++)
            totalResp[i] = 0;

        //Alpha is annealed in. Alpha=0 during ML steps
        double currAlpha = 0;
        

    	////////////
        //Plot the initial pi & priors if plotting
    	////////////
        if(plotEM && plotter!=null){
        	String outName = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EM_r"+trainingRound+"_t0" : null;
        	String outNameZoom = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EMzoom_r"+trainingRound+"_t0" : null;
        	int trimLeft = (model.getWidth()/2)-50;
        	int trimRight = (model.getWidth()/2)-50;
//        	fullImages.add(plotter.plotCompositeEM(outName, composite, model, mu, pi, trainingRound, 0, 0,0));
//        	zoomImages.add(plotter.plotCompositeEM(outNameZoom, composite, model, mu, pi, trainingRound, 0, trimLeft, trimRight));
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
    		//M-step: maximize mu (positions)
    		/////////////////////
			//Maximize mu: calculate maximization sums assuming no events shared across conditions
    		for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePosition()){
				//Define window to look for new component position
				int start=Math.max(mu[j]-config.EM_MU_UPDATE_WIN, 0);
    			int end = Math.min(model.getWidth(), mu[j]+config.EM_MU_UPDATE_WIN);

    			//Score the current window
    			double currScore=0, maxScore=-Double.MAX_VALUE;
    			int maxPos = 0;
    			
    			for(int x=start; x<end; x++){
    				currScore=0;
        			for(int i=0;i<hitNum;i++){
        				int dist = hitPos[i]-x;
        				currScore+=(r[j][i]*hitCounts[i]) * model.getAllComponents().get(j).getTagDistribution().logProbability(dist, hitPlusStr[i]);
        			}

    				
    				if(currScore>maxScore){
    					maxPos=x;
    					maxScore=currScore;
    				}
    			}
    			newMu[j] = maxPos;
    		}}
			
    		//Update mu values
			for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePosition()){
				mu[j] = newMu[j];
			}}
			

    		//Maximize mu follow-up: Resolve duplicate positions (combine & delete one copy)
    		HashMap<Integer, Integer> pos2index = new HashMap<Integer, Integer>(); //Position to array index map 
    		for(int j=0;j<numComp;j++){ if(pi[j]>0){
    			if(model.getAllComponents().get(j).hasUpdatablePosition()){
	    			if(pos2index.containsKey(mu[j])){ 
	    				int orig = pos2index.get(mu[j]);
	    				//Combine
	    				pi[orig]+=pi[j];
	                   	for(int i=0; i<hitNum;i++)
	                   		r[orig][i] += r[j][i];
	                   	//Delete
	                   	pi[j]=0.0;
	                   	for(int i=0; i<hitNum;i++)
	                   		r[j][i] = 0;

	    			}else{
	    				pos2index.put(mu[j], j);
	    			}
    			}
    		}}
    		
        		
    		/////////////////////
    		//M-step: maximize pi
    		/////////////////////
    		boolean componentEliminated=false;
    		double[] sumR=new double[numComponents];
    		for(int j=0;j<numComp;j++){ if(pi[j]>0){
    			for(int i=0;i<hitNum;i++)
    				sumR[j] += r[j][i]*hitCounts[i];
            }}
    		int minIndex=0; double minVal=Double.MAX_VALUE;
    		for(int j=0;j<numComp;j++){ if(pi[j]>0){
    			if(sumR[j]<minVal){ minVal=sumR[j]; minIndex=j;}
    		}}                
            if(minVal>currAlpha){
                // No component to be eliminated, update pi(j)
            	for(int j=0;j<numComp;j++){ if(pi[j]>0 && model.getAllComponents().get(j).hasUpdatablePi()){
            		pi[j]=Math.max(0, sumR[j]-currAlpha); 
            	}}
            }else{
                // Eliminate worst binding component
                // Responsibilities will be redistributed in the E step
               	pi[minIndex]=0.0; sumR[minIndex]=0.0;
               	for(int i=0; i<hitNum;i++)
               		r[minIndex][i] = 0;
               	//I discussed this bit with Chris, and we decided that the best thing to do is
               	//to re-estimate pi values for non-eliminated components using the current responsibility assignments
               	for(int j=0;j<numComp;j++){ 
               		if(j!=minIndex)
               			pi[j]=Math.max(0, sumR[j]); 
            	}
               	componentEliminated=true;
            }
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
            //Special case for CS component:
			if(pi[csCompIndex]<config.MIN_CS_PI){
				pi[csCompIndex]=config.MIN_CS_PI;
				//Renormalize
				double totPi=0;
				for(int j=0;j<numComp;j++){ if(pi[j]>0 && j!=csCompIndex && model.getAllComponents().get(j).hasUpdatablePi()){ 
					totPi +=  pi[j];
				}}
				for(int j=0;j<numComp;j++){ if(pi[j]>0 && j!=csCompIndex && model.getAllComponents().get(j).hasUpdatablePi()){
	        			pi[j]/=totPi;
	        			pi[j]*=1-totalNonUpdateablePi-config.MIN_CS_PI;
	        	}}
			}
        	
        	/////////////
        	//Anneal alpha
        	//////////////
    		if (t >config.EM_ML_ITER && t <= config.ALPHA_ANNEALING_ITER)
    			currAlpha = alphaMax * (t-config.EM_ML_ITER)/(config.ALPHA_ANNEALING_ITER-config.EM_ML_ITER);
    		else if(t > config.ALPHA_ANNEALING_ITER)
    			currAlpha = alphaMax;
            
        	
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
	            	
	            LP+=-(currAlpha*sum_log_pi);

	            LAP = LL+LP;
	
	            System.out.println("EM: "+t+"\t"+LAP+"\t("+nonZeroComps+" non-zero components).");
            }
        	
        	////////////
            //Plot the current pi & priors if plotting
        	////////////
        	if(plotEM && plotter!=null){
            	String outName = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EM_r"+trainingRound+"_t"+(t+1) : null;
            	String outNameZoom = config.getWriteSinglePlots() ? config.getOutputImagesDir()+File.separator+"EMzoom_r"+trainingRound+"_t"+(t+1) : null;
            	int trimLeft = (model.getWidth()/2)-50;
            	int trimRight = (model.getWidth()/2)-50;
//            	fullImages.add(plotter.plotCompositeEM(outName, composite, model, mu, pi, trainingRound, t, 0,0));
//            	zoomImages.add(plotter.plotCompositeEM(outNameZoom, composite, model, mu, pi, trainingRound, t, trimLeft, trimRight));
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
            if (nonZeroComps>0 && ( t<=config.ALPHA_ANNEALING_ITER || (config.CALC_LL && Math.abs(LAP-lastLAP)>config.EM_CONVERGENCE) || stateEquivCount<config.EM_STATE_EQUIV_ROUNDS)){
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
        
        //Initialize responsibilities
    	double [] totalResp = new double[hitNum];
        for(int i=0;i<hitNum;i++)
            totalResp[i] = 0;
        
    	
    	//////////////////////////////////////////////////////////////////////////////////////////
        //Run E step once
    	//////////////////////////////////////////////////////////////////////////////////////////
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
    }//end of responsibilityAssignment method
    
    /**
     * Set responsibility profile for each component (for kernel update)
     * @param bindComponents
     * @param signals
     * @param responsibilities
     */
    private void setComponentResponsibilityProfiles(double[][] responsibilities) {
		for(int j=0;j<numComponents;j++){
        	CompositeModelComponent comp = model.getAllComponents().get(j);
        	int width = comp.getTagDistribution().getWinSize();
        	int center = -1*comp.getTagDistribution().getLeft();
        	
        	// store binding profile (read responsibilities) of this component
			double[] profile_plus = new double[width];
			double[] profile_minus = new double[width]; 
			for (int i=0; i< width ; i++){
				profile_plus[i]=0; profile_minus[i]=0;
			}
			
			for(int i=0;i<hitNum;i++){
				int offset = hitPos[i]-comp.getPosition()+center;
				if(offset>=0 && offset<width)
					if (hitPlusStr[i])
						profile_plus[offset]=+responsibilities[j][i]*hitCounts[i];
					else
						profile_minus[offset]=+responsibilities[j][i]*hitCounts[i];
			}
        	comp.setTagProfile(profile_plus,  true);
			comp.setTagProfile(profile_minus, false);
		}
	}//end of setComponentResponsibilities method
	
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
    
	/**
	 * Print responsibilities for each base in the composite to the file
	 * @param cond
	 * @param filename
	 */
	private void printResponsibilitiesToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			fout.write("#Pos\tTotal");
			for(int j=0;j<numComponents;j++){ if(pi[j]>0){
				fout.write("\t"+j);
			}}fout.write("\n");
			
			for(int i=0;i<hitNum;i++){
				fout.write(hitPos[i]+"\t"+hitCounts[i]);
    			for(int j=0;j<numComponents;j++){ if(pi[j]>0){
    				double resp = hitPlusStr[i] ? r[j][i]*hitCounts[i] : -r[j][i]*hitCounts[i]; 
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
