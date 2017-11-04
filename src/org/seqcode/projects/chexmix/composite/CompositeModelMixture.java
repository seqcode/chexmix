package org.seqcode.projects.chexmix.composite;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fitting.GaussianFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.viz.compositeplot.TagProfile;
import org.seqcode.viz.compositeplot.TagProfilePaintable;


/**
 * CompositeModelMixture: defines a protein-DNA interaction mixture model in a composite tag distribution.
 * 
 *  This mixture assumes that the same protein-DNA interaction model is valid across all examined conditions. 
 * 
 * TODO: support multi-rep TagDistributions
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class CompositeModelMixture {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig config;
	protected ExperimentCondition condition;
	protected CompositeTagDistribution compositeDistrib; //The composite tag distribution under investigation
	protected CompositeTagDistribution controlCompositeDistrib; //Optional composite tag distribution from matching control
	protected ProteinDNAInteractionModel model; //The model to train
	protected List<CompositeModelSiteAssignment> siteAssignments;
	protected CompositeModelEM EMtrainer; //Training method
	
	protected TagProbabilityDensity initBackDistrib, initCSDistrib, initXLDistrib;
	protected double initNoisePi;
	protected List<Region> regionsToPlot;
	protected int trainingRound=0;
	protected boolean trainedModel=false;
	
	public CompositeModelMixture(CompositeTagDistribution tagDist, CompositeTagDistribution ctrlDist, GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon, ExperimentCondition cond){
		gconfig = gcon;
		econfig = econ;
		config = ccon;
		condition = cond;
		compositeDistrib = tagDist;
		controlCompositeDistrib = ctrlDist;
		regionsToPlot = config.getRegionsToPlot();

		EMtrainer = new CompositeModelEM(compositeDistrib, config, condition);

		initializeTagDistributions();
		
		model = new ProteinDNAInteractionModel(config, compositeDistrib, initXLDistrib, initCSDistrib, initBackDistrib, initNoisePi);
		siteAssignments = new ArrayList<CompositeModelSiteAssignment>();
	}
	
	//Accessors
	public ProteinDNAInteractionModel getModel(){return model;}
	public void setModel(ProteinDNAInteractionModel mod){model=mod; trainedModel=true;} //Useful if not training model, but still want to use ML on a saved model
	public List<CompositeModelSiteAssignment> getSiteAssignments(){return siteAssignments;}
	
	
	/**
	 * Train the mixture model with EM  
	 * 
	 */
	public void trainEM(){		
		try{
			//Run EM outer loops
			boolean converged=false;
			double[] kl = new double[config.getMaxModelUpdateRounds()];
			for(trainingRound=0; trainingRound<config.getMaxModelUpdateRounds() && !converged; trainingRound++){
				//Run EM inner loops
				model = EMtrainer.train(model, trainingRound, true);
				kl[trainingRound]=updateTagDistributions();
			
				System.out.println("TrainingRound: "+trainingRound+":\n"+model.toString()+"\n\tModelKL="+kl[trainingRound]);
				
				//Check for convergence
				if(trainingRound>=config.getMinModelUpdateRounds()-1 && 
						kl[trainingRound]<config.getModelConvergenceKL())
					converged=true;
			}
			if(config.getPlotEM())
				EMtrainer.makeGifs();
			trainedModel = true;
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Assign component responsibilities to site tag distributions using a trained model. 
	 */
	public void assignML(boolean updatePi){
		try{
			if(trainedModel){
				//TODO: multithreading
				
				//Run ML
				CompositeModelMLAssign MLassigner = new CompositeModelMLAssign(compositeDistrib, model, config, condition);
				for(int s=0; s<compositeDistrib.getPoints().size(); s++){
					
					CompositeModelSiteAssignment currAssignment = MLassigner.assignToPointFromComposite(s, updatePi);
					siteAssignments.add(currAssignment);
				}
				
			}else{
				System.err.println("CompositeModelMixture Error: trying to use assignML with an untrained model...");
				System.exit(1);
			}
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	
     
	/**
	 * Initialize the CS, XL, and background tag distributions
	 */
	protected void initializeTagDistributions(){
		//CS (empirical)
		/*if(config.getDefaultCSModel()!=null)
			initCSDistrib = config.getDefaultCSModel();
		else
			initCSDistrib = new TagProbabilityDensity(TagProbabilityDensity.defaultChipSeqEmpiricalDistribution, null);
		*/
		initCSDistrib = new TagProbabilityDensity(compositeDistrib.getWinSize());
		initCSDistrib.loadGaussianDistrib(-150, 100, 150, 100);
		
		//XO
		initXLDistrib = new TagProbabilityDensity(200);
		initXLDistrib.loadGaussianDistrib(-config.getXLDistribOffset(), config.getXLDistribSigma(),config.getXLDistribOffset(), config.getXLDistribSigma());
		
		//Background
		if(controlCompositeDistrib==null){
			initBackDistrib = new TagProbabilityDensity(compositeDistrib.getWinSize());
			initBackDistrib.loadFlatDistrib();
			//estimate noise rate from the tails of the composite distribution
			double sum=0, count=0;
			double[] wprobs=compositeDistrib.getCompositeWatson(), cprobs=compositeDistrib.getCompositeCrick();
			for(int x=0; x<10; x++){
				sum+=wprobs[wprobs.length-x-1];
				sum+=cprobs[x];
				count+=2;
			}

			initNoisePi=Math.min(config.NOISE_EMISSION_MAX, Math.max(config.NOISE_EMISSION_MIN, (sum/count)*(compositeDistrib.getWinSize())));
		}else{
			//Make a paired list from the summed counts in the control tag distribution
			double[] wcounts = new double[controlCompositeDistrib.getWinSize()];
			double[] ccounts = new double[controlCompositeDistrib.getWinSize()];
			double sum=0;
			double[] wtmp = controlCompositeDistrib.getCompositeWatson();
			double[] ctmp = controlCompositeDistrib.getCompositeWatson();
			for(int x=0; x<controlCompositeDistrib.getWinSize(); x++){
				wcounts[x]+=wtmp[x];
				ccounts[x]+=ctmp[x];
				sum+=wtmp[x]+ctmp[x];
			}

			initBackDistrib = new TagProbabilityDensity(compositeDistrib.getWinSize());
			try {
				initBackDistrib.loadData(wcounts, ccounts);
			} catch (Exception e) {
				e.printStackTrace();
			}
			initNoisePi=sum/(controlCompositeDistrib.getWinSize()*2);
		}
	}
	
	/**
     * Update tag distributions given the current component responsibilities. 
     * 
     * @return log KL values
     */
    public Double updateTagDistributions(){
    	double logKL = 0;
    	
    	//CS component (empirical or gaussian fit)
    	TagProbabilityDensity CSdistrib = model.getCSTagDistribution();
    	double[] oldCSModelW = CSdistrib.getWatsonProbabilities();
    	double[] oldCSModelC = CSdistrib.getCrickProbabilities();
		double[] csW = model.getCSComponent().getTagProfile(true);
    	double[] csC = model.getCSComponent().getTagProfile(false);
    	//smooth with arbitrary gaussian
    	csW = StatUtil.gaussianSmoother(csW, 5);
    	csC = StatUtil.gaussianSmoother(csC, 5);
    	/*//This part is for an empirical distribution fit
    	double[] newCSModelW=new double[CSdistrib.getWinSize()];
		double[] newCSModelC=new double[CSdistrib.getWinSize()];
		for(int i=0; i<CSdistrib.getWinSize(); i++){
			newCSModelW[i]=0; newCSModelC[i]=0;
		}int i=0;
    	for(int x=offsetLeft; x<=offsetRight; x++){
    		newCSModelW[i]=csW[x];
    		newCSModelC[i]=csC[x];
    		i++;
    	}
    	StatUtil.mutate_normalize(newCSModelW);
    	StatUtil.mutate_normalize(newCSModelC);
    	try {
			CSdistrib.loadData(newCSModelW, newCSModelC);
		} catch (Exception e) {
			e.printStackTrace();
		}*/
    	GaussianFitter fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
    	for(int i=0; i<CSdistrib.getWinSize(); i++){
    		fitter.addObservedPoint((double)(i+CSdistrib.getLeft()), csW[i]+csC[CSdistrib.getWinSize()-i-1]);
    	}
    	double[] parameters = fitter.fit();;
    	double newOffset = parameters[1];
    	double newSigma = parameters[2];
    	System.out.println("CSGaussianFit:\t"+newOffset+"\t"+newSigma);
    	CSdistrib.loadGaussianDistrib(newOffset, newSigma, -newOffset, newSigma); //Symmetric CS component
    	double[] newCSModelW = CSdistrib.getWatsonProbabilities();
    	double[] newCSModelC = CSdistrib.getCrickProbabilities();
    	//Calc KL
    	logKL += StatUtil.log_KL_Divergence(oldCSModelW, newCSModelW);
    	logKL += StatUtil.log_KL_Divergence(oldCSModelC, newCSModelC);
    	
    	
    	//XL components (fit gaussian)
    	for(CompositeModelComponent xlComp : model.getXLComponents()){
    		if(xlComp.isNonZero()){
	    		TagProbabilityDensity XLdistrib = model.getXLTagDistribution(xlComp);
	    		if(!XLdistrib.isGaussian())
	    			System.err.println("CompositeModelMixture: Warning: replacing an empirical distribution with a Gaussian in XL component");
	    		
	    		double[] oldXLModelW = XLdistrib.getWatsonProbabilities();
		    	double[] oldXLModelC = XLdistrib.getCrickProbabilities();
				double[] xlW = xlComp.getTagProfile(true);
	    		double[] xlC = xlComp.getTagProfile(false);
	    		
	    		if(config.fixedXLOffset()){
	    			//Mean stays constant - calculate sigmas 
	    			double newSigmaW=config.getXLDistribSigma(), newSigmaC=config.getXLDistribSigma(),
	    					offsetW = XLdistrib.getGaussOffsetW(), offsetC = XLdistrib.getGaussOffsetC(),
	    					sumDiffW=0, sumDiffC=0, totW=0, totC=0;
	    			int l=XLdistrib.getLeft();
	    			
	    			if(config.XL_DISTRIB_SYMMETRIC){
	    				for(int i=0; i<XLdistrib.getWinSize(); i++){
		    				//Note that the mean (i.e. position) is at position zero w.r.t XL distribution 
		    				sumDiffW += (xlW[i]*(l+i-offsetW)*(l+i-offsetW)) + 
		    						    (xlC[XLdistrib.getWinSize()-i-1]*(l+i-offsetW)*(l+i-offsetW)); //note... -oldOffsetW is correct here
		    				totW+=xlW[i]+xlC[XLdistrib.getWinSize()-i-1];
		    			}
	    				sumDiffC = sumDiffW; totC = totW;
	    			}else{
		    			for(int i=0; i<XLdistrib.getWinSize(); i++){
		    				//Note that the mean (i.e. position) is at position zero w.r.t XL distribution 
		    				sumDiffW += (xlW[i]*(l+i-offsetW)*(l+i-offsetW));
		    				totW+=xlW[i];
		    				sumDiffC += (xlC[i]*(l+i-offsetC)*(l+i-offsetC));
		    				totC+=xlC[i];
		    			}
	    			}
	    			if(sumDiffW>0)
	    				newSigmaW = Math.sqrt(sumDiffW/totW);
	    			if(sumDiffC>0)
	    				newSigmaC = Math.sqrt(sumDiffC/totC);
	    			System.out.println("XLGaussianFit:\tOffset:"+config.getXLDistribOffset()+"\tSigmaW:"+newSigmaW+"\tSigmaC:"+newSigmaC);
			    	XLdistrib.loadGaussianDistrib(-config.getXLDistribOffset(), newSigmaW, config.getXLDistribOffset(), newSigmaC);
	    		}else{
	    			//Calculate new mean and sigma (symmetric)
			    	fitter = new GaussianFitter(new LevenbergMarquardtOptimizer());
			    	for(int i=0; i<XLdistrib.getWinSize(); i++){
			    		fitter.addObservedPoint((double)(i+XLdistrib.getLeft()), xlW[i]+xlC[XLdistrib.getWinSize()-i-1]);
			    	}
			    	parameters = fitter.fit();;
			    	newOffset = parameters[1];
			    	newSigma = parameters[2];
			    	System.out.println("XLGaussianFit:\tOffset:"+newOffset+"\tSigma:"+newSigma);
			    	XLdistrib.loadGaussianDistrib(newOffset, newSigma, -newOffset, newSigma); //symmetric
	    		}
		    	double[] newXLModelW = XLdistrib.getWatsonProbabilities();
		    	double[] newXLModelC = XLdistrib.getCrickProbabilities();
		    	//Calc KL
		    	logKL += StatUtil.log_KL_Divergence(oldXLModelW, newXLModelW);
		    	logKL += StatUtil.log_KL_Divergence(oldXLModelC, newXLModelC);
    		}
    	}
    	
		return logKL;
	}
    
    /**
	 * Print per-site component responsibilities to a file
	 * @param cond
	 * @param filename
	 */
	public void printPerSiteComponentResponsibilitiesToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			fout.write("#Point\tTotalTags");
			for(CompositeModelComponent nz : model.getNonZeroComponents()){
				int pos = nz.getPosition() - model.getCenterOffset();
				String compLabel = nz.getLabel()+","+pos;
				fout.write("\t"+compLabel);
			}fout.write("\n");
			
			for(CompositeModelSiteAssignment cmsa : siteAssignments){
				fout.write(cmsa.toString());
    			fout.write("\n");
    		}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

    /**
	 * Print component profiles to a file
	 * @param cond
	 * @param filename
	 */
	public void printComponentProfilesToFile(ExperimentCondition cond, String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			// Print all component responsibility profiles
			String out="";
			for (CompositeModelComponent comp :  model.getAllComponents()){
				if(comp.isNonZero() && comp.getTagProfile(true) !=null){
					double[] watson = comp.getTagProfile(true);
					double[] crick = comp.getTagProfile(false);
					int width = watson.length;
					out=out+"#TagProfile,"+comp.getPosition()+","+width+",\n";	
					out = out+"#TagProfileWatson,";
					for(int x=0; x<width; x++)
						out=out+watson[x]+",";
					out = out+"\n";
					out = out+"#TagProfileCrick,";
					for(int x=0; x<width; x++)
						out=out+crick[x]+",";
					out = out+"\n";
				}
			}
			fout.write(out);
			fout.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void saveCompositePlots(String filename){
		try {
			//double[] Ymax = new double[manager.getNumConditions()];
			double Ymax=0;
			//full composites
			TagProfile full = new TagProfile(compositeDistrib.getCompositeWatson(), 
					compositeDistrib.getCompositeCrick(), 
					compositeDistrib.getCenterOffset());
			TagProfilePaintable fullPainter = new TagProfilePaintable(full);
			fullPainter.autoYmax(true);
			fullPainter.setFilledColumns(true);
			fullPainter.setWatsonColor(new Color(122,88,143));
			fullPainter.setCrickColor(new Color(154,114,179));
			Ymax = Math.max(Ymax,  fullPainter.getYmax());
			String fullImageFileName = filename+"_compositeFull."+".png";
			fullPainter.saveImage(new File(fullImageFileName), 1060, 500, true);
			fullPainter.setProfileLeftLimit(-50);
			fullPainter.setProfileRightLimit(+49);
			String fullZoomImageFileName = filename+"_compositeFullZoom."+".png";
			fullPainter.saveImage(new File(fullZoomImageFileName), 1060, 500, true);

			
			//Background component responsibilities 
			TagProfile back = new TagProfile(model.getBackgroundComponent().getTagProfile(true), 
					model.getBackgroundComponent().getTagProfile(false), 
					compositeDistrib.getCenterOffset());
			TagProfilePaintable backPainter = new TagProfilePaintable(back);
			backPainter.autoYmax(false);
			backPainter.setYmax(Ymax);
			backPainter.setFilledColumns(true);
			backPainter.setWatsonColor(new Color(145,145,145));
			backPainter.setCrickColor(new Color(181,181,181));
			String backImageFileName = filename+"_backCompResponsibilities.png";
			backPainter.saveImage(new File(backImageFileName), 1060, 500, true);
			backPainter.setProfileLeftLimit(-50);
			backPainter.setProfileRightLimit(+49);
			String backZoomImageFileName = filename+"_backCompResponsibilitiesZoom.png";
			backPainter.saveImage(new File(backZoomImageFileName), 1060, 500, true);
			
			//CS component responsibilities 
			TagProfile cs = new TagProfile(model.getCSComponent().getTagProfile(true), 
					model.getCSComponent().getTagProfile(false), 
					compositeDistrib.getCenterOffset());
			TagProfilePaintable csPainter = new TagProfilePaintable(cs);
			csPainter.autoYmax(false);
			csPainter.setYmax(Ymax);
			csPainter.setFilledColumns(true);
			csPainter.setWatsonColor(new Color(179,59,61));
			csPainter.setCrickColor(new Color(218,87,88));
			String csImageFileName = filename+"_csCompResponsibilities.png";
			csPainter.saveImage(new File(csImageFileName), 1060, 500, true);
			csPainter.setProfileLeftLimit(-50);
			csPainter.setProfileRightLimit(+49);
			String csZoomImageFileName = filename+"_csCompResponsibilitiesZoom.png";
			csPainter.saveImage(new File(csZoomImageFileName), 1060, 500, true);
			
			//XL component responsibilities
			for(CompositeModelComponent xlComp : model.getXLComponents()){
				if(xlComp.isNonZero() && xlComp.getTagProfile(true) !=null){
					double[] rW = new double[compositeDistrib.getWinSize()];
		        	double[] rC = new double[compositeDistrib.getWinSize()];
		        	for(int r=0; r<compositeDistrib.getWinSize(); r++){ rW[r]=0; rC[r]=0;}
		        	double[] xlW=xlComp.getTagProfile(true);
		        	double[] xlC=xlComp.getTagProfile(false);
		        	int pos = xlComp.getPosition()-compositeDistrib.getCenterOffset();
		        	for(int i=xlComp.getTagDistribution().getLeft(); i<xlComp.getTagDistribution().getRight(); i++){ 
		        		int iNew = pos+i + compositeDistrib.getCenterOffset();
		        		int iOld = i-xlComp.getTagDistribution().getLeft();
		        		if(iNew>0 && iNew<compositeDistrib.getWinSize()){
		        			rW[iNew]= xlW[iOld];
		        			rC[iNew]= xlC[iOld];
		        		}
		        	}
		        	
					TagProfile currXL = new TagProfile(rW, rC,
							compositeDistrib.getCenterOffset());
					TagProfilePaintable xlPainter = new TagProfilePaintable(currXL);
					xlPainter.autoYmax(false);
					xlPainter.setYmax(Ymax);
					xlPainter.setFilledColumns(true);
					xlPainter.setWatsonColor(new Color(62,115,165));
					xlPainter.setCrickColor(new Color(88,147,204));
					xlPainter.setPointOfInterest(xlComp.getPosition()-compositeDistrib.getCenterOffset());
					String xlImageFileName = filename+"_XL"+pos+"Responsibilities.png";
					xlPainter.saveImage(new File(xlImageFileName), 1060, 500, true);
					xlPainter.setProfileLeftLimit(-50);
					xlPainter.setProfileRightLimit(+49);
					String xlZoomImageFileName  = filename+"_XL"+pos+"ResponsibilitiesZoom.png";
					xlPainter.saveImage(new File(xlZoomImageFileName), 1060, 500, true);
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}
    
}
