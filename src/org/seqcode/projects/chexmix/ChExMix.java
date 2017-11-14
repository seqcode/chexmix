package org.seqcode.projects.chexmix;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.math.diff.Normalization;
import org.seqcode.projects.chexmix.composite.CompositeModelComponent;
import org.seqcode.projects.chexmix.composite.CompositeModelMixture;
import org.seqcode.projects.chexmix.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.ProteinDNAInteractionModel;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.composite.XLAnalysisConfig;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingModel;
import org.seqcode.projects.chexmix.events.EnrichmentSignificance;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.framework.OutputFormatter;
import org.seqcode.projects.chexmix.framework.PotentialRegionFilter;
import org.seqcode.projects.chexmix.mixturemodel.BindingMixture;
import org.seqcode.projects.chexmix.utilities.EventsPostAnalysis;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;

public class ChExMix {
	
	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected ChExMixConfig gpsconfig;
	protected XLAnalysisConfig xlconfig;
	protected BindingManager bindingManager;
	protected CompositeModelMixture modelmix;
	protected BindingMixture mixtureModel;
	protected PotentialRegionFilter potentialFilter;
	protected OutputFormatter outFormatter;
	protected Normalization normalizer;
	protected Map<ControlledExperiment, List<TagProbabilityDensity>> repBindingModels;
	protected Map<ControlledExperiment, List<BindingModel>> repUnstrandedBindingModels;
	
	public ChExMix(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, ChExMixConfig c, XLAnalysisConfig xc, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		manager = eMan;
		gpsconfig = c;
		xlconfig = xc;
		xlconfig.makeXLAnalysisOutputDirs(true);
		gpsconfig.makeChExMixOutputDirs(true);
		outFormatter = new OutputFormatter(gpsconfig);
		bindingManager = new BindingManager(evconfig, manager);
		
		//Initialize binding models & binding model record
		repBindingModels = new HashMap<ControlledExperiment, List<TagProbabilityDensity>>();
		repUnstrandedBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
		for(ControlledExperiment rep : manager.getReplicates())
			repBindingModels.put(rep,  new ArrayList<TagProbabilityDensity>());
		
		List<TagProbabilityDensity> tagProbDensities = new ArrayList<TagProbabilityDensity>();	
		boolean strandedModelSet =false;
		if (!gpsconfig.getInitialClustPoints().isEmpty()){
			List<List<StrandedPoint>> initialClustPoints = gpsconfig.getInitialClustPoints();
			for (List<StrandedPoint> points : initialClustPoints){
				List<StrandedPoint> compositePoints = new ArrayList<StrandedPoint>();
				compositePoints.addAll(points);
				//Build the composite distribution(s)
				CompositeTagDistribution signalComposite = new CompositeTagDistribution(compositePoints, manager.getConditions().get(0), gpsconfig.MAX_BINDINGMODEL_WIDTH,true);
				TagProbabilityDensity currDensity=new TagProbabilityDensity(signalComposite.getWinSize()-1);

				try {
					currDensity.loadData(signalComposite.getCompositeWatson(), signalComposite.getCompositeCrick());
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				tagProbDensities.add(currDensity);
			}	
			strandedModelSet=true;
		}else if (xlconfig.getModelFilename()!= null){
			ProteinDNAInteractionModel model = ProteinDNAInteractionModel.loadFromFile(xlconfig, new File(xlconfig.getModelFilename()));
			List<ProteinDNAInteractionModel> models = new ArrayList<ProteinDNAInteractionModel>();
			models.add(model);
			TagProbabilityDensity density = model.makeTagProbabilityDensityFromAllComponents();
			tagProbDensities.add(density);
			for(ControlledExperiment rep : manager.getReplicates())
				bindingManager.setProteinDNAInteractionModel(rep, models);
			strandedModelSet=true;
		}		
		// Set tag probability density models to binding manager
		if (strandedModelSet){
			for(ControlledExperiment rep : manager.getReplicates()){
				bindingManager.setBindingModel(rep, tagProbDensities);
				repBindingModels.get(rep).addAll(bindingManager.getBindingModel(rep));	
			}
		}
		
		// Set unstranded binding models
		for(ControlledExperiment rep : manager.getReplicates()){		
			if(evconfig.getDefaultBindingModel()!=null){
				bindingManager.setUnstrandedBindingModel(rep, evconfig.getDefaultBindingModel());
				repUnstrandedBindingModels.put(rep, new ArrayList<BindingModel>());
				repUnstrandedBindingModels.get(rep).add(bindingManager.getUnstrandedBindingModel(rep));
			}else{
				bindingManager.setUnstrandedBindingModel(rep, new BindingModel(BindingModel.defaultChipExoEmpiricalDistribution));
				repUnstrandedBindingModels.put(rep, new ArrayList<BindingModel>());
				repUnstrandedBindingModels.get(rep).add(bindingManager.getUnstrandedBindingModel(rep));
			}
		}	
		
		// Testing only
		if (strandedModelSet){printInitialDistribution();}
		List<Integer> motifIndexes = new ArrayList<Integer>();
		motifIndexes.add(-1);
		for(ExperimentCondition cond : manager.getConditions()){
			bindingManager.updateMaxInfluenceRange(cond, true);
			bindingManager.setMotifIndexes(cond, motifIndexes);
		}
		
		//Find potential binding regions
		System.err.println("Finding potential binding regions.");
		potentialFilter = new PotentialRegionFilter(evconfig, gpsconfig, econfig, manager, bindingManager);
		
		
		List<Region> potentials = null;
		if (gpsconfig.getInitialPos()!=null)
			potentials = potentialFilter.executeInitialPositionFiler();
		else
			potentials = potentialFilter.execute();
		System.err.println(potentials.size()+" potential regions found. Total length: "+potentialFilter.getPotRegionLengthTotal());
		if(potentials.size()==0){
			System.err.println("No potential regions - exiting.");
			System.exit(1);
		}
		potentialFilter.printPotentialRegionsToFile();		
	}
	
	/**
	 * Make initial tag probability density to do binding event finding
	 */
	public TagProbabilityDensity makeInitialTagProbabilityDensity(){
		
		int modelWidth = xlconfig.getCompositeWinSize();		
		// Initial density distribution for each component
		TagProbabilityDensity initBackDistrib, initCSDistrib, initXLDistrib;
		//CS
		initCSDistrib = new TagProbabilityDensity(modelWidth);
//		initCSDistrib.loadGaussianDistrib(-150, 100, 150, 100);
		initCSDistrib.loadGaussianDistrib(-100, 75, 100, 75);
		
		//XO
		initXLDistrib = new TagProbabilityDensity(200);
		initXLDistrib.loadGaussianDistrib(-xlconfig.getXLDistribOffset(), xlconfig.getXLDistribSigma(),xlconfig.getXLDistribOffset(), xlconfig.getXLDistribSigma());
		//Background
		initBackDistrib = new TagProbabilityDensity(modelWidth);
		initBackDistrib.loadFlatDistrib();
		
		// Make CompositeModelComponents
		//Background component
		CompositeModelComponent backgroundComponent = new CompositeModelComponent(initBackDistrib, modelWidth/2, 0, "Back",  false, true);
		//ChIP-seq component
		CompositeModelComponent CSComponent = new CompositeModelComponent(initCSDistrib, modelWidth/2, 1, "CS", false, true);
		//XL components
		CompositeModelComponent XLComponent = new CompositeModelComponent(initXLDistrib, modelWidth/2, 3, "XL", true, true);
		
		//All components
		List<CompositeModelComponent> allComponents = new ArrayList<CompositeModelComponent>();
		allComponents.add(backgroundComponent);
		allComponents.add(CSComponent);
		allComponents.add(XLComponent);
		
		// Setting pi(s)
		double noisePi = gpsconfig.NOISE_EMISSION_MIN;
		backgroundComponent.setPi(noisePi);
		double bindingPi = 1-noisePi;
		double initCS = 0.95; // Arbitrary choice for initial CS pi
		CSComponent.setPi(initCS);
		XLComponent.setPi(bindingPi*(1-initCS));
		
		//Sum all components to create TagProbabilityDensity
		double[] watsons = new double [modelWidth];
		double[] cricks = new double [modelWidth];
		for (int i=0; i < modelWidth; i++){
			watsons[i]=0;
			cricks[i]=0;
		}
		
		for (CompositeModelComponent comp : allComponents){
			if (comp.getPi() >0){
				double[] currCompWatsonProb = comp.getTagDistribution().getWatsonProbabilities();
				double[] currCompCrickProb = comp.getTagDistribution().getCrickProbabilities();
				for (int i=0; i < currCompWatsonProb.length ; i++){
					int index = comp.getPosition()-comp.getTagDistribution().getWinSize()/2+i;
					if (index >=0 && index < modelWidth){
						watsons[index] += currCompWatsonProb[i]*comp.getPi();
						cricks[index] += currCompCrickProb[i]*comp.getPi();
					}					
				}
			}
		}
			
		List<Pair<Integer,Double>> empiricalWatson = new ArrayList<Pair<Integer,Double>>(); 
		List<Pair<Integer,Double>> empiricalCrick = new ArrayList<Pair<Integer,Double>>(); 
		for (int i=0; i<modelWidth; i++){			
			empiricalWatson.add(new Pair<Integer,Double>(i-modelWidth/2,watsons[i]));
			empiricalCrick.add(new Pair<Integer,Double>(i-modelWidth/2,cricks[i]));
		}
		return (new TagProbabilityDensity(empiricalWatson,empiricalCrick));
	}
	
	/**
	 * Print initial density to file : testing only
	 */
	public void printInitialDistribution(){
		// write tag probability density of initial distribution
		String distribFilename = gpsconfig.getOutputIntermediateDir()+File.separator+gpsconfig.getOutBase()+"_t"+0+"_ReadDistrib_CompositeComponents";
		for(ControlledExperiment rep : manager.getReplicates()){
			int i=0;
			for (TagProbabilityDensity probs : bindingManager.getBindingModel(rep)){
				probs.printDensityToFile(distribFilename+"_"+i);
				i++;
			}
		}
	}
	

	/**
	 * Run the mixture model to find binding events. 
	 * @throws Exception 
	 */
	public void runMixtureModel() throws Exception {
		
		Double[] kl;
		System.err.println("Initialzing mixture model");
		mixtureModel = new BindingMixture(gconfig, econfig, evconfig, gpsconfig,xlconfig, manager, bindingManager, potentialFilter);
		
		int round = 0;
		boolean converged = false;
		
		System.err.println("\n============================ Round "+round+" ============================");
		//Execute the MultiGPS mixture model
		mixtureModel.execute(true, true, false); //EM
		
		//Update noise models
        mixtureModel.updateGlobalNoise();
		
        //Print current components
        mixtureModel.printActiveComponentsToFile();
        
        //Do ML assignment and enrichment analysis after the first round of multiGPS style peak calling
        //ML quantification of events
        mixtureModel.execute(false, false, true); //ML
        bindingManager.setBindingEvents(mixtureModel.getBindingEvents());
        //Update sig & noise counts in each replicate
        bindingManager.estimateSignalVsNoiseFractions(bindingManager.getBindingEvents());
        //Statistical analysis: Enrichment over controls 
        EnrichmentSignificance testerR0 = new EnrichmentSignificance(evconfig, manager, bindingManager, evconfig.getMultiGPSMinEventFoldChange(), econfig.getMappableGenomeLength());
        testerR0.execute();
        
        mixtureModel.setActiveComponents(bindingManager.getComponentsFromEnrichedEvents());
        		
        while (!converged){
        	
        	//Update motifs
            mixtureModel.updateMotifs();
            
            //Update binding models
            String distribFilename = gpsconfig.getOutputIntermediateDir()+File.separator+gpsconfig.getOutBase()+"_t"+round;
            kl = mixtureModel.updateBindingModelUsingMotifs(distribFilename);
//               kl = mixtureModel4.updateBindingModelUsingClustering(distribFilename);
            
            if (!gpsconfig.getInitialClustPoints().isEmpty() && round ==0){
            	// Use provided initial cluster points
            	System.out.println("round "+round + "use provided read density from clusters");
            }else{
            	kl = mixtureModel.updateBindingModelUsingReadDistributions(distribFilename);
            }

            //Add new binding models to the record
            for(ControlledExperiment rep : manager.getReplicates())           	
    			repBindingModels.get(rep).addAll(bindingManager.getBindingModel(rep));
            mixtureModel.updateAlphas();	
        	
            System.err.println("\n============================ Round "+round+" ============================");
            
            //Execute the mixture model
            mixtureModel.execute(true, false, false); //EM
            
            //Update noise models
            mixtureModel.updateGlobalNoise();
            
            //Print current components
            mixtureModel.printActiveComponentsToFile();          
            
            round++;
            
            //Check for convergence
            if(round>gpsconfig.getMaxModelUpdateRounds()){
            	converged=true;
            }else{
            	converged = true;
            	for(int l=0; l<kl.length; l++)
            		converged = converged && (kl[l]<-5 || kl[l].isNaN());
            }
        }
       
        outFormatter.plotAllReadDistributions(repBindingModels);
        
        //ML quantification of events
        System.err.println("\n============================ ML read assignment ============================");
        mixtureModel.execute(false, false, false); //ML        
        bindingManager.setBindingEvents(mixtureModel.getBindingEvents());
        bindingManager.updateNumBindingTypes();
        //Update sig & noise counts in each replicate
        bindingManager.estimateSignalVsNoiseFractions(bindingManager.getBindingEvents());
        System.err.println("ML read assignment finished.");
                
        System.err.println("\n============================= Post-processing ==============================");
        
        //Align motifs to get relative offsets
        if(gpsconfig.getFindingMotifs())
        	mixtureModel.getMotifFinder().alignMotifs();
        
        //Statistical analysis: Enrichment over controls 
        EnrichmentSignificance tester = new EnrichmentSignificance(evconfig, manager, bindingManager, evconfig.getMinEventFoldChange(), econfig.getMappableGenomeLength());
		tester.execute();
        
		//Write the replicate counts to a file (needed before EdgeR differential enrichment)
		bindingManager.writeReplicateCounts(gpsconfig.getOutputParentDir()+File.separator+gpsconfig.getOutBase()+".replicates.counts");
        
        // Print final events to files
		bindingManager.writeBindingEventFiles(gpsconfig.getOutputParentDir()+File.separator+gpsconfig.getOutBase(), evconfig.getQMinThres(), evconfig.getRunDiffTests(), evconfig.getDiffPMinThres());
        System.err.println("Binding event detection finished!\nBinding events are printed to files in "+gpsconfig.getOutputParentDir()+" beginning with: "+gpsconfig.getOutName());
        
        
        //Post-analysis of peaks
        EventsPostAnalysis postAnalyzer = new EventsPostAnalysis(gconfig,evconfig, gpsconfig, manager, bindingManager, bindingManager.getBindingEvents(), mixtureModel.getMotifFinder());
        postAnalyzer.execute(400);
        
    }
	
	
	/**
	 * Main driver method for xogps
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception{
		System.setProperty("java.awt.headless", "true");
		System.err.println("ChExMix version "+ChExMixConfig.version+"\n\n");
		GenomeConfig gcon = new GenomeConfig(args);
		EventsConfig evconfig = new EventsConfig(gcon, args);
		ChExMixConfig config = new ChExMixConfig(gcon, args);
		XLAnalysisConfig xlconfig = new XLAnalysisConfig(gcon, args);
		
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if (!config.useReadFilter())
			econ.setPerBaseReadFiltering(false);	
		econ.setLoadRead2(false);
		
		if(config.helpWanted()){
			System.err.println(config.getArgsList());
		}else{
			
			ExperimentManager manager = new ExperimentManager(econ);
			
			//Just a test to see if we've loaded all conditions
			if(manager.getConditions().size()==0){
				System.err.println("No experiments specified. Use --expt or --design options."); System.exit(1);
			}
			
			ChExMix gps = new ChExMix(gcon, econ, evconfig, config, xlconfig, manager);
			gps.runMixtureModel();
			
			manager.close();
		}
	}
	
	/**
	 * returns a string describing the arguments for the public version of ChExMix. 
	 * @return String
	 */
	public static String getxogpsArgsList(){
		return(new String("" +
				"Copyright (C) Shaun Mahony 2017\n" +
				"<http://mahonylab.org/software/chexmix>\n" +
				"\n" +
				"ChExMix comes with ABSOLUTELY NO WARRANTY.  This is free software, and you\n"+
				"are welcome to redistribute it under certain conditions.  See the MIT license \n"+
				"for details.\n"+
				"\n OPTIONS:\n" +
				" General:\n"+
				"\t--out <output file prefix>\n" +
				"\t--threads <number of threads to use (default=1)>\n" +
				"\t--verbose [flag to print intermediate files and extra output]\n" +
				"\t--config <config file: all options here can be specified in a name<space>value text file, over-ridden by command-line args>\n" +
				" Genome:\n" +
				"\t--geninfo <genome info file> AND --seq <fasta seq directory reqd if using motif prior>\n" +
				" Loading Data:\n" +
				"\t--expt <file name> AND --format <SAM/BED/IDX>\n" +
				"\t--ctrl <file name (optional argument. must be same format as expt files)>\n" +
				"\t--design <experiment design file name to use instead of --expt and --ctrl; see website for format>\n"+
				"\t--fixedpb <fixed per base limit (default: estimated from background model)>\n" +
				"\t--poissongausspb <filter per base using a Poisson threshold parameterized by a local Gaussian sliding window>\n" +
				"\t--nonunique [flag to use non-unique reads]\n" +
				"\t--mappability <fraction of the genome that is mappable for these experiments (default=0.8)>\n" +
				"\t--nocache [flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)]\n" +
				"Scaling control vs signal counts:\n" +
				"\t--noscaling [flag to turn off auto estimation of signal vs control scaling factor]\n" +
				"\t--medianscale [flag to use scaling by median ratio (default = scaling by NCIS)]\n" +
				"\t--regressionscale [flag to use scaling by regression (default = scaling by NCIS)]\n" +
				"\t--sesscale [flag to use scaling by SES (default = scaling by NCIS)]\n" +
				"\t--fixedscaling <multiply control counts by total tag count ratio and then by this factor (default: NCIS)>\n" +
				"\t--scalewin <window size for scaling procedure (default=10000)>\n" +
				"\t--plotscaling [flag to plot diagnostic information for the chosen scaling method]\n" +
				" Running MultiGPS:\n" +
				"\t--d <binding event read distribution file>\n" +
				"\t--r <max. model update rounds, default=3>\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update (default=500)>\n" +
				"\t--nomodelsmoothing [flag to turn off binding model smoothing]\n" +
				"\t--splinesmoothparam <spline smoothing parameter (default=30)>\n" +
				"\t--gaussmodelsmoothing [flag to turn on Gaussian model smoothing (default = cubic spline)]\n" +
				"\t--gausssmoothparam <Gaussian smoothing std dev (default=3)>\n" +
				"\t--jointinmodel [flag to allow joint events in model updates (default=do not)]\n" +
				"\t--fixedmodelrange [flag to keep binding model range fixed to inital size (default: vary automatically)]\n" +
				"\t--prlogconf <Poisson log threshold for potential region scanning(default=-6)>\n" +
				"\t--alphascale <alpha scaling factor(default=1.0>\n" +
				"\t--fixedalpha <impose this alpha (default: set automatically)>\n" +
				"\t--mlconfignotshared [flag to not share component configs in the ML step]\n" +
				"\t--exclude <file of regions to ignore>\n" +
				" MultiGPS priors:\n"+
				"\t--noposprior [flag to turn off inter-experiment positional prior (default=on)]\n" +
				"\t--probshared <probability that events are shared across conditions (default=0.9)>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default=6)>\n"+
				"\t--mememaxw <maxw arg for MEME (default=18)>\n"+
				"\t--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n"+
				" Reporting binding events:\n" +
				"\t--q <Q-value minimum (default=0.001)>\n" +
				"\t--minfold <minimum event fold-change vs scaled control (default=1.5)>\n" +
				"\t--nodifftests [flag to turn off differential enrichment tests]\n" +
				"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
				"\t--edgerod <EdgeR overdispersion parameter (default=0.15)>\n" +
				"\t--diffp <minimum p-value for reporting differential enrichment (default=0.01)>\n" +
				""));
	}

}
