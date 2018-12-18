package org.seqcode.projects.chexmix;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.math.diff.Normalization;

import org.seqcode.projects.chexmix.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.ProteinDNAInteractionModel;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.composite.XLAnalysisConfig;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingModel;
import org.seqcode.projects.chexmix.events.BindingSubtype;
import org.seqcode.projects.chexmix.events.EnrichmentSignificance;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.framework.OutputFormatter;
import org.seqcode.projects.chexmix.framework.PotentialRegionFilter;
import org.seqcode.projects.chexmix.mixturemodel.BindingMixture;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
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
	protected ChExMixConfig chexconfig;
	protected XLAnalysisConfig xlconfig;
	protected ShapeAlignConfig shapeconfig;
	protected BindingManager bindingManager;
	protected BindingMixture mixtureModel;
	protected PotentialRegionFilter potentialFilter;
	protected OutputFormatter outFormatter;
	protected Normalization normalizer;
	protected Map<ControlledExperiment, List<TagProbabilityDensity>> repBindingModels;
	protected Map<ControlledExperiment, List<BindingModel>> repUnstrandedBindingModels;
	protected Map<ExperimentCondition, List<BindingSubtype>> prevModels;
	
	public ChExMix(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, ChExMixConfig c, XLAnalysisConfig xc, ShapeAlignConfig sc, ExperimentManager eMan){
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		manager = eMan;
		chexconfig = c;
		xlconfig = xc;
		xlconfig.makeXLAnalysisOutputDirs(true);
		shapeconfig = sc;
		chexconfig.makeChExMixOutputDirs(true);
		outFormatter = new OutputFormatter(chexconfig);
		bindingManager = new BindingManager(evconfig, manager, chexconfig.getReplicateConsistencyMode());
		
		//Initialize binding models & binding model record
		repBindingModels = new HashMap<ControlledExperiment, List<TagProbabilityDensity>>();
		repUnstrandedBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();		
		for(ControlledExperiment rep : manager.getReplicates())
			repBindingModels.put(rep,  new ArrayList<TagProbabilityDensity>());
		prevModels = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		
		List<TagProbabilityDensity> tagProbDensities = new ArrayList<TagProbabilityDensity>();	
		if (chexconfig.getInitialClustPoints()!=null){
			List<List<StrandedPoint>> initialClustPoints = chexconfig.getInitialClustPoints();
			for (List<StrandedPoint> points : initialClustPoints){
				List<StrandedPoint> compositePoints = new ArrayList<StrandedPoint>();
				compositePoints.addAll(points);
				//Build the composite distribution(s)
				CompositeTagDistribution signalComposite = new CompositeTagDistribution(compositePoints, manager.getConditions().get(0), chexconfig.MAX_BINDINGMODEL_WIDTH,true);
				TagProbabilityDensity currDensity=new TagProbabilityDensity(signalComposite.getWinSize()-1);

				try {
					currDensity.loadData(signalComposite.getCompositeWatson(), signalComposite.getCompositeCrick());
				} catch (Exception e) {
					e.printStackTrace();
				}
				tagProbDensities.add(currDensity);
			}	
		}else if (xlconfig.getModelFilename()!= null){
			ProteinDNAInteractionModel model = ProteinDNAInteractionModel.loadFromFile(xlconfig, new File(xlconfig.getModelFilename()));
			tagProbDensities.add(model.makeTagProbabilityDensityFromAllComponents());
		}else if (chexconfig.getDistribA()!=null){
			File pFile = new File(chexconfig.getDistribA());
			if(!pFile.isFile()){
				System.err.println("\nCannot find read distribution file: "+chexconfig.getDistribA());
				System.exit(1);
			}
			tagProbDensities.add(new TagProbabilityDensity(pFile));
			
			if (chexconfig.getDistribB()!=null){
				File pbFile = new File(chexconfig.getDistribB());
				if(!pFile.isFile()){
					System.err.println("\nCannot find read distribution file: "+chexconfig.getDistribB());
					System.exit(1);
				}
				tagProbDensities.add(new TagProbabilityDensity(pbFile));				
			}		
		}else{	//make default model using default unstranded model
			List<Pair<Integer,Double>> watsonModel= null;
			if(evconfig.getDefaultBindingModel()!=null)
				watsonModel=evconfig.getDefaultBindingModel().getEmpiricalDistribution();
			else
				watsonModel= BindingModel.defaultChipExoEmpiricalDistribution;
			List<Pair<Integer,Double>> crickModel = new ArrayList<Pair<Integer,Double>>();
			List<Integer> pos = new ArrayList<Integer>();
			List<Double> prob = new ArrayList<Double>();
			for (Pair<Integer,Double> data : watsonModel){
				pos.add(data.car());
				prob.add(data.cdr());
			}
			List<Double> revProb = new ArrayList<Double>(prob);
			Collections.reverse(revProb);
			for (int i=0; i < revProb.size(); i++)
				crickModel.add(new Pair<Integer,Double>(pos.get(i), revProb.get(i)));
			tagProbDensities.add(new TagProbabilityDensity(watsonModel, crickModel));
		}	
		// Set tag probability density models to binding manager
		for(ExperimentCondition cond : manager.getConditions()){
			List<BindingSubtype> initSubtypes = new ArrayList<BindingSubtype>();
			for (TagProbabilityDensity currDensity : tagProbDensities)
				initSubtypes.add(new BindingSubtype(cond, currDensity, 0));
			bindingManager.setBindingSubtypes(cond, initSubtypes);
		}
		
		// Set unstranded binding models
		for(ControlledExperiment rep : manager.getReplicates()){		
			if(evconfig.getDefaultBindingModel()!=null)
				bindingManager.setUnstrandedBindingModel(rep, evconfig.getDefaultBindingModel());
			else
				bindingManager.setUnstrandedBindingModel(rep, new BindingModel(BindingModel.defaultChipExoEmpiricalDistribution));
			repUnstrandedBindingModels.put(rep, new ArrayList<BindingModel>());
			repUnstrandedBindingModels.get(rep).add(bindingManager.getUnstrandedBindingModel(rep));
		}	
		
		// Testing only
		printInitialDistribution();
		
		//Find potential binding regions
		System.err.println("Finding potential binding regions.");
		potentialFilter = new PotentialRegionFilter(evconfig, chexconfig, econfig, manager, bindingManager);
		
		
		List<Region> potentials = null;
		if (chexconfig.getInitialPos()!=null)
			potentials = potentialFilter.executeInitialPositionFiler();
		else
			potentials = potentialFilter.execute();
		System.err.println(potentials.size()+" potential regions found. Total length: "+potentialFilter.getPotRegionLengthTotal());
		if(potentials.size()==0){
			System.err.println("No potential regions - exiting.");
	//		System.exit(1);
		}
		potentialFilter.printPotentialRegionsToFile();		
	}
	
	/**
	 * Print initial density to file : testing only
	 */
	public void printInitialDistribution(){
		// write tag probability density of initial distribution
		String distribFilename = chexconfig.getOutputIntermediateDir()+File.separator+chexconfig.getOutBase()+"_t"+0+"_initReadDistrib";
		for (ExperimentCondition cond : manager.getConditions()){
			int i=0;
			for (BindingSubtype sub : bindingManager.getBindingSubtype(cond)){
				sub.getBindingModel(0).printDensityToFile(distribFilename+"_"+i+".txt");
				i++;
			}
		}
	}
	
	/**
	 * Accessor for testing
	 */	
	public BindingManager getBindingManager(){return bindingManager;}
	public PotentialRegionFilter getPotentialFilter(){return potentialFilter;}
	
	

	/**
	 * Run the mixture model to find binding events. 
	 * @throws Exception 
	 */
	public void runMixtureModel() throws Exception {
		
		System.err.println("Initializing mixture model");
		mixtureModel = new BindingMixture(gconfig, econfig, evconfig, chexconfig,xlconfig,shapeconfig, manager, bindingManager, potentialFilter);
		
		int round = 0;
		boolean converged = false;
		
		System.err.println("\n============================ Round "+round+" ============================");
		//Execute the ChExMix mixture model
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
        testerR0.execute(chexconfig.getModelRange());
        
        mixtureModel.setActiveComponents(bindingManager.getComponentsFromEnrichedEvents(potentialFilter.getPotentialRegions()));
        
        round++;
        		
        while (!converged){
        	
        	//Update motifs
            mixtureModel.updateBindingModelUsingMotifs();
            
            //Update binding models
            String distribFilename = chexconfig.getOutputIntermediateDir()+File.separator+chexconfig.getOutBase()+"_t"+round;
            if (round <= 1){
                mixtureModel.doReadDistributionClustering();
            }else{
            	mixtureModel.updateBindingModelUsingReadDistributions(distribFilename);
            }
             
            // Merge similar binding models
            mixtureModel.consolidateBindingModels(); 
            
            // Update alpha
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
            if(round>chexconfig.getMaxModelUpdateRounds())
            	converged=true;

        }  
        
        //Add binding models to the record
        for(ControlledExperiment rep : manager.getReplicates())
        	for (BindingSubtype sub : bindingManager.getBindingSubtype(rep.getCondition()))
        		repBindingModels.get(rep).add(sub.getBindingModel(0));
       
        outFormatter.plotAllReadDistributions(repBindingModels);
        
        //ML quantification of events
        System.err.println("\n============================ ML read assignment ============================");
        mixtureModel.clearBindingEvents(); // Clear binding events found in initial multiGPS style peak calling
        mixtureModel.execute(false, false, false); //ML        
        bindingManager.setBindingEvents(mixtureModel.getBindingEvents());
        bindingManager.updateNumBindingTypes();
        //Update sig & noise counts in each replicate
        bindingManager.estimateSignalVsNoiseFractions(bindingManager.getBindingEvents());
        System.err.println("ML read assignment finished.");
                
        System.err.println("\n============================= Post-processing ==============================");
        
        //Align motifs to get relative offsets
        if(chexconfig.getFindingMotifs())
        	mixtureModel.getMotifFinder().alignMotifs();
        
        //Statistical analysis: enrichment over controls 
        EnrichmentSignificance tester = new EnrichmentSignificance(evconfig, manager, bindingManager, evconfig.getMinEventFoldChange(), econfig.getMappableGenomeLength());
		tester.execute();
        
		//Write the replicate counts to a file (needed before EdgeR differential enrichment)
		bindingManager.writeReplicateCounts(chexconfig.getOutputParentDir()+File.separator+chexconfig.getOutBase()+".replicates.counts");
        
		// Print per-replicate events to files (replicate consistency mode only)
		if(chexconfig.getReplicateConsistencyMode())
			bindingManager.writePerReplicateBindingEventFiles(chexconfig.getOutputParentDir()+File.separator+chexconfig.getOutBase(), evconfig.getQMinThres(), evconfig.getRunDiffTests(), evconfig.getDiffPMinThres());
        
		// Print final events to files
		bindingManager.writeBindingEventFiles(chexconfig.getOutputParentDir()+File.separator+chexconfig.getOutBase(), evconfig.getQMinThres(), evconfig.getRunDiffTests(), evconfig.getDiffPMinThres());
        System.err.println("Binding event detection finished!\nBinding events are printed to files in "+chexconfig.getOutputParentDir()+" beginning with: "+chexconfig.getOutName());
        
        
        //Post-analysis of peaks
        EventsPostAnalysis postAnalyzer = new EventsPostAnalysis(gconfig,evconfig, chexconfig, manager, bindingManager, bindingManager.getBindingEvents(), mixtureModel.getMotifFinder());
        postAnalyzer.execute(400);
        
    }
	
	
	/**
	 * Main driver method for ChExMix
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
		ShapeAlignConfig sc = new ShapeAlignConfig(args);
		sc.setWindowSize(config.getClusterWindowSize());
		
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if (!config.useReadFilter())
			econ.setPerBaseReadFiltering(false);	
		econ.setLoadRead2(false);
		
		if(config.helpWanted()){
			System.err.println(ChExMix.getChExMixArgsList());
		}else{
			
			ExperimentManager manager = new ExperimentManager(econ);
			
			//Just a test to see if we've loaded all conditions
			if(manager.getConditions().size()==0){
				System.err.println("No experiments specified. Use --expt or --design options."); System.exit(1);
			}
			
			ChExMix gps = new ChExMix(gcon, econ, evconfig, config, xlconfig,sc, manager);
			gps.runMixtureModel();
			
			manager.close();
		}
	}
	
	/**
	 * returns a string describing the arguments for the public version of ChExMix. 
	 * @return String
	 */
	public static String getChExMixArgsList(){
		return(new String("" +
				"Copyright (C) Naomi Yamada 2018\n" +
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
				"\t--geninfo <genome info file> AND --seq <fasta seq directory reqd if finding motif>\n" +
				"\t--back <Markov background model file reqd if finding motif>\n"+
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
				" Running ChExMix:\n" +
				"\t--round <max. model update rounds (default=3)>\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update (default=100)>\n" +
				"\t--prlogconf <Poisson log threshold for potential region scanning (default=-6)>\n" +
				"\t--alphascale <alpha scaling factor (default=1.0)>\n" +
				"\t--betascale <beta scaling factor (default=0.05)>\n" +
				"\t--epsilonscale <epsilon scaling factor (default=0.2)>\n" +
				"\t--fixedalpha <impose this alpha (default: set automatically)>\n" +
				"\t--mlconfignotshared [flag to not share component configs in the ML step]\n" +
				"\t--exclude <file of regions to ignore> OR --excludebed <file of regions to ignore in bed format>\n" +
				"\t--peakf <file of peaks to initialize component positions>\n" +
				"\t--motfile <file of motifs in transfac format to initialize subtype motifs>\n" +
				"\t--norepcon [flag to turn off replicate consistency mode (reproduces binding event reporting behavior in v0.2 and below)]\n" +				
				"\t--galaxyhtml [flag to produce a html output appropreate for galaxy]\n" +
				" Finding ChExMix subtypes using motif:\n"+
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default=6)>\n"+
				"\t--mememaxw <maxw arg for MEME (default=18)>\n"+
				"\t--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n"+
				"\t--minroc <minimum motif ROC value (default=0.7)>\n"+
				"\t--minmodelupdaterefs <minimum number of motif reference to support an subtype distribution update (default=50)>\n"+
				"\t--seqrmthres <Filter out sequences with motifs below this threshold for recursively finding motifs (default=0.1)>\n" +
				" Finding ChExMix subtypes using read distribution clustering:\n"+
				"\t--noclustering [flag to turn off read distribution clustering]\n" +
				"\t--pref <preference value for read distribution clustering (default=-0.1)>\n"+
				"\t--numcomps <number of components to cluster (default=500)>\n"+
				"\t--win <window size of read profiles (default=150)>\n"+
				" Reporting binding events:\n" +
				"\t--q <Q-value minimum (default=0.01)>\n" +
				"\t--minfold <minimum event fold-change vs scaled control (default=1.5)>\n" +
				""));
	}

}
