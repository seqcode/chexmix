package org.seqcode.projects.chexmix.mixturemodel;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.deepseq.stats.PoissonBackgroundModel;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.gseutils.Pair;
import org.seqcode.gseutils.RealValuedHistogram;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.ml.clustering.Clusterable;
import org.seqcode.ml.clustering.affinitypropagation.APCluster;
import org.seqcode.ml.clustering.affinitypropagation.MatrixSimilarityMeasure;
import org.seqcode.ml.clustering.affinitypropagation.SimilarityMeasure;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.composite.XLAnalysisConfig;
import org.seqcode.projects.chexmix.events.BindingEvent;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingSubtype;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.framework.PotentialRegionFilter;
import org.seqcode.projects.chexmix.motifs.MotifPlatform;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignmentTesting;


/**
 * BindingMixture: defines a mixture model over binding events and subtypes.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class BindingMixture {

	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected EventsConfig evconfig;
	protected ChExMixConfig config;
	protected XLAnalysisConfig mixconfig;
	protected ShapeAlignConfig shapeconfig;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected PotentialRegionFilter potRegFilter;
	protected List<Region> testRegions;
	protected HashMap<Region, List<List<BindingSubComponents>>> activeComponents; //Components active after a round of execute()
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds=new HashMap<ExperimentCondition, BackgroundCollection>(); //Genomic Background models for each condition -- used to set alpha values in sparse prior
	protected HashMap<ControlledExperiment, BackgroundCollection> replicateBackgrounds=new HashMap<ControlledExperiment, BackgroundCollection>(); //Genomic Background models for each condition -- used to set alpha values in sparse prior in lenient mode
	protected List<BindingEvent> bindingEvents;
	protected List<Region> regionsToPlot;
	protected int trainingRound=0;
	protected double noisePerBase[];      //Defines global noise 
	protected double relativeCtrlNoise[]; //Defines global noise
	protected HashMap<Region, Double[]> noiseResp = new HashMap<Region, Double[]>(); //noise responsibilities after a round of execute(). Hashed by Region, indexed by condition
	protected MotifPlatform motifFinder;
	
	public BindingMixture(GenomeConfig gcon, ExptConfig econ, EventsConfig evcon, ChExMixConfig c,XLAnalysisConfig mixcon,ShapeAlignConfig sc, ExperimentManager eMan, BindingManager bMan, PotentialRegionFilter filter){
		gconfig = gcon;
		econfig = econ;
		evconfig = evcon;
		config = c;
		mixconfig = mixcon;
		shapeconfig = sc;
		manager = eMan;
		bindingManager = bMan;
		potRegFilter=filter;
		testRegions = filter.getPotentialRegions();
		List<Region> expTestRegions = new ArrayList<Region>();
		for(Region r : testRegions)
			expTestRegions.add(r.expand(evcon.SEQPLOTWIN*2, evcon.SEQPLOTWIN*2));
				
		if(config.getFindingMotifs())
			motifFinder = new MotifPlatform(gconfig, config, manager, bindingManager, expTestRegions);
		else 
			motifFinder=null;
		regionsToPlot = config.getRegionsToPlot();
		bindingEvents = new ArrayList<BindingEvent>();
		BindingEvent.setExperimentManager(manager);
		BindingEvent.setConfig(evconfig);
		BindingEvent.setGPSConfig(config);		
		
		activeComponents = new HashMap<Region, List<List<BindingSubComponents>>>();
		
		//Setting initial alpha value
		for(ExperimentCondition cond : manager.getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
			conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalCount()*(1-cond.getTotalSignalVsNoiseFrac()), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), config.getModelRange(), '.', 1, true));
			double condAlf = config.getFixedAlpha()>0 ? config.getFixedAlpha() : (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
			double alf=condAlf;
			
			//In lenient mode, we have to allow alpha to be small enough to catch events significant in individual replicates
			double repAlf = Double.MAX_VALUE;	
			if(config.getFixedAlpha()<=0 && (config.isLenientMode() || config.isLenientPlusMode())){
				for(ControlledExperiment rep : cond.getReplicates()){
					replicateBackgrounds.put(rep, new BackgroundCollection());
					replicateBackgrounds.get(rep).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), rep.getSignal().getHitCount()*(1-rep.getSignalVsNoiseFraction()), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), config.getModelRange(), '.', 1, true));
					repAlf = Math.min(repAlf, (double)replicateBackgrounds.get(rep).getMaxThreshold('.'));
				}
				alf = Math.min(condAlf,  repAlf);
			}
			
			System.err.println("Alpha "+cond.getName()+"\tRange="+config.getModelRange()+"\t"+alf);
			bindingManager.setAlpha(cond,alf);
		}
		
		noisePerBase = new double[manager.getNumConditions()];
		relativeCtrlNoise = new double[manager.getNumConditions()];
		initializeGlobalNoise();
		
	}
	
	
	/**
	 * Initialize the threads and execute the binding mixture
	 * Relies on the following being set: testRegions, bindingModels (from exptMan)
	 * uniformBindingSubComponentss: if true, components are initialized at uniform prob and spacing, 
	 * 					otherwise, components are initialized from activeComponents
	 * EM: if true, run EM, otherwise run ML assignment   
	 * bindingEvents : if true, run EM over binding events, otherwise run EM over cross-linking points
	 */
	public void execute(boolean EM, boolean uniformBindingSubComponents, boolean multiGPSML){
		if (EM)
			trainingRound++;
		
		//Have to split the test regions up by chromosome in order to maintain compatibility with experiment file cache loading
		//There will be some performance hit here, as all threads have to finish in a given chromosome before moving on to the next one. 
		Iterator<Region> chroms = new ChromosomeGenerator().execute(gconfig.getGenome());
		while (chroms.hasNext()) {
			Region currChr = chroms.next();
			List<Region> currChrTestReg = new ArrayList<Region>();
			for(Region r : testRegions)
				if(currChr.overlaps(r))
					currChrTestReg.add(r);
			
			if(currChrTestReg.size()>0){
				int numThreads = config.getMaxThreads()>currChrTestReg.size() ?  currChrTestReg.size() : config.getMaxThreads(); 
				Thread[] threads = new Thread[numThreads];
		        ArrayList<Region> threadRegions[] = new ArrayList[numThreads];
		        int i = 0;
		        for (i = 0 ; i < threads.length; i++) {
		            threadRegions[i] = new ArrayList<Region>();
		        }i=0;
		        for(Region r : currChrTestReg){
		            threadRegions[(i++) % numThreads].add(r);
		        }
		
		        for (i = 0 ; i < threads.length; i++) {
		            Thread t = new Thread(new BindingMixtureThread(threadRegions[i], EM, uniformBindingSubComponents, multiGPSML));
		            t.start();
		            threads[i] = t;
		        }
		        boolean anyrunning = true;
		        while (anyrunning) {
		            anyrunning = false;
		            try {
		                Thread.sleep(5000);
		            } catch (InterruptedException e) { }
		            for (i = 0; i < threads.length; i++) {
		                if (threads[i].isAlive()) {
		                    anyrunning = true;
		                    break;
		                }
		            }
		        }
			}
		}
	}
	
	/**
	 * Return the active components, flattening the hash map first. 
	 * @return
	 */
	public List<List<BindingSubComponents>> getBindingSubComponents(){
		List<List<BindingSubComponents>> comps = new ArrayList<List<BindingSubComponents>>();
		for(int c=0; c<manager.getNumConditions(); c++)
			comps.add(new ArrayList<BindingSubComponents>());
		for(Region r : activeComponents.keySet()){
			for(int c=0; c<manager.getNumConditions(); c++){
				comps.get(c).addAll(activeComponents.get(r).get(c));
			}
		}
		return comps;
	}
	
	/**
	 * Clear binding events.
	 * @return
	 */
	public void clearBindingEvents(){bindingEvents = new ArrayList<BindingEvent>();}
	
	/**
	 * Return the discovered binding events. Call after ML assignment.
	 * @return
	 */
	public List<BindingEvent> getBindingEvents(){return bindingEvents;}

	
	/**
	 * Return the initialized motif-finder
	 * @return
	 */
	public MotifPlatform getMotifFinder(){return motifFinder;}

	
	/**
	 * Set activeComponents after an initial enrichment test
	 */	
	public void setActiveComponents(HashMap<Region, List<List<BindingSubComponents>>> activeComponents){this.activeComponents = activeComponents;}
	
	/**
	 * Return minimum score for KL divergence comparison between two distribution with sliding window
	 */		
	public double getMinKLDivergenceScore(TagProbabilityDensity densityA, TagProbabilityDensity densityB){
		// Measure KL divergence with a sliding window
		double minLogKL = Double.MAX_VALUE;
		int minOffset= 0;
		boolean minReverse=false;
		for (int offset=-config.SLIDING_WINDOW/2; offset<=config.SLIDING_WINDOW/2; offset++){
			double currLogKL=0;
			//copy current window
			double[] currW = new double[config.MAX_BINDINGMODEL_WIDTH];
			double[] currC = new double[config.MAX_BINDINGMODEL_WIDTH];
			for (int w=0; w< config.MAX_BINDINGMODEL_WIDTH; w++){
				currW[w]=0; currC[w]=0;
			}
			for (int w=0; w< config.MAX_BINDINGMODEL_WIDTH; w++){
				if ((offset+w) >=0 && (offset+w) < config.MAX_BINDINGMODEL_WIDTH){
					currW[w]=densityB.getWatsonProbabilities()[offset+w];
					currC[w]=densityB.getCrickProbabilities()[offset+w];
				}
			}
			//Calc KL        		        				
			currLogKL += StatUtil.log_KL_Divergence(densityA.getWatsonProbabilities(), currW) + StatUtil.log_KL_Divergence(currW, densityA.getWatsonProbabilities());
			currLogKL += StatUtil.log_KL_Divergence(densityA.getCrickProbabilities(), currC) + StatUtil.log_KL_Divergence(currC, densityA.getCrickProbabilities()); 
	
			if (currLogKL < minLogKL) {
				minLogKL=currLogKL;
				minOffset=offset;
				minReverse=false;
			}
        				
        	// calculate divergence in reversed tags
        	currLogKL = 0;		
        	//copy current window
        	double[] rcurrW = new double[config.MAX_BINDINGMODEL_WIDTH];
        	double[] rcurrC = new double[config.MAX_BINDINGMODEL_WIDTH];
        	for (int w=0; w< config.MAX_BINDINGMODEL_WIDTH; w++){
        		rcurrW[w]=currC[config.MAX_BINDINGMODEL_WIDTH-w-1];
        		rcurrC[w]=currW[config.MAX_BINDINGMODEL_WIDTH-w-1];
        	}							       							
        	//Calc KL        		        				
        	currLogKL += StatUtil.log_KL_Divergence(densityA.getWatsonProbabilities(), rcurrW) + StatUtil.log_KL_Divergence(rcurrW, densityA.getWatsonProbabilities());
        	currLogKL += StatUtil.log_KL_Divergence(densityA.getCrickProbabilities(), rcurrC) + StatUtil.log_KL_Divergence(rcurrC, densityA.getCrickProbabilities());  

        	if (currLogKL < minLogKL) { 
        		minLogKL=currLogKL;
        		minOffset=offset;
        		minReverse=true;
        	}
		}
		return minLogKL;
	}
	
	
	public void doReadDistributionClustering() throws Exception{
		// Execute affinity propagation clustering
		if (config.getClusteringReads()){
		
			List<List<StrandedPoint>> initClustPoints = new ArrayList<List<StrandedPoint>>();		
			if (config.getInitialClustPoints()!=null){		
				initClustPoints = config.getInitialClustPoints();		
			}else{	
				
				System.out.println("Performing read distribution clustering");				
				// Choose which ones to include
				//Sum read profiles if there are enough binding components
				List<BindingSubComponents> currComps = new ArrayList<BindingSubComponents>();
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			double currAlpha = calcAlpha(cond);
	 
	    			//Choose which components to include
	    			for(Region r : activeComponents.keySet()){
	    				for(BindingSubComponents bc : activeComponents.get(r).get(cond.getIndex())){
	    					//1) Component must not be at the edge of the region 
	    					if((bc.getPosition()-r.getStart()>config.getModelRange()/2) && (r.getEnd()-bc.getPosition()>config.getModelRange()/2)){
	    						//2) Arbitrary minimum read support for BM components
	    						if(bc.getSumResponsibility()>(config.getMinComponentReadFactorForBM()*currAlpha))
	    							currComps.add(bc);
	    					}
	    				}
	    			}	
	    		}
            	Collections.sort(currComps, new Comparator<BindingSubComponents>(){
            		public int compare(BindingSubComponents o1, BindingSubComponents o2) {return o1.compareByResp(o2);}
				});
				Collections.reverse(currComps);		
							
				List<StrandedPoint> topPos = new ArrayList<StrandedPoint>();
				for (int i=0; i< Math.min(currComps.size(), config.getNumClusteringComps());i++)
					topPos.add(new StrandedPoint(currComps.get(i).getCoord(), '+'));
				
				if (topPos.size()> config.getMinComponentsForBMUpdate()){			
					shapeconfig.setStrandedPoints(topPos);			
					//Run the aligner
					ShapeAlignmentTesting profile = new ShapeAlignmentTesting(shapeconfig, gconfig, manager); 	
					try {
						profile.execute();
					} catch (FileNotFoundException | UnsupportedEncodingException e) {
						e.printStackTrace();
					}

					//Get the alignment matrix
					double[][] alignmentScores = profile.getSimilarityMatrix();
					List<StrandedRegion> regs = profile.getStrandedRegions();
					List<String> regNames= new ArrayList<String>();
					for(StrandedRegion sr : regs)
						regNames.add(sr.getLocationString());
			
					//Run AffinityPropagation
					MatrixSimilarityMeasure<Clusterable> msm = new MatrixSimilarityMeasure<Clusterable>(regNames, 
							alignmentScores, config.getPreferenceValue());
					double netsim = APCluster.cluster(msm.objects(), msm, 0.9, 200, 2000); //higher damping factor leads to slower convergence
					List<SimilarityMeasure<Clusterable>.APExemplar> exemplars = msm.getExemplars();
					List<SimilarityMeasure<Clusterable>.APAssignment> assignments = msm.getAssignments();
			
					int numExemplars = exemplars.size();
					Map<Integer, List<Integer>> clusters = 
							new HashMap<Integer, List<Integer>>();
					for(SimilarityMeasure<Clusterable>.APExemplar e : exemplars){
						clusters.put(e.index, new ArrayList<Integer>());
						clusters.get(e.index).add(e.index);
					}
			
					for(SimilarityMeasure<Clusterable>.APAssignment a : assignments){
						//System.out.println(a.index+"\t"+a.name+"\t"+a.exemplar.index);
						clusters.get(a.exemplar.index).add(a.index);
					}
			
					System.out.println("Affinity Propagation exemplars:");
					for(SimilarityMeasure<Clusterable>.APExemplar e : exemplars){
						System.out.println(e.index+"\t"+e.name+"\t"+clusters.get(e.index).size()+" members");
					}
			
			
					//OUTPUT
					//Get composites for each cluster
					int window = shapeconfig.getWindowSize();
					for(Integer c : clusters.keySet()){
						for (ExperimentCondition condition : manager.getConditions()){		
							for (ControlledExperiment rep: condition.getReplicates()){	
								double[][] composite = profile.getExemplarBasedAlignment(rep, c, clusters.get(c));
						
								try {
									File compFlie = new File(config.getOutputIntermediateDir()+File.separator+"cluster"+c+".composite.txt");
									FileWriter fout = new FileWriter(compFlie);
									fout.write("#Cluster"+c+"\n");
									for(int i=0; i<=window; i++){
										fout.write(composite[i][0]+"\t"+composite[i][1]+"\n");
									}
									fout.close();
								} catch (IOException e1) {
									e1.printStackTrace();
								}					
							}
						}
					}
			
					//Get aligned points for each cluster
					for(Integer c : clusters.keySet()){
						List<StrandedPoint> spts = profile.getExemplarBasedAlignedPoints(c, clusters.get(c));
				
						try {
							File compFlie = new File(config.getOutputIntermediateDir()+File.separator+"cluster"+c+".points");
							FileWriter fout = new FileWriter(compFlie);
							fout.write("#Cluster"+c+"\n");
							for(StrandedPoint s :spts)
								fout.write(s.getLocationString()+"\n");
							fout.close();
						} catch (IOException e1) {
							e1.printStackTrace();
						}			
						initClustPoints.add(spts);
					}	
				}
			}	
			
			// Test on clustered points
			List<List<List<StrandedPoint>>> clustPoints = new ArrayList<List<List<StrandedPoint>>>();
			for (ExperimentCondition cond : manager.getConditions())
				clustPoints.add(initClustPoints);
			
			// Make binding subtype
			for (ExperimentCondition cond : manager.getConditions()){
				List<BindingSubtype> subtypes = new ArrayList<BindingSubtype>();
				for (List<StrandedPoint> modelRefs : clustPoints.get(cond.getIndex())){
					BindingSubtype currType = new BindingSubtype(cond, modelRefs, config.MAX_BINDINGMODEL_WIDTH);
					currType.setClusteredProfile(true);
					subtypes.add(currType);
				}
				bindingManager.addPotentialBindingSubtypes(cond, subtypes);
			}
		}
	}
	
	
	/**
     * Update binding models for each replicate given the discovered binding components. 
     * @param left 
     * @param right
     * @param String filename for new distribution file
     * @return double array of log KL values
	 * @throws IOException 
     */
	public void updateBindingModelUsingReadDistributions(String distribFilename) throws IOException{
    	
    	int width = config.MAX_BINDINGMODEL_WIDTH;
    	int left = config.MAX_BINDINGMODEL_WIDTH/2;
		
		if(config.doBMUpdate()){
		
		// Choose which ones to include
		//Sum read profiles if there are enough binding components
    	for(ExperimentCondition cond : manager.getConditions()){
    		int numPeaks = 0;
    		List<List<BindingSubComponents>> currComps = new ArrayList<List<BindingSubComponents>>();
    		//Initialize
    		for (int bt=0; bt < bindingManager.getNumBindingType(cond); bt++)
    			currComps.add(new ArrayList<BindingSubComponents>());
			double currAlpha = calcAlpha(cond);
    		//Choose which components to include
    		for(Region r : activeComponents.keySet()){
    			for(BindingSubComponents bc : activeComponents.get(r).get(cond.getIndex())){
    				//1) Component must not be at the edge of the region 
    				if((bc.getPosition()-r.getStart()>config.getModelRange()/2) && (r.getEnd()-bc.getPosition()>config.getModelRange()/2)){
    					//2) Arbitrary minimum read support for BM components
    					if(bc.getSumResponsibility()>(config.getMinComponentReadFactorForBM()*currAlpha)){
    						currComps.get(bc.getMaxType()).add(bc);
    						numPeaks++;
    					}}}}			
    		
    		// Update read profile using assigned reads to a binding component per subtype group
    		List<double[][]> newModelList = new ArrayList<double[][]>();
    		List<Integer> eventCounter = new ArrayList<Integer>();
    		for (List<BindingSubComponents> compGroups : currComps){	//Iterate each subtype group
    			// Check to see if there are enough components to be used for BM updates
    			if (compGroups.size() < Math.max(config.getMinSubtypeFraction()*numPeaks, config.getMinComponentsForBMUpdate())){
    				System.out.println("Number of binding components is too few for distribution update ( "+compGroups.size()+" < "+Math.max((int) config.getMinSubtypeFraction()*numPeaks, config.getMinComponentsForBMUpdate())+" )");
    			}else{
    				eventCounter.add(compGroups.size());
    				double[] newModel_plus=new double[width];
    				double[] newModel_minus=new double[width];
    				for (int i=0;i<width;i++){
    					newModel_plus[i]=1;
    					newModel_minus[i]=1;
    				}

    				for (ControlledExperiment rep : cond.getReplicates()){
        				int x = rep.getIndex();
        				for(BindingSubComponents comp : compGroups){
        					double[] currProfile_plus = comp.getReadProfile_plus(rep.getIndex());
        		    		double[] currProfile_minus = comp.getReadProfile_minus(rep.getIndex());
        					if (comp.getMaxTypeStrand()=='+'){
        						for (int i=0; i< width; i++){
        							newModel_plus[i]+=currProfile_plus[i];
        							newModel_minus[i]+=currProfile_minus[i];
        						}
        					}else{
        						for (int i=0; i< width; i++){
        							newModel_plus[i]+=currProfile_minus[width-i-1];
        							newModel_minus[i]+=currProfile_plus[width-i-1];
        						}}}		
			    		mutate_normalize(newModel_plus,newModel_minus);
        			}
    				double[][] currModel=new double[2][];
    				currModel[0]=newModel_plus; currModel[1]=newModel_minus;
    				newModelList.add(currModel);
    			}
    		} // end of iterating each subtype group
    		
    		List<BindingSubtype> subtypes = new ArrayList<BindingSubtype>();
    		for (int index=0; index< newModelList.size(); index++){
	    			
	    		double[] newModel_plus = newModelList.get(index)[0];
	    		double[] newModel_minus = newModelList.get(index)[1];    			

	    		List<Pair<Integer,Double>> empiricalWatson = new ArrayList<Pair<Integer, Double>>(); 
	    		List<Pair<Integer,Double>> empiricalCrick = new ArrayList<Pair<Integer, Double>>();
	    		for (int i=0; i<width;i++){
	    			double data_watson = newModel_plus[i];
	    			double data_crick = newModel_minus[i];
	    			data_watson = data_watson >=0? data_watson:2.0E-300;
	    			data_crick = data_crick >=0? data_crick:2.0E-300; 
	    			Pair<Integer, Double> p_watson = new Pair<Integer, Double>(i-left, data_watson);
	    			Pair<Integer, Double> p_crick = new Pair<Integer, Double>(i-left, data_crick);
	    			empiricalWatson.add(p_watson);
	    			empiricalCrick.add(p_crick);
	    		}
		
	    		TagProbabilityDensity model = new TagProbabilityDensity(empiricalWatson,empiricalCrick);
	    		String outFile = distribFilename+"_"+cond.getName()+"_"+index+".txt";
	    		model.printDensityToFile(outFile);
	    		subtypes.add(new BindingSubtype(cond, model,eventCounter.get(index)));	    			

	    		System.err.println("Updated read distribution from " + eventCounter.get(index) +" binding events.");
	    	}  		
    		bindingManager.addPotentialBindingSubtypes(cond, subtypes);
  
    	} // end of condition loop
    	
		}else{
			System.err.println("Read distribution updates turned off.");
		}		
	}
	
	public void consolidateBindingModels(){
		
		// Merge similar binding models
		for (ExperimentCondition cond : manager.getConditions()){	
			List<BindingSubtype> currSubtypes = bindingManager.getPotentialBindingSubtypes(cond);
			if (!currSubtypes.isEmpty()){
				int numTypes = currSubtypes.size();
				System.out.println("number of binding subtype is "+numTypes);
				if (numTypes > 1){
					double[][] klScores = new double[numTypes][numTypes];
					for (int i=0; i < numTypes ; i++)
						for (int j=0; j < numTypes; j++)
							klScores[i][j]=0;
					for (int a=0; a < numTypes; a++)
						for (int b=1; b < numTypes; b++)
							klScores[a][b]=getMinKLDivergenceScore(currSubtypes.get(a).getBindingModel(0),currSubtypes.get(b).getBindingModel(0));
					// Identify similar binding model one by one
					List<Integer> model2remove = new ArrayList<Integer>();
					boolean checkSimilarity=true;
					do{
						double minScore=Double.MAX_VALUE;
						int indexA=0; int indexB=0;
						for (int i=0; i < numTypes; i++){
							for (int j=i+1; j < numTypes; j++){
								if (!model2remove.contains(i) && !model2remove.contains(j)){
									if (klScores[i][j] < minScore ){
										minScore =klScores[i][j];
										indexA=i; indexB=j;								
									}}}}
						if (minScore <config.getKLDivergenceThres()){
							BindingSubtype subA = currSubtypes.get(indexA);
							BindingSubtype subB = currSubtypes.get(indexB);
							if (subA.hasMotif() && subB.hasMotif()){
								double maxscore = motifFinder.motifAlignMaxScore(subA.getFreqMatrix(), subB.getFreqMatrix());
								System.out.println("max score is "+maxscore/config.getMinMotifLength());
								if (maxscore/config.getMinMotifLength() > config.getMotifPCCThres()){
									if (subA.getNumEvents() < subB.getNumEvents())
										model2remove.add(indexA);
									else
										model2remove.add(indexB);						
								}else{
									klScores[indexA][indexB]=Double.MAX_VALUE;
								}
							}else if (subA.hasMotif()){ // only subtype A has motif
								model2remove.add(indexB);
							}else if (subB.hasMotif()){ // only subtype B has motif
								model2remove.add(indexA);
							}else{						// subtype A and subtype B do not have motif
								if (subA.getNumEvents() < subB.getNumEvents())
									model2remove.add(indexA);
								else
									model2remove.add(indexB);	
							}
						}else{
							checkSimilarity=false;
						}							
					}while (checkSimilarity); // End of while loop
				
					System.out.println("models to be removed "+model2remove.toString());
					List<BindingSubtype> subs2remove = new ArrayList<BindingSubtype>();
					for (Integer i : model2remove){
						subs2remove.add(currSubtypes.get(i));
					}
					currSubtypes.removeAll(subs2remove);
					System.out.println("Models are consolidated from "+numTypes+" to "+currSubtypes.size());
				}
				bindingManager.setBindingSubtypes(cond, currSubtypes);	
				bindingManager.clearPotentialBindingSubtypes(cond);
			}
		}
	}
	
	// this method will mutate the input array
	public static void mutate_normalize(double[] w_dist, double[] c_dist){
		double total=0;
		for (int i=0;i<w_dist.length;i++){
			if (w_dist[i]<=0)
				w_dist[i]=1e-20;
			total += w_dist[i];
		}
		for (int i=0;i<c_dist.length;i++){
			if (c_dist[i]<=0)
				c_dist[i]=1e-20;
			total += c_dist[i];
		}
		for(int i=0;i<w_dist.length;i++){
			w_dist[i] /=total;
		}
		for(int i=0;i<c_dist.length;i++){
			c_dist[i] /=total;
		}
	}

    /**
     * Update condition and replicate backgrounds for alphas 
     */
    public void updateAlphas(){
    	for(ExperimentCondition cond : manager.getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
			conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), cond.getTotalSignalCount()*(1-cond.getTotalSignalVsNoiseFrac()), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), config.getModelRange(), '.', 1, true));
			double condAlf = config.getFixedAlpha()>0 ? config.getFixedAlpha() : (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
			double alf=condAlf;
			
			//In lenient mode, we have to allow alpha to be small enough to catch events significant in individual replicates
			double repAlf = Double.MAX_VALUE;	
			if(config.getFixedAlpha()<=0 && (config.isLenientMode() || config.isLenientPlusMode())){
				for(ControlledExperiment rep : cond.getReplicates()){
					replicateBackgrounds.put(rep, new BackgroundCollection());
					replicateBackgrounds.get(rep).addBackgroundModel(new PoissonBackgroundModel(-1, config.getSigLogConf(), rep.getSignal().getHitCount()*(1-rep.getSignalVsNoiseFraction()), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), config.getModelRange(), '.', 1, true));
					repAlf = Math.min(repAlf, (double)replicateBackgrounds.get(rep).getMaxThreshold('.'));
				}
				alf = Math.min(condAlf,  repAlf);
			}
			
			System.err.println("Alpha "+cond.getName()+"\tRange="+config.getModelRange()+"\t"+alf);
			bindingManager.setAlpha(cond,alf);
		}
    }
    
    /**
     * Calculate alpha value
     */
    private double calcAlpha(ExperimentCondition cond) {
    	double condAlf = config.getFixedAlpha()>0 ? config.getFixedAlpha() : (double)conditionBackgrounds.get(cond).getMaxThreshold('.');
		double alf=condAlf;	    			
		//In lenient mode, we have to allow alpha to be small enough to catch events significant in individual replicates
		double repAlf = Double.MAX_VALUE;	
		if(config.getFixedAlpha()<=0 && (config.isLenientMode() || config.isLenientPlusMode())){
			for(ControlledExperiment rep : cond.getReplicates()){
				repAlf = Math.min(repAlf, (double)replicateBackgrounds.get(rep).getMaxThreshold('.'));
			}
			alf = Math.min(condAlf,  repAlf);
		}
		return alf;
    }
    
    /**
     * Run motif-finding, given the current BindingSubComponentss
     * @throws Exception 
     */
    public void updateBindingModelUsingMotifs() throws Exception{
    	if(config.getFindingMotifs()){
    		motifFinder.updateSubtypesUsingMotifs(activeComponents, trainingRound);
    		
    		//Print progress
    		for(ExperimentCondition cond : manager.getConditions()){
    			boolean noMotif=true;
    			for (BindingSubtype sub : bindingManager.getBindingSubtype(cond)){
    				if (sub.hasMotif()){
    					System.err.println(cond.getName()+"\t"+WeightMatrix.getConsensus(sub.getFreqMatrix())+"\toffset:"+sub.getMotifOffset());
    					noMotif=false; break;
    				}
    			}
    			if (noMotif) {System.err.println(cond.getName()+"\tNOMOTIF");}
    		}
    	}
    }

    /**
     * Initialize the global noise parameters. Inferred either from:
     *  - non-potential region read counts, or
     *  - noise proportion of reads from SES
     */
    protected void initializeGlobalNoise(){
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getIndexedCondition(e);
    		
    		// Part that deals with read counts from non-potential regions... calculate values anyway whether using them or not
    		double potRegLengthTotal = potRegFilter.getPotRegionLengthTotal();
    		double nonPotRegLengthTotal = config.getGenome().getGenomeLength() - potRegLengthTotal;
    		//Combine control channel counts (avoiding duplication)
    		double potRegCountsSigChannel=potRegFilter.getPotRegCountsSigChannel(cond);
    		double nonPotRegCountsSigChannel=potRegFilter.getNonPotRegCountsSigChannel(cond); 
    		double potRegCountsCtrlChannel=potRegFilter.getPotRegCountsCtrlChannel(cond);
    		double nonPotRegCountsCtrlChannel=potRegFilter.getNonPotRegCountsCtrlChannel(cond); 
    		List<Sample> ctrls = cond.getControlSamples();
    		
    		//relativeCtrlNoise just tells us if there is a systemic over/under representation of reads in potential regions (in the control)
    		//NOTE: not used for anything right now. 
    		relativeCtrlNoise[e] = (potRegCountsCtrlChannel==0 && nonPotRegCountsCtrlChannel==0) ? 
    				1 : (potRegCountsCtrlChannel/potRegLengthTotal)/(nonPotRegCountsCtrlChannel/nonPotRegLengthTotal);

    		
    		noisePerBase[e] = nonPotRegCountsSigChannel/nonPotRegLengthTotal;  //Signal channel noise per base
    		System.err.println("Global noise per base initialization for "+cond.getName()+" = "+String.format("%.4f", noisePerBase[e]));
    	}
    }
    
    /**
     * Update the global noise parameters, using both non-potential region counts and assigned noise responsibilities
     */
    public void updateGlobalNoise(){
    	for(int e=0; e<manager.getNumConditions(); e++){
    		ExperimentCondition cond = manager.getIndexedCondition(e);
    		
    		//Don't need to examine noise reads in the update
    		double noiseReads=potRegFilter.getNonPotRegCountsSigChannel(cond); 
    		
    		for(Region r : noiseResp.keySet())
    			noiseReads+=noiseResp.get(r)[e];
    		noisePerBase[e] = noiseReads/config.getGenome().getGenomeLength();  //Signal channel noise per base
    	}
    }
    

       
    /**
     * Print reference points to a file.
     * TESTING ONLY 
     */
    public void printPointsToFile(String fileName, List<StrandedPoint> points){
    	try {
			FileWriter fout = new FileWriter(fileName);
			for (StrandedPoint p : points)
				fout.write(p.toString()+"\n");
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }  
	
    /**
     * Print all components active at the current time to a file.
     * TESTING ONLY 
     */
    public void printActiveComponentsToFile(){
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_t"+trainingRound+".components";
			FileWriter fout = new FileWriter(filename);
			for(Region rr : activeComponents.keySet()){
	    		List<List<BindingSubComponents>> comps = activeComponents.get(rr);
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			for(BindingSubComponents comp : comps.get(cond.getIndex())){
	    				fout.write(rr.getLocationString()+"\t"+cond.getName()+"\t"+comp.toString()+"\n");			
	    			}
	    		}
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
    /**
     * Print all components active at the current time.
     * TESTING ONLY 
     */
    public void printActiveComponents(){
    	for(Region rr : activeComponents.keySet()){
    		List<List<BindingSubComponents>> comps = activeComponents.get(rr);
    		for(ExperimentCondition cond : manager.getConditions()){
    			for(BindingSubComponents comp : comps.get(cond.getIndex())){
    				System.err.println(rr.getLocationString()+"\t"+cond.getName()+"\t"+comp.toString());			
    			}
    		}
    	}
    }
    

    /**
     * Print all reads relative to reference region (ie motif).
     * TESTING ONLY 
     */
    public void printReadsAroundRegionsToFile(){
    	int length = 100; // output histogram length
    	List<StrandedRegion> refRegion = config.getMotifRegions();
    	List<StrandedRegion> expandedRegion = new ArrayList<StrandedRegion>();
    	for (StrandedRegion r : refRegion)
    		expandedRegion.add(r.expand(length/2, length/2));   	
    	int histoLength=expandedRegion.get(0).getWidth();       	
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_readsAroundRefRegion.txt";
    		FileWriter fout = new FileWriter(filename);
    		fout.write("#reads within region distance histograms\n\n");
    		// data is indexed by rep, hashed by region
    		List<Map<Region,List<StrandedBaseCount>>> repRegionData = new ArrayList<Map<Region,List<StrandedBaseCount>>>();
    		for(ExperimentCondition cond : manager.getConditions()){
    			for(ControlledExperiment rep : cond.getReplicates()){
    				repRegionData.add(new HashMap<Region,List<StrandedBaseCount>>());
    				for (StrandedRegion r : expandedRegion)
    					repRegionData.get(rep.getIndex()).put(r, new ArrayList<StrandedBaseCount>());
    			}
    		}
    		// add base counts
    		for(ExperimentCondition cond : manager.getConditions()){
    			for(ControlledExperiment rep : cond.getReplicates()){
    				for (StrandedRegion r : expandedRegion)
    					repRegionData.get(rep.getIndex()).get(r).addAll(rep.getSignal().getBases(r));
    			}
    		}
    		// create histogram separately for positive and negative counts
    		for(ExperimentCondition cond : manager.getConditions()){
    			fout.write("#Condition:"+cond.getName()+"\n");
    			for(ControlledExperiment rep : cond.getReplicates()){
    				//add positive counts in histogram
    				RealValuedHistogram poshisto = new RealValuedHistogram(0, histoLength, histoLength);
    				RealValuedHistogram neghisto = new RealValuedHistogram(0, histoLength, histoLength);
    				for (StrandedRegion r : expandedRegion){
    					for(StrandedBaseCount b :repRegionData.get(rep.getIndex()).get(r)){
    						int index=0;
    						if (r.getStrand()=='+'){ // if region strand is positive
    							index=b.getCoordinate()-r.getStart();
    							if (b.getStrand()=='+')
    								poshisto.addValue(index, b.getCount());
    							else
    								neghisto.addValue(index, b.getCount());
    						}else{ // if region strand is negative
    							index=r.getEnd()-b.getCoordinate();
    							if (b.getStrand()=='+')
    								neghisto.addValue(index, b.getCount());
    							else
    								poshisto.addValue(index, b.getCount());
    						}
    					}
    				}
    				// write to file
    				fout.write("#Positive hits"+"\n");
    				fout.write(poshisto.contentsToString()+"\n");
    				fout.write("#Negative hits"+"\n");
    				fout.write(neghisto.contentsToString()+"\n");   				
    			}
    		}
    		fout.close();
    	} catch (IOException e) {
			e.printStackTrace();
		}
    }
    
	/**
	 * BindingMixtureThread: run binding mixtures over a subset of regions
	 * uniformComponents: if true, components are initialized at uniform prob and spacing, 
	 * 					otherwise, components are initialized from activeComponents
	 * bindingEvents: if true, run EM over binding events, otherwise run EM over cross-linking points
	 * @author Naomi Yamada
	 * @version	%I%, %G%
	 */
	class BindingMixtureThread implements Runnable{
		private Collection<Region> regions;
		private int numBindingComponents=1;	//Assuming that the total number of components (active+inactive) is the same in every condition makes coding easier in the BindingEM class.  
		private boolean runEM = true;
		private boolean uniformBindingComponents=false;
		private boolean runMultiGPSML = false;
		
		public BindingMixtureThread(Collection<Region> regs, boolean EM, boolean uniformBindingComponents, boolean multiGPSML){
			regions = regs;	
			this.uniformBindingComponents = uniformBindingComponents;
			runEM=EM;
			runMultiGPSML = multiGPSML;
		}
		
		/**
		 * Run the binding mixture over each test region
		 */
		public void run() {
			//For each region in the test set
        	for (Region rr : regions) {
        		try{
        			//Initialize array of binding component lists, indexed by condition

        			List<List<BindingSubComponents>> currComps = new ArrayList<List<BindingSubComponents>>();
        			for(int e=0; e<manager.getNumConditions(); e++)
        				currComps.add(new ArrayList<BindingSubComponents>());
                    
                    if(runEM ){
		        		//Run EM
	                    Double[] noiseRSums = new Double[manager.getNumConditions()];
	                    for(int e=0; e<manager.getNumConditions(); e++){ noiseRSums[e]=0.0;}
		        		Pair<List<NoiseComponent>, List<List<BindingSubComponents>>> wComps = analyzeWindowEM(rr);
		        		for(int e=0; e<manager.getNumConditions(); e++){
		        			noiseRSums[e] += wComps.car().get(e).getSumResponsibility();
		        			currComps.get(e).addAll(wComps.cdr().get(e));
	                    }
		        		
		        		//Only non-zero components are returned by analyzeWindow, so add them to the recorded active components
		        		synchronized(activeComponents){	activeComponents.put(rr, currComps);}
		        		
		        		//Add the sum of noise responsibilities to this region
		        		synchronized(noiseResp){ noiseResp.put(rr, noiseRSums);}
		        		
                    }else{
                    	//Run ML assignment
                    	List<BindingEvent> windowBindingEvents = new ArrayList<BindingEvent>(); 
                    	windowBindingEvents.addAll(  analyzeWindowML(rr) ); 

                    	synchronized(bindingEvents){bindingEvents.addAll(windowBindingEvents);}
                    }
        		} catch(Exception e){
                    System.err.println("ERROR: Exception when analyzing region "+rr.toString());
                    e.printStackTrace(System.err);
                    System.exit(-1);
                }
            }
		}
		
		/**
		 * Train a BindingComponent EM over a given window
		 * 
		 * In the old implementation, there was a choice here to run EM or scan for the ML soln
         * (if we could be sure there was only one binding event). However, we can no longer scan because
         * the solution can depend on multiple conditions. So we always do EM.
         * 
         * In addition, since we are now using components that can change their position along the genome, 
         * there is no need for the component resolution updates. 
         * 
         * Components are now semi-independently estimated for each condition, so they have become lists-of-lists.
         * 
         * We also now initialize noise components, which are position-less and have a fixed (estimated) emission probability.
         *  
		 * @param w
		 * @return Pair of component lists (noise components and binding components) indexed by condition
		 * @throws FileNotFoundException 
		 */
		private Pair<List<NoiseComponent>, List<List<BindingSubComponents>>> analyzeWindowEM(Region w) throws FileNotFoundException{
			BindingEM EM = new BindingEM(config, manager, bindingManager, conditionBackgrounds, replicateBackgrounds, potRegFilter.getPotentialRegions().size());
			MultiGPSBindingEM multiGPSEM = new MultiGPSBindingEM(config, manager, bindingManager, conditionBackgrounds, replicateBackgrounds,  potRegFilter.getPotentialRegions().size());
			
			List<List<BindingSubComponents>> bindingComponents=null;
			List<NoiseComponent> noiseComponents=null;
			List<List<BindingSubComponents>> nonZeroComponents = new ArrayList<List<BindingSubComponents>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				nonZeroComponents.add(new ArrayList<BindingSubComponents>());
			
			//Check if this is a region to plot (and if so, calculate any reqd offsets)
			Region plotSubReg=null;
			for(Region p : regionsToPlot)
				if(w.overlaps(p))
					plotSubReg = p;
			
			//Load signal data
			List<List<StrandedBaseCount>> signals = loadSignalData(w);
            if (signals==null)
                return new Pair<List<NoiseComponent>, List<List<BindingSubComponents>>>(noiseComponents, nonZeroComponents);
            //Load control data
            List<List<StrandedBaseCount>> controls = loadControlData(w);
            
            //Initialize noise components
            noiseComponents = initializeNoiseComponents(w, signals, controls);

            //Initialize binding components
            if(uniformBindingComponents)
            	bindingComponents = initializeBindingComponentsUniformly(w, noiseComponents);
            else
            	bindingComponents = initializeBindingComponentsFromAllConditionActive(w, noiseComponents, config.getAddFlankingComponents());
            
            //Motif prior
            String seq = config.getFindingMotifs() ? motifFinder.getSeq(w):null;
            double[][][] forMotifPrior = config.getFindingMotifs() ? motifFinder.scanStrandedRegionWithMotifs(w, seq, true) : null;
            double[][][] revMotifPrior = config.getFindingMotifs() ? motifFinder.scanStrandedRegionWithMotifs(w, seq, false) : null;
                        
            //EM learning: resulting binding components list will only contain non-zero components   
            if (uniformBindingComponents && mixconfig.getModelFilename()== null)
            	nonZeroComponents = multiGPSEM.train(signals, w, noiseComponents, bindingComponents, numBindingComponents, trainingRound, plotSubReg);
            else
            	nonZeroComponents = EM.train(signals, w, noiseComponents, bindingComponents, numBindingComponents, forMotifPrior, revMotifPrior, trainingRound, plotSubReg);
           
            return new Pair<List<NoiseComponent>, List<List<BindingSubComponents>>>(noiseComponents, nonZeroComponents);
        }//end of analyzeWindowEM method
		

		/**
		 * Assign BindingComponents over a given window with ML solution
		 *  
		 * @param w
		 * @return Pair of component lists (noise components and binding components) indexed by condition
		 */
		private List<BindingEvent> analyzeWindowML(Region w){
			BindingMLAssignment ML = new BindingMLAssignment(econfig, evconfig, config, manager,bindingManager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
			MultiGPSMLAssignment GPSML = new MultiGPSMLAssignment(econfig, evconfig, config, manager,bindingManager, conditionBackgrounds, potRegFilter.getPotentialRegions().size());
			List<BindingSubComponents> bindingComponents=null;
			List<NoiseComponent> noiseComponents=null;
			List<BindingEvent> currEvents = new ArrayList<BindingEvent>(); 
			
			//Load signal data
			List<List<StrandedBaseCount>> signals = loadSignalData(w);
            if (signals==null)
                return currEvents;
            //Load control data
            List<List<StrandedBaseCount>> controls = loadControlData(w);
            
            //Initialize noise components
            noiseComponents = initializeNoiseComponents(w, signals, controls);
            
            //Configurations seen in another condition
            ArrayList<ComponentConfiguration> seenConfigs = new ArrayList<ComponentConfiguration>();
    		//Assign reads to components
            for(ExperimentCondition cond : manager.getConditions()){
            	//Initialize binding components: condition-specific  
            	if (runMultiGPSML)
            		bindingComponents = initializeBindingComponentsFromOneConditionActive(w, noiseComponents.get(cond.getIndex()), cond.getIndex());
            	else
            		bindingComponents = initializeBindingSubComponentsFromOneConditionActive(w, noiseComponents.get(cond.getIndex()), cond.getIndex());
            	
            	int numComp = bindingComponents.size();
            	
            	//Construct configuration
    			ComponentConfiguration currCC = new ComponentConfiguration(bindingComponents, cond.getIndex());
    		
    			//Have we already seen this configuration?
    			boolean ccFound=false;
    			for(ComponentConfiguration cc : seenConfigs){
    				if(currCC.isSameAs(cc)){//If so, add another valid condition to the observed binding events
    					int parent = cc.getParentCondition();
    					for(BindingEvent be : currEvents)
    						if(be.isFoundInCondition(parent))
    							be.setIsFoundInCondition(cond.getIndex(),true);
    					ccFound=true;
    					break;
    				}
    			}
    			
    			if(!ccFound){
    				//Add configuration to seen
    				seenConfigs.add(currCC);
    				
    				//ML assignment
    				List<BindingEvent> condEvents = null;
    				if (runMultiGPSML){
    					condEvents = GPSML.assign(signals, controls, w, noiseComponents, bindingComponents, numComp, cond);  					
    				}else{
    					//Make ML assignment condition specific or hack to get some read assignment
        				condEvents = ML.assign(signals, controls, w, noiseComponents, bindingComponents, numComp, cond);
    				}
    				
    				for(BindingEvent be : condEvents)
    					be.setIsFoundInCondition(cond.getIndex(),true);
    				currEvents.addAll(condEvents);
    			}
            }
            
            //If we haven't used shared component ML, we need to edit and consolidate binding events
            // 1) Consolidate, because otherwise you can have duplicate binding events
            // 2) Edit - set counts to zero at conditions where the event is not active
            currEvents = consolidateBindingEvents(currEvents);

            
            //Add in sequences and final motif scores here
            if(evconfig.isAddingSequences()){
	            String seq = config.getFindingMotifs() ? motifFinder.getSeq(w):null;
	            Pair<Double[][][], String[][][]> motifForScores = config.getFindingMotifs() ? motifFinder.scanRegionWithMotifsGetSeqs(w, seq, true) : null;
	            Pair<Double[][][], String[][][]> motifRevScores = config.getFindingMotifs() ? motifFinder.scanRegionWithMotifsGetSeqs(w, seq, false) : null;
	            if(seq!=null && motifForScores!=null && motifRevScores!=null){
		            for(ExperimentCondition cond : manager.getConditions()){
		            	int numBindingType = bindingManager.getNumBindingType(cond);
		            	for(BindingEvent b: currEvents){
		            		double[][] scores = new double[numBindingType][2];
			            	String[][] seqs = new String[numBindingType][2];
			            	if (b.isFoundInCondition(cond)){
		            			for (int bt=0; bt < numBindingType; bt++){
		            				BindingSubtype subtype = bindingManager.getBindingSubtype(cond).get(bt);
		            				scores[bt][0]=0; seqs[bt][0] = ""; scores[bt][1]=0; seqs[bt][1]="";
		            				if (subtype.hasMotif() && b.getTypePoints(cond)!=null){
		            					if (b.getTypePoints(cond)[bt][0].getLocation()-w.getStart()>0){
		            						scores[bt][0] = motifForScores.car()[cond.getIndex()][bt][b.getTypePoints(cond)[bt][0].getLocation()-w.getStart()];
		            						seqs[bt][0] = motifForScores.cdr()[cond.getIndex()][bt][b.getTypePoints(cond)[bt][0].getLocation()-w.getStart()];
		            					}
		            					if (b.getTypePoints(cond)[bt][1].getLocation()-w.getStart()>0){
		            						scores[bt][1] = motifRevScores.car()[cond.getIndex()][bt][b.getTypePoints(cond)[bt][1].getLocation()-w.getStart()];
		            						seqs[bt][1] = motifRevScores.cdr()[cond.getIndex()][bt][b.getTypePoints(cond)[bt][1].getLocation()-w.getStart()];
		            					}
		            				}
		            			}
		            		}
		            		b.setMotifScore(cond,scores);
		            		b.setSequence(cond, seqs);	          			
		            	}
		            }
	            }
            }

            return currEvents;
        }//end of analyzeWindowML method

		/**
		 * Set the conditions in which a given binding event is still active after EM training
		 * @param b
		 * @param currReg
		 */
		private void setFoundInConditions(BindingEvent b, Region currReg){
			for(int e=0; e<manager.getNumConditions(); e++)
        		for(BindingSubComponents comp : activeComponents.get(currReg).get(e)){
        			if(comp.getPosition() == b.getPoint().getLocation()){
        				b.setIsFoundInCondition(e, true);
        			}
        		}
		}
		
		/**
		 * Load all signal read hits in a region by condition. 
		 * 
		 * @param w
		 * @return List of List of StrandedBaseCounts, indexed by replicate index
		 */
		private List<List<StrandedBaseCount>> loadSignalData(Region w){
			List<List<StrandedBaseCount>> data = new ArrayList<List<StrandedBaseCount>>();
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedBaseCount>());
			}
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					data.get(rep.getIndex()).addAll(rep.getSignal().getBases(w));
				}
			}
			return data;
		}
		
		/**
		 * Load all control read hits in a region by condition. 
		 * 
		 * @param w
		 * @return List of List of StrandedBaseCounts, indexed by replicate index
		 */
		private List<List<StrandedBaseCount>> loadControlData(Region w){
			List<List<StrandedBaseCount>> data = new ArrayList<List<StrandedBaseCount>>();
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates())
					data.add(new ArrayList<StrandedBaseCount>());
			}
			for(ExperimentCondition cond : manager.getConditions()){
				for(ControlledExperiment rep : cond.getReplicates()){
					if(rep.hasControl())
						data.get(rep.getIndex()).addAll(rep.getControl().getBases(w));
				}
			}
			return data;
		}
		
		/**
         * Initializes the components uniformly: i.e. space them evenly along the region.
         *
         * @param currReg
         */
        private List<List<BindingSubComponents>> initializeBindingComponentsUniformly(Region currReg, List<NoiseComponent> noise){
        	List<List<BindingSubComponents>> components = new ArrayList<List<BindingSubComponents>>();
			for(int e=0; e<manager.getNumConditions(); e++)
				components.add(new ArrayList<BindingSubComponents>());
			
            //Place components along region
			int componentSpacing =  config.getInitialCompSpacing();
			if(componentSpacing >= currReg.getWidth())
				System.err.println("Error:  region width less than component spacing in "+currReg.getLocationString());

			numBindingComponents = currReg.getWidth()/componentSpacing;

            //Set up the components
			for(int e=0; e<manager.getNumConditions(); e++){
	            double numC=0;
	            for(int i=0; i<numBindingComponents; i++){
	            	Point pos = new Point(config.getGenome(), currReg.getChrom(), currReg.getStart()+(i*componentSpacing));
	            	BindingSubComponents currComp = new BindingSubComponents(pos, manager.getReplicates().size());
	                currComp.setIndex(i);
	                numC++;
	                components.get(e).add(currComp);
	            }
	            //Initialize normalized mixing probabilities (subtracting the noise emission probability)
	            double emission = (1-noise.get(e).getPi())/numC;
	            for(BindingSubComponents b : components.get(e)){
	                b.uniformInit(emission);	                
	            }
			}
            return components; 
        }//end of initializeComponents method

        /**
         * Initializes the components from all active components in all conditions: 
         * 		Uses all active component locations from the last round of training in each condition,
         * 		i.e. not just the conditions in which those active components were active. Also adds in 
         * 		extra components flanking the active locations in case the binding distribution update
         * 		has made more joint events separable. If no components exist, a rescue component is added. 
         *
         * @param currReg
         */
        private List<List<BindingSubComponents>> initializeBindingComponentsFromAllConditionActive(Region currReg, List<NoiseComponent> noise, boolean addFlanking){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(int e=0; e<manager.getNumConditions(); e++){
        		for(BindingSubComponents comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd())
        				componentPositions.add(comp.getPosition());

        			if(addFlanking){
        				if(!componentPositions.contains(comp.getPosition()-config.getAddFlankingComponentSpacing())
        						 && comp.getPosition()-config.getAddFlankingComponentSpacing()>=currReg.getStart())
        					componentPositions.add(comp.getPosition()-config.getAddFlankingComponentSpacing());
        				
        				if(!componentPositions.contains(comp.getPosition()+config.getAddFlankingComponentSpacing())
        						&& comp.getPosition()+config.getAddFlankingComponentSpacing()<currReg.getEnd())
        					componentPositions.add(comp.getPosition()+config.getAddFlankingComponentSpacing());						       				
        			}
        		}
        	}

        	numBindingComponents = componentPositions.size();
        	
        	//If no components exist in region, add one to the center to allow rescues
        	if(numBindingComponents==0 && addFlanking){
        		componentPositions.add(currReg.getMidpoint().getLocation());
        		numBindingComponents++;     		
        	}
        	  
        	//Make new components with these locations
        	List<List<BindingSubComponents>> components = new ArrayList<List<BindingSubComponents>>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		components.add(new ArrayList<BindingSubComponents>());
        	
        	//Set up the components
        	for(int e=0; e<manager.getNumConditions(); e++){
        		double numC=(double)numBindingComponents;
        		double emission = (1-noise.get(e).getPi())/numC;
        		        		
	    		for(int index=0; index<componentPositions.size(); index++){
	    			Point pos = new Point(config.getGenome(), currReg.getChrom(),componentPositions.get(index));  
	    			BindingSubComponents currComp = new BindingSubComponents(pos, manager.getReplicates().size());
	    			currComp.setIndex(index);
	    			
    				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
        			currComp.uniformInit(emission);
    				components.get(e).add(currComp);
    			}
    		}
        	return components; 
        }//end of initializeComponents method
        
        /**
         * Initializes components from active components in a single condition: 
         * 		Uses active component locations from the last round of training in one condition,
         * 		No flanking components or resuce components added here, since resulting components will only be used
         * 		in ML assignment.  
         *
         * @param currReg
         */
        private List<BindingSubComponents> initializeBindingComponentsFromOneConditionActive(Region currReg, NoiseComponent noise, int conditionIndex){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	for(BindingSubComponents comp : activeComponents.get(currReg).get(conditionIndex)){
        		if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd()){
        			componentPositions.add(comp.getPosition());
        		}
        	}

        	numBindingComponents = componentPositions.size();

        	//Make new components with these locations
        	List<BindingSubComponents> components = new ArrayList<BindingSubComponents>();
        	
        	//Set up the components
        	double numC=(double)numBindingComponents; 
    		double emission = (1-noise.getPi())/numC;
    		for(int index=0; index<componentPositions.size(); index++){
    			Point pos = new Point(config.getGenome(), currReg.getChrom(),componentPositions.get(index));  
    			BindingSubComponents currComp = new BindingSubComponents(pos, manager.getReplicates().size());
    			currComp.setIndex(index);  	
				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
    			currComp.uniformInit(emission);
    			components.add(currComp);
			}
    		
        	return components; 
        }//end of initializeComponents method

        
        /**
         * Initializes the sub components from all active components in all conditions: 
         * 		Uses all active component locations from the last round of training in each condition,
         * 		i.e. not just the conditions in which those active components were active. 
         *
         * @param currReg
         */
        private List<List<BindingSubComponents>> initializeBindingSubComponentsFromAllConditionActive(Region currReg, List<NoiseComponent> noise){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	List<int[][]> subComponentPositions = new ArrayList<int[][]>();
        	List<double[][]> subComponentTaus = new ArrayList<double[][]>();
        	for(int e=0; e<manager.getNumConditions(); e++){
        		for(BindingSubComponents comp : activeComponents.get(currReg).get(e)){
        			if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd()){
        				componentPositions.add(comp.getPosition());
        				subComponentPositions.add(comp.getPositions());
        				subComponentTaus.add(comp.getTau());
        	}}}		

        	numBindingComponents = componentPositions.size();
  
        	//Make new components with these locations
        	List<List<BindingSubComponents>> components = new ArrayList<List<BindingSubComponents>>();
        	for(int e=0; e<manager.getNumConditions(); e++)
        		components.add(new ArrayList<BindingSubComponents>());
        	
        	//Set up the components
        	for(int e=0; e<manager.getNumConditions(); e++){
        		double numC=(double)numBindingComponents;
        		double emission = (1-noise.get(e).getPi())/numC;
	    		for(int index=0; index<componentPositions.size(); index++){
	    			Point pos = new Point(config.getGenome(), currReg.getChrom(),componentPositions.get(index));  
	    			BindingSubComponents currComp = new BindingSubComponents(pos, manager.getReplicates().size());
	    			currComp.setIndex(index);
	    			currComp.setPositions(subComponentPositions.get(index));
	    			currComp.setTau(subComponentTaus.get(index));
	    			
    				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
        			currComp.uniformInit(emission);
    				components.get(e).add(currComp);
    			}
    		}
        	return components; 
        }//end of initializeComponents method

        /**
         * Initializes components from active components in a single condition: 
         * 		Uses active component locations from the last round of training in one condition,
         * 		No flanking components or resuce components added here, since resulting components will only be used
         * 		in ML assignment.  
         *
         * @param currReg
         */
        private List<BindingSubComponents> initializeBindingSubComponentsFromOneConditionActive(Region currReg, NoiseComponent noise, int conditionIndex){
        	//Initialize component positions with active locations
        	List<Integer> componentPositions = new ArrayList<Integer>();
        	List<int[][]> subComponentPositions = new ArrayList<int[][]>();
        	List<double[][]> subComponentTaus = new ArrayList<double[][]>();
        	for(BindingSubComponents comp : activeComponents.get(currReg).get(conditionIndex)){
        		if(!componentPositions.contains(comp.getPosition()) && comp.getPosition()>=currReg.getStart() && comp.getPosition()<currReg.getEnd()){
        			componentPositions.add(comp.getPosition());
        			subComponentPositions.add(comp.getPositions());
        			subComponentTaus.add(comp.getTau());
        		}
        	}

        	numBindingComponents = componentPositions.size();

        	//Make new components with these locations
        	List<BindingSubComponents> components = new ArrayList<BindingSubComponents>();
        	
        	//Set up the components
        	double numC=(double)numBindingComponents; 
    		double emission = (1-noise.getPi())/numC;
    		for(int index=0; index<componentPositions.size(); index++){
    			Point pos = new Point(config.getGenome(), currReg.getChrom(),componentPositions.get(index));  
    			BindingSubComponents currComp = new BindingSubComponents(pos, manager.getReplicates().size());
    			currComp.setIndex(index);  	
    			currComp.setPositions(subComponentPositions.get(index));
    			currComp.setTau(subComponentTaus.get(index));
				//Initialize normalized mixing probabilities (subtracting the noise emission probability)
    			currComp.uniformInit(emission);
    			components.add(currComp);
			}
    		
        	return components; 
        }//end of initializeComponents method
        
        
        /**
         * Initializes the noise components.
         * Noise distributions are set per replicate, so control channel reads are not combined into a single list. 
         *
         * @param currReg
         */
        private List<NoiseComponent> initializeNoiseComponents(Region currReg, List<List<StrandedBaseCount>> sigHits, List<List<StrandedBaseCount>> ctrlHits){
        	List<NoiseComponent> noise = new ArrayList<NoiseComponent>();
        	int numReps = manager.getReplicates().size();
        	double [] localSigRepCounts=new double [numReps];
        	double [] localCtrlRepCounts=new double [numReps];
        	
        	//Calculate expected noise distributions
        	double [][] distribs=new double[numReps][];
        	for(ExperimentCondition cond : manager.getConditions())
        		for(ControlledExperiment rep : cond.getReplicates()){
	    			if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
	            		distribs[rep.getIndex()] = smoothNoiseDistribs(currReg, ctrlHits.get(rep.getIndex()));
	            		localCtrlRepCounts[rep.getIndex()]=0;
	            		for(StrandedBaseCount b : ctrlHits.get(rep.getIndex()))
	            			localCtrlRepCounts[rep.getIndex()]+=b.getCount();
	    			}else
	    				distribs[rep.getIndex()] = null;
	    		}
        	
        	//Initialize the noise component
        	for(int e=0; e<manager.getNumConditions(); e++){
        		ExperimentCondition cond = manager.getIndexedCondition(e);
        		double emission = 0;
        		
        		//Sum signal reads & set local region experiment counts
        		double sigCounts=0;
        		for(ControlledExperiment rep : cond.getReplicates()){
        			localSigRepCounts[rep.getIndex()]=0;
        			for(StrandedBaseCount b : sigHits.get(rep.getIndex())){
        				sigCounts+=b.getCount(); localSigRepCounts[rep.getIndex()]+=b.getCount();
        			}
        		}
        		
        		//Calculate a local noise factor to check for expected over-representation of noise reads, as specified in the control.
        		//We assume that only noise generates the control channel. Therefore, the first two terms in the localNoiseFactor calculation
        		//test for over-representation in the observed control read counts in the local window. 
        		//The last term in the calculation is a local replicate weight, to account for the fact that some replicate signal channels are 
        		//producing more of the signal reads (and of course, the noise we are actually accounting for is in the signal channel). 
        		double localNoiseFactor = 0;
        		for(ControlledExperiment rep : cond.getReplicates()){
        			if(rep.hasControl() && ctrlHits.get(rep.getIndex()).size()>0){
        				localNoiseFactor+=(localCtrlRepCounts[rep.getIndex()]/(double)currReg.getWidth()) /
        								  (rep.getControl().getHitCount()/(double)currReg.getWidth())     *
        								  (localSigRepCounts[rep.getIndex()]/sigCounts); //over-rep x weight
        			}else{
        				localNoiseFactor+=(localSigRepCounts[rep.getIndex()]/sigCounts); //1 x weight
        			}
        		}
        		
        		//Calculate expected noise emission = number of expected noise reads over the total signal reads in this region
        		// If local noise factor is above 1, account for it. Otherwise, meh.
        		emission = (noisePerBase[e] * (double)currReg.getWidth()) / sigCounts;
        		if(localNoiseFactor>1)
        			emission*=localNoiseFactor;
        		if(emission>config.NOISE_EMISSION_MAX)
        			emission = config.NOISE_EMISSION_MAX;
        		if(emission<config.NOISE_EMISSION_MIN)
        			emission = config.NOISE_EMISSION_MIN;
        		
        		//Add the noise component
        		NoiseComponent n = new NoiseComponent(emission, distribs, currReg, numReps);
        		noise.add(n);
        	}
        	return noise;
        }//end of initializeNoiseComponents method
      
        /**
         * Smooth the distribution of the control reads over the window to make a probability distribution of noise over the region
         * @param currReg
         * @param ctrlHits
         * @return double array the same width as the region, containing probabilities normalized to sum to 1
         */
        private double[] smoothNoiseDistribs(Region currReg, List<StrandedBaseCount> ctrlHits){
        	double [] distrib = new double[currReg.getWidth()];
        	double [] counts = new double[currReg.getWidth()];
        	//Pseudocounts for distrib
        	for(int d=0; d<currReg.getWidth(); d++)
        		counts[d]=1;
        	//Add in count weights
        	for(StrandedBaseCount hit : ctrlHits){
        		int index = hit.getCoordinate()-currReg.getStart();
        		if(index>=0 && index<currReg.getWidth())
        			counts[index]+=hit.getCount();
        	}
        	
        	//Smooth
        	for(int i=0; i<currReg.getWidth(); i++){
        		double sum=0, num=0;
        		for(int s=i-(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s<i+(config.NOISE_DISTRIB_SMOOTHING_WIN/2); s++){
        			if(s>=0 && s<currReg.getWidth()){
        				num++; sum+=counts[s];
        			}
        		}
        		distrib[i] = sum/num;
        	}
        	
        	//Normalize
        	double total = 0;
        	for(int d=0; d<currReg.getWidth(); d++)
        		total+=distrib[d];
        	for(int d=0; d<currReg.getWidth(); d++)
        		distrib[d] = distrib[d]/total;
        	
        	return distrib;
        }
	}
	
	/**
	 * Consolidate and edit binding events.
	 * Follows a strict definition of binding event quantification - 
	 * if the event is not present in the condition, it gets a zero count assigned.
	 * Also merges positional duplicate events. 
	 * @param ev
	 * @return
	 */
	private List<BindingEvent> consolidateBindingEvents(List<BindingEvent> ev){
		List<BindingEvent> newEvents = new ArrayList<BindingEvent>();
		HashMap<Point, Integer> eventMap = new HashMap<Point, Integer>();
		int count=0;
		for(BindingEvent be : ev){
			if(!eventMap.containsKey(be.getPoint())){
				eventMap.put(be.getPoint(), count);
				count++;
				
				newEvents.add(be);
				//First time an event is added, clear out inactive events
				for(ExperimentCondition cond : manager.getConditions()){
					if(!be.isFoundInCondition(cond)){
						be.setCondSigHits(cond, 0.0);
		        		for(ControlledExperiment rep : cond.getReplicates()){
		        			be.setRepSigHits(rep, 0.0);
		        		}if(evconfig.CALC_EVENTS_LL)
			            	be.setLLd(cond, 0);
		        }	}
			}else{
				int index = eventMap.get(be.getPoint());
				//For events that are already in the list, just update the active condition read counts
				for(ExperimentCondition cond : manager.getConditions()){
					if(be.isFoundInCondition(cond)){
						newEvents.get(index).setCondSigHits(cond, be.getCondSigHits(cond));
						newEvents.get(index).setCondCtrlHits(cond, be.getCondCtrlHits(cond));
		        		for(ControlledExperiment rep : cond.getReplicates()){
		        			newEvents.get(index).setRepSigHits(rep, be.getRepSigHits(rep));
							newEvents.get(index).setRepCtrlHits(rep, be.getRepCtrlHits(rep));
		        		}if(evconfig.CALC_EVENTS_LL)
		        			newEvents.get(index).setLLd(cond, be.getLLd(cond));
		        }	}
			}
		}
		return newEvents;
	}
	
	/**
	 * ComponentConfiguration: represents a configuration of binding components as an array of positions.
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	protected class ComponentConfiguration{
		int [] positions=null;
		int parentCondIndex;
		//Constructor
		public ComponentConfiguration(List<BindingSubComponents> comps, int parentCondition){
			Collections.sort(comps);
			positions = new int[comps.size()];
			for(int p=0; p<comps.size(); p++)
				positions[p]=comps.get(p).getPosition();
			parentCondIndex = parentCondition;
		}
		//Return the index of the condition that initialized this configuration
		public int getParentCondition(){return parentCondIndex;}
		//Compare two configurations
		public boolean isSameAs(ComponentConfiguration cc){
			if(positions.length != cc.positions.length)
				return false;
			boolean isEqual=true;
			for(int p=0; p<positions.length; p++)
				isEqual = isEqual && positions[p]==cc.positions[p];
			return isEqual;
		}
	}
}
