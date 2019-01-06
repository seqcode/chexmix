package org.seqcode.projects.chexmix.framework;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.stats.BackgroundCollection;
import org.seqcode.deepseq.stats.PoissonBackgroundModel;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromosomeGenerator;
import org.seqcode.gseutils.RealValuedHistogram;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingModel;
import org.seqcode.projects.chexmix.events.EventsConfig;


/**
 * PotentialRegionFilter: Find a set of regions that are above a threshold in at least one condition. 
 * 		A region the size of the model span (i.e. 2x model range) potentially contains a binding site if 
 * 		it passes the Poisson threshold in at least one condition.
 * 		The Poisson thresholds are based on the model span size to keep consistent with the final used thresholds. 
 * Alternatively, pre-defined potential regions can be supplied as gff file. 
 * Overall counts for reads in potential regions and outside potential regions are maintained to assist noise model initialization.  
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class PotentialRegionFilter {

	protected ExperimentManager manager; 
	protected BindingManager bindingManager;
	protected EventsConfig evconfig;
	protected ChExMixConfig config;
	protected ExptConfig econfig;
	protected Genome gen;
	protected float maxBinWidth=0, binStep, winExt;
	protected boolean loadControl=true; 
	protected boolean stranded=false;
	protected List<Region> potentialRegions = new ArrayList<Region>();
	protected double potRegionLengthTotal=0;
	protected HashMap<ExperimentCondition, BackgroundCollection> conditionBackgrounds=new HashMap<ExperimentCondition, BackgroundCollection>(); //Background models for each replicate
	protected HashMap<ExperimentCondition, Double> potRegCountsSigChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> nonPotRegCountsSigChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> potRegCountsCtrlChannel = new HashMap<ExperimentCondition, Double>();
	protected HashMap<ExperimentCondition, Double> nonPotRegCountsCtrlChannel = new HashMap<ExperimentCondition, Double>();	
	protected HashMap<ControlledExperiment, Double> potRegCountsSigChannelByRep = new HashMap<ControlledExperiment, Double>();
	protected HashMap<ControlledExperiment, Double> nonPotRegCountsSigChannelByRep = new HashMap<ControlledExperiment, Double>();
	
	public PotentialRegionFilter(EventsConfig ec, ChExMixConfig c, ExptConfig econ, ExperimentManager eman, BindingManager bman){
		manager = eman;
		bindingManager = bman;
		evconfig = ec;
		config = c; 
		econfig = econ;
		gen = config.getGenome();
		maxBinWidth = c.getModelRange();
		//Initialize background models
		for(ExperimentCondition cond : manager.getConditions()){
			conditionBackgrounds.put(cond, new BackgroundCollection());
    			
    		//global threshold
    		conditionBackgrounds.get(cond).addBackgroundModel(new PoissonBackgroundModel(-1, config.getPRLogConf(), cond.getTotalSignalCount(), config.getGenome().getGenomeLength(), econfig.getMappableGenomeProp(), maxBinWidth, '.', 1, true));
    		//local windows won't work since we are testing per condition and we don't have a way to scale signal vs controls at the condition level (at least at this stage of execution)
    		
    		double thres = conditionBackgrounds.get(cond).getGenomicModelThreshold();
    		System.err.println("PotentialRegionFilter: genomic threshold for "+cond.getName()+" with bin width "+maxBinWidth+" = "+thres);
    			
    		//Initialize counts
    		potRegCountsSigChannel.put(cond, 0.0);
    		nonPotRegCountsSigChannel.put(cond, 0.0);
    		potRegCountsCtrlChannel.put(cond, 0.0);
    		nonPotRegCountsCtrlChannel.put(cond, 0.0);
    		for(ControlledExperiment rep : cond.getReplicates()){
    			potRegCountsSigChannelByRep.put(rep, 0.0);
        		nonPotRegCountsSigChannelByRep.put(rep, 0.0);
    		}
    	}
		binStep = config.POTREG_BIN_STEP;
		if(binStep>maxBinWidth/2)
			binStep=maxBinWidth/2;
		winExt = maxBinWidth/2;
	}
	
	//Accessors for read counts
	public Double getPotRegCountsSigChannel(ExperimentCondition e){ return potRegCountsSigChannel.get(e);}
	public Double getNonPotRegCountsSigChannel(ExperimentCondition e){ return nonPotRegCountsSigChannel.get(e);}
	public Double getPotRegCountsCtrlChannel(ExperimentCondition e){ return potRegCountsCtrlChannel.get(e);}
	public Double getNonPotRegCountsCtrlChannel(ExperimentCondition e){ return nonPotRegCountsCtrlChannel.get(e);}
	public Double getPotRegCountsSigChannelByRep(ControlledExperiment e){ return potRegCountsSigChannelByRep.get(e);}
	public Double getNonPotRegCountsSigChannelByRep(ControlledExperiment e){ return nonPotRegCountsSigChannelByRep.get(e);}
	public List<Region> getPotentialRegions(){return potentialRegions;}
	public List<Region> getSeqPotentialRegions(){return potentialRegions;}
	public double getPotRegionLengthTotal(){return potRegionLengthTotal;}
	
	/**
	 * Find list of potentially enriched regions 
	 * (windows that contain the minimum number of reads needed to pass the Poisson backgrounds).
	 * @param testRegions
	 */
	public List<Region> execute(){
		//TODO: check config for defined subset of regions
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(config.getGenome());
		
		//Threading divides analysis over entire chromosomes. This approach is not compatible with file caching. 
		int numThreads = econfig.getCacheAllData() ? config.getMaxThreads() : 1;
				
		Thread[] threads = new Thread[numThreads];
        ArrayList<Region> threadRegions[] = new ArrayList[numThreads];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }i=0;
        while(testRegions.hasNext()){
        	Region r = testRegions.next(); 
            threadRegions[(i++) % numThreads].add(r);
        }

        for (i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new PotentialRegionFinderThread(threadRegions[i]));
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
        
        //Initialize signal & noise counts based on potential region calls
        for(ExperimentCondition cond : manager.getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			if(rep.getSignalVsNoiseFraction()==0) //Only update if not already initialized
    				rep.setSignalVsNoiseFraction(potRegCountsSigChannelByRep.get(rep)/(potRegCountsSigChannelByRep.get(rep)+nonPotRegCountsSigChannelByRep.get(rep)));
    		}
        }
        
        for(Region r : potentialRegions)
        	potRegionLengthTotal+=(double)r.getWidth();
        
     	return potentialRegions;
	}
	
	/**
	 * Use potential regions provided by a user as a gff file
	 * @param testRegions
	 */
	public List<Region> executeInitialPositionFiler(){
				
		List<Region> geneTrackRegions = new ArrayList<Region>();
		// List of potential regions defined by running gene track
		for (Point p : config.getInitialPos()){
			geneTrackRegions.add(p.expand((int) config.getWindowExtension()/2));
		}
		// sort and merge regions
		List<Region> meregedRegion = Region.mergeRegions(geneTrackRegions);
		
		// Regions should not start with zero because it produces an error in sequence generation
		List<Region> chrStartExcluded = new ArrayList<Region>();
		for (Region r : meregedRegion){
			if (r.getStart()==0)
				chrStartExcluded.add(new Region(r.getGenome(), r.getChrom(),1, r.getEnd()));
			else
				chrStartExcluded.add(r);
		}
		
		potentialRegions = filterExcluded(chrStartExcluded);
		
		// signal and control counts from potential regions
		countReadsInRegionsNoLandscape(potentialRegions);
		countReadsInRegionsByRepNoLandscape(potentialRegions);
		
		//Initialize signal & noise counts based on potential region calls
        for(ExperimentCondition cond : manager.getConditions()){
    		for(ControlledExperiment rep : cond.getReplicates()){
    			if(rep.getSignalVsNoiseFraction()==0) //Only update if not already initialized
    				rep.setSignalVsNoiseFraction(potRegCountsSigChannelByRep.get(rep)/(potRegCountsSigChannelByRep.get(rep)+nonPotRegCountsSigChannelByRep.get(rep)));
    		}
        }
		
		// total length of potential regions
        for(Region r : potentialRegions)
        	potRegionLengthTotal+=(double)r.getWidth();
		
        //Note: it looks like currPotRegions and threadPotentials are redundant in the above, but they are not.
        //currPotRegions is only used to count sig/noise reads in the current section. threadPotentials stores regions over the entire run. 	
		return potentialRegions;
	}
	
	/**
	 * Count the total reads within potential regions (by condition). Assumes regions are sorted. This does not make a hit landscape.
	 * We  ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
	 * (and that thus will not be accounted for by the global noise model)  
	 * @param regs
	 */
	protected void countReadsInRegionsNoLandscape(List<Region> regs){
		List<Region> notPotentials = new ArrayList<Region>();
		String currChrm = null; Region prevReg = null ;	
		//Iterate through regions
		for (Region  reg : regs){
			Region notPotRegionPrevChrom = null; Region notPotRegionCurrChrom = null;
			if (!reg.getChrom().equals(currChrm)){	//new chromosome
				currChrm = reg.getChrom().toString();	//update chromosome only when you encounter the new chromosome
				if (prevReg == null){	// 1st iteration
					if (reg.getStart()!=1) //if region starts is not beginning of chromosome
						notPotRegionCurrChrom = new Region(gen,reg.getChrom(),1, reg.getStart()-1);
				}else{
					if (prevReg.getEnd()!=gen.getChromLength(prevReg.getChrom())) //if the prevReg end is not chromosome end
						notPotRegionPrevChrom = new Region(gen,prevReg.getChrom(),prevReg.getEnd()+1, gen.getChromLength(prevReg.getChrom()));
					if (reg.getStart()!=1) //if region starts is not beginning of chromosome
						notPotRegionCurrChrom = new Region(gen,reg.getChrom(),1, reg.getStart()-1);
				}
			}else{	//same chromosome
				if (reg.getStart() > prevReg.getEnd()+1) // potential regions are not non-adjacent
					notPotRegionCurrChrom = new Region(gen,reg.getChrom(),prevReg.getEnd()+1, reg.getStart()-1);				
			}
			// add regions to the list
			if (notPotRegionPrevChrom!=null)
				notPotentials.add(notPotRegionPrevChrom);
			if (notPotRegionCurrChrom!=null)
				notPotentials.add(notPotRegionCurrChrom);
				
			// update previous region
			prevReg = reg;
		}
		if (prevReg.getEnd()!=gen.getChromLength(prevReg.getChrom()))
			notPotentials.add(new Region(gen,prevReg.getChrom(),prevReg.getEnd()+1, gen.getChromLength(prevReg.getChrom())));
			
		//Iterate through experiments
		for(ExperimentCondition cond : manager.getConditions()){	
			double currPotWeightSig=0, currNonPotWeightSig=0, currPotWeightCtrl=0, currNonPotWeightCtrl=0;
			// count hits of notPotRegion and potRegion from signal and control experiments
			for (ControlledExperiment rep : cond.getReplicates()){
				for (Region  reg : regs){
					currPotWeightSig += rep.getSignal().countHits(reg);
					if (rep.hasControl())
						currPotWeightCtrl += rep.getControl().countHits(reg);
				}
				for (Region reg : notPotentials){
					currNonPotWeightSig += rep.getSignal().countHits(reg);
					if (rep.hasControl())
						currNonPotWeightCtrl += rep.getControl().countHits(reg);
				}
			}
			
			potRegCountsSigChannel.put(cond, currPotWeightSig);
			nonPotRegCountsSigChannel.put(cond, currNonPotWeightSig);
			potRegCountsCtrlChannel.put(cond, currPotWeightCtrl);
			nonPotRegCountsCtrlChannel.put(cond, currNonPotWeightCtrl);
		}
	}
	/**
	 * Count the total reads within potential regions (by replicate). Assumes regions are sorted. This does not make a hit landscape.
	 * We  ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
	 * (and that thus will not be accounted for by the global noise model)  
	 * @param regs
	 */
	protected void countReadsInRegionsByRepNoLandscape(List<Region> regs){
		List<Region> notPotentials = new ArrayList<Region>();
		String currChrm = null; Region prevReg = null ;	
		//Iterate through regions
		for (Region  reg : regs){
			Region notPotRegionPrevChrom = null; Region notPotRegionCurrChrom = null;
			if (!reg.getChrom().equals(currChrm)){	//new chromosome
				currChrm = reg.getChrom().toString();	//update chromosome only when you encounter the new chromosome
				if (prevReg == null){	// 1st iteration
					if (reg.getStart()!=1) //if region starts is not beginning of chromosome
						notPotRegionCurrChrom = new Region(gen,reg.getChrom(),1, reg.getStart()-1);
				}else{
					if (prevReg.getEnd()!=gen.getChromLength(prevReg.getChrom())) //if the prevReg end is not chromosome end
						notPotRegionPrevChrom = new Region(gen,prevReg.getChrom(),prevReg.getEnd()+1, gen.getChromLength(prevReg.getChrom()));
					if (reg.getStart()!=1) //if region starts is not beginning of chromosome
						notPotRegionCurrChrom = new Region(gen,reg.getChrom(),1, reg.getStart()-1);
				}
			}else{	//same chromosome
				if (reg.getStart() > prevReg.getEnd()+1) // potential regions are not non-adjacent
					notPotRegionCurrChrom = new Region(gen,reg.getChrom(),prevReg.getEnd()+1, reg.getStart()-1);				
			}
			// add regions to the list
			if (notPotRegionPrevChrom!=null)
				notPotentials.add(notPotRegionPrevChrom);
			if (notPotRegionCurrChrom!=null)
				notPotentials.add(notPotRegionCurrChrom);
				
			// update previous region
			prevReg = reg;
		}
		if (prevReg.getEnd()!=gen.getChromLength(prevReg.getChrom()))
			notPotentials.add(new Region(gen,prevReg.getChrom(),prevReg.getEnd()+1, gen.getChromLength(prevReg.getChrom())));
				
		//Iterate through experiments
		for(ControlledExperiment rep : manager.getReplicates()){
			double currPotWeightSig=0, currNonPotWeightSig=0;							
			// count hits of notPotRegion and potRegion from signal and control experiments				
			for (Region  reg : regs)
				currPotWeightSig += rep.getSignal().countHits(reg);	
			for (Region reg : notPotentials){
				currNonPotWeightSig += rep.getSignal().countHits(reg);			
			potRegCountsSigChannelByRep.put(rep, currPotWeightSig);
			nonPotRegCountsSigChannelByRep.put(rep, currNonPotWeightSig);
			}
		}
	}	
	
	//Filter out pre-defined regions to ignore (e.g. tower regions)
    protected List<Region> filterExcluded(List<Region> testRegions) {
		List<Region> filtered = new ArrayList<Region>();
		if(config.getRegionsToIgnore().size()==0)
			return testRegions;
		
		for(Region t : testRegions){
			boolean ignore = false;
			int x=0;
			while(x<config.getRegionsToIgnore().size() && ignore==false){
				Region i = config.getRegionsToIgnore().get(x);
				if(t.overlaps(i))
					ignore = true;
			}
			if(!ignore)
				filtered.add(t);
		}
		return filtered;
	}
	
	/**
     * Print potential regions to a file.
     * TESTING ONLY 
     */
    public void printPotentialRegionsToFile(){
    	try {
    		String filename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".potential.regions";
			FileWriter fout = new FileWriter(filename);
			for(Region r : potentialRegions){
	    		fout.write(r.getLocationString()+"\n");			
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
	
    class PotentialRegionFinderThread implements Runnable {
        private Collection<Region> regions;
        private float[][] landscape=null;
        private float[][] starts=null;
        private List<Region> threadPotentials = new ArrayList<Region>();
        
        public PotentialRegionFinderThread(Collection<Region> r) {
            regions = r;
        }
        
        public void run() {
        	int expansion = (int)(winExt + maxBinWidth/2);
        	for (Region currentRegion : regions) {
            	Region lastPotential=null;
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
                    int y = (int) (x+config.MAXSECTION+(expansion)); //Leave a little overhang to handle enriched regions that may hit the border. Since lastPotential is defined above, a region on the boundary should get merged in.
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
                    
                    List<Region> currPotRegions = new ArrayList<Region>();
                    List<List<StrandedBaseCount>> ipHits = new ArrayList<List<StrandedBaseCount>>();
                    List<List<StrandedBaseCount>> backHits = new ArrayList<List<StrandedBaseCount>>();
                    List<List<StrandedBaseCount>> ipHitsByRep = new ArrayList<List<StrandedBaseCount>>();
                    
                    synchronized(manager){
	                    //Initialize the read lists
                    	for(ExperimentCondition cond : manager.getConditions()){
                    		ipHits.add(new ArrayList<StrandedBaseCount>());
                			backHits.add(new ArrayList<StrandedBaseCount>());
                    		for(ControlledExperiment rep : cond.getReplicates())
                    			ipHitsByRep.add(new ArrayList<StrandedBaseCount>());
                    	}
                    	//Load signal reads by condition and by replicate, so that signal proportion estimates can be assigned to each replicate 
                    	for(ExperimentCondition cond : manager.getConditions()){
                    		for(ControlledExperiment rep : cond.getReplicates()){
                    			ipHits.get(cond.getIndex()).addAll(rep.getSignal().getBases(currSubRegion));
                    			ipHitsByRep.get(rep.getIndex()).addAll(rep.getSignal().getBases(currSubRegion));
                    		}for(Sample ctrl : cond.getControlSamples())
                    			backHits.get(cond.getIndex()).addAll(ctrl.getBases(currSubRegion));
                    		Collections.sort(ipHits.get(cond.getIndex()));
                    		Collections.sort(backHits.get(cond.getIndex()));
                    	}
                    }
                    
            		int numStrandIter = stranded ? 2 : 1;
                    for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
                        //If stranded peak-finding, run over both strands separately
                        char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					 
                        makeHitLandscape(ipHits, currSubRegion, maxBinWidth, binStep, str);
                        float ipHitCounts[][] = landscape.clone();
                        float ipBinnedStarts[][] = starts.clone();
                        float backBinnedStarts[][] = null;
                        if (loadControl) {
                            makeHitLandscape(backHits, currSubRegion, maxBinWidth, binStep, str);
                            backBinnedStarts = starts.clone();
                        }
					
                        //Scan regions
                        int currBin=0;
                        for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)maxBinWidth; i+=(int)binStep){
                        	boolean regionPasses=false;
                        	for(ExperimentCondition cond : manager.getConditions()){
                        		double ipWinHits=ipHitCounts[cond.getIndex()][currBin];
                        		//First Test: is the read count above the genome-wide thresholds?
                        		//If there is a fixed alpha, we should use that as the only threshold
                        		if(config.getFixedAlpha()>0){
                        			if(ipWinHits>config.getFixedAlpha()){
                        				regionPasses=true;
                        				break;
                        			}
                        		}else if(conditionBackgrounds.get(cond).passesGenomicThreshold((int)ipWinHits, str)){
                        			//Second Test: refresh all thresholds & test again
                        			conditionBackgrounds.get(cond).updateModels(currSubRegion, i-x, ipBinnedStarts[cond.getIndex()], backBinnedStarts==null ? null : backBinnedStarts[cond.getIndex()], binStep);
                        			if(conditionBackgrounds.get(cond).passesAllThresholds((int)ipWinHits, str)){
                        				//If the region passes the thresholds for one condition, it's a potential
                        				regionPasses=true;
                        				break;
		                            }
		                        }
                        	}
                        	if(regionPasses){
                        		Region currPotential = new Region(gen, currentRegion.getChrom(), Math.max(i-expansion, 1), Math.min((int)(i-1+expansion), currentRegion.getEnd()));
                        		if(lastPotential!=null && currPotential.overlaps(lastPotential)){
                        			lastPotential = lastPotential.expand(0, currPotential.getEnd()-lastPotential.getEnd());
                        		}else{
                        			//Add the last recorded region to the list
                        			if(lastPotential!=null){
                        				if(lastPotential.getWidth()<=config.getBMAnalysisWindowMax()){
                        					currPotRegions.add(lastPotential);
                        					threadPotentials.add(lastPotential);
                        				}else{
                        					//Break up long windows
                        					List<Region> parts = breakWindow(lastPotential, ipHits, config.getBMAnalysisWindowMax(), str);
                        					for(Region p : parts){
                        						currPotRegions.add(p);
                        						threadPotentials.add(p);
                        					}
                        				}
                        			}lastPotential = currPotential;
                        		}
                        	}
                            currBin++;
                        }
					}
                    //Count all "signal" reads overlapping the regions in currPotRegions (including the lastPotential)
                    if(lastPotential!=null)
                    	currPotRegions.add(lastPotential);
                    currPotRegions = filterExcluded(currPotRegions);
                    countReadsInRegions(currPotRegions, ipHits, backHits, y==currentRegion.getEnd() ? y : y-expansion);
                    countReadsInRegionsByRep(currPotRegions, ipHitsByRep, y==currentRegion.getEnd() ? y : y-expansion);
                    //Note: it looks like currPotRegions and threadPotentials are redundant in the above, but they are not.
                    //currPotRegions is only used to count sig/noise reads in the current section. threadPotentials stores regions over the entire run.
                }
                //Add the final recorded region to the list
                if(lastPotential!=null){
                	if(lastPotential.getWidth()<=config.getBMAnalysisWindowMax()){
                		threadPotentials.add(lastPotential);
                	}else{
                		//Break up long windows
                		List<List<StrandedBaseCount>> ipHits = new ArrayList<List<StrandedBaseCount>>();
                		for(ExperimentCondition cond : manager.getConditions())
                		    ipHits.add(new ArrayList<StrandedBaseCount>());
                		for(ExperimentCondition cond : manager.getConditions())
                			for(ControlledExperiment rep : cond.getReplicates())
                				ipHits.get(cond.getIndex()).addAll(rep.getSignal().getBases(lastPotential));
    					List<Region> parts = breakWindow(lastPotential, ipHits, config.getBMAnalysisWindowMax(), '.');
    					for(Region p : parts)
    						threadPotentials.add(p);
                	}
                }
                threadPotentials = filterExcluded(threadPotentials);
            }
        	if(threadPotentials.size()>0){
        		synchronized(potentialRegions){
        			potentialRegions.addAll(threadPotentials);
        		}
        	}	
        }
        
        //Break up a long window into parts
        //For now, we just choose the break points as the bins with the lowest total signal read count around the desired length.
        //TODO: improve?
        protected List<Region> breakWindow(Region lastPotential, List<List<StrandedBaseCount>> ipHits, int preferredWinLen, char str) {
			List<Region> parts = new ArrayList<Region>();
			makeHitLandscape(ipHits, lastPotential, maxBinWidth, binStep, str);
            float ipHitCounts[][] = landscape.clone();
            
            int currPartStart = lastPotential.getStart();
            double currPartTotalMin=Double.MAX_VALUE; int currPartTotalMinPos = -1;
            int currBin=0;
            for(int i=lastPotential.getStart(); i<lastPotential.getEnd()-(int)maxBinWidth; i+=(int)binStep){
            	if(lastPotential.getEnd()-currPartStart < (preferredWinLen*1.5))
            		break;
            	float currBinTotal=0;
            	for(ExperimentCondition cond : manager.getConditions())
                	currBinTotal+=ipHitCounts[cond.getIndex()][currBin];
            	
            	if(i>(currPartStart+preferredWinLen-1000) && i<(currPartStart+preferredWinLen+1000)){ 
            		if(currBinTotal<currPartTotalMin){
            			currPartTotalMin=currBinTotal;
            			currPartTotalMinPos=i;
            		}
            	}
            	//Add a new part
            	if(i>=(currPartStart+preferredWinLen+1000)){
            		parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, currPartTotalMinPos));
            		currPartStart = currPartTotalMinPos+1;
            		currPartTotalMin=Double.MAX_VALUE; currPartTotalMinPos = -1;
            	}
            	currBin++;
            }
            parts.add(new Region(lastPotential.getGenome(), lastPotential.getChrom(), currPartStart, lastPotential.getEnd()));
            
			return parts;
		}

		//Makes integer arrays corresponding to the read landscape over the current region.
        //Reads are semi-extended out to bin width to account for the bin step
        //No needlefiltering here as that is taken care of during read loading (i.e. in Sample)
    	protected void makeHitLandscape(List<List<StrandedBaseCount>> hits, Region currReg, float binWidth, float binStep, char strand){
    		int numBins = (int)(currReg.getWidth()/binStep);
    		landscape = new float[hits.size()][numBins+1];
    		starts = new float[hits.size()][numBins+1];
    		float halfWidth = binWidth/2;

    		for(ExperimentCondition cond : manager.getConditions()){
        		List<StrandedBaseCount> currHits = hits.get(cond.getIndex());
    			for(int i=0; i<=numBins; i++){landscape[cond.getIndex()][i]=0; starts[cond.getIndex()][i]=0; }
	    		for(StrandedBaseCount r : currHits){
	    			if(strand=='.' || r.getStrand()==strand){
	    				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
	    				int binoff = inBounds((int)(offset/binStep), 0, numBins);
	    				starts[cond.getIndex()][binoff]+=r.getCount();
	    				int binstart = inBounds((int)((double)(offset-halfWidth)/binStep), 0, numBins);
	    				int binend = inBounds((int)((double)(offset)/binStep), 0, numBins);
	    				for(int b=binstart; b<=binend; b++)
	    					landscape[cond.getIndex()][b]+=r.getCount();
	    			}
            	}
    		}
    	}
    	protected final int inBounds(int x, int min, int max){
    		if(x<min){return min;}
    		if(x>max){return max;}
    		return x;
    	}
    	
    	/**
    	 * Count the total reads within potential regions via semi binary search (by condition).
    	 * Assumes both regs and ipHits are sorted.
    	 * We don't have to check chr String matches, as the hits were extracted from the chromosome
    	 * EndCoord accounts for the extra overhang added to some wide regions
    	 * We also ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
    	 * (and that thus will not be accounted for by the global noise model)  
    	 * @param regs
    	 * @param ipHits
    	 * @param ctrlHits
    	 * @param endCoord
    	 */
    	protected void countReadsInRegions(List<Region> regs, List<List<StrandedBaseCount>> ipHits, List<List<StrandedBaseCount>> ctrlHits, int endCoord){
    		//Iterate through experiments
    		for(ExperimentCondition cond : manager.getConditions()){
    			double currPotWeightSig=0, currNonPotWeightSig=0, currPotWeightCtrl=0, currNonPotWeightCtrl=0;
    			//Iterate through signal hits
    			for(StrandedBaseCount hit : ipHits.get(cond.getIndex())){
    				if(regs.size()==0)
    					currNonPotWeightSig+=hit.getCount();
    				else{
    					//Binary search for closest region start
        				int hpoint = hit.getCoordinate();
        				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightSig+=hit.getCount(); inPot=true; break;
	        	            	}
	        	            }
	        	            if(!inPot)
	        	            	currNonPotWeightSig+=hit.getCount();
        				}
    				}
    			}
    			//Iterate through control hits
    			for(StrandedBaseCount hit : ctrlHits.get(cond.getIndex())){
    				if(regs.size()==0)
    					currNonPotWeightCtrl+=hit.getCount();
    				else{
        				//Binary search for closest region start
        				int hpoint = hit.getCoordinate();
        				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightCtrl+=hit.getCount(); inPot=true; break;
	        	            	}
	        	            }
	        	            if(!inPot)
	        	            	currNonPotWeightCtrl+=hit.getCount();
        				}
    				}
    			}
    			synchronized(potRegCountsSigChannel){
    				potRegCountsSigChannel.put(cond, potRegCountsSigChannel.get(cond)+currPotWeightSig);
    			}
    			synchronized(nonPotRegCountsSigChannel){
    				nonPotRegCountsSigChannel.put(cond, nonPotRegCountsSigChannel.get(cond)+currNonPotWeightSig);
    			}
    			synchronized(potRegCountsCtrlChannel){
    				potRegCountsCtrlChannel.put(cond, potRegCountsCtrlChannel.get(cond)+currPotWeightCtrl);
    			}
    			synchronized(nonPotRegCountsCtrlChannel){
    				nonPotRegCountsCtrlChannel.put(cond, nonPotRegCountsCtrlChannel.get(cond)+currNonPotWeightCtrl);
    			}
    		}
	    }
    }
    
    /**
	 * Count the total reads within potential regions via semi binary search (by replicate).
	 * Assumes both regs and ipHitsByRep are sorted.
	 * We don't have to check chr String matches, as the hits were extracted from the chromosome
	 * EndCoord accounts for the extra overhang added to some wide regions
	 * We also ignore strandedness here -- the object is to count ALL reads that will be loaded for analysis later
	 * (and that thus will not be accounted for by the global noise model)  
	 * @param regs
	 * @param ipHitsByRep
	 * @param ctrlHits
	 * @param endCoord
	 */
	protected void countReadsInRegionsByRep(List<Region> regs, List<List<StrandedBaseCount>> ipHitsByRep, int endCoord){
		//Iterate through experiments
		for(ExperimentCondition cond : manager.getConditions()){
			for(ControlledExperiment rep : cond.getReplicates()){
				double currPotWeightSig=0, currNonPotWeightSig=0;
				//Iterate through signal hits
				for(StrandedBaseCount hit : ipHitsByRep.get(rep.getIndex())){
					if(regs.size()==0)
						currNonPotWeightSig+=hit.getCount();
					else{
						//Binary search for closest region start
	    				int hpoint = hit.getCoordinate();
	    				if(hpoint<endCoord){ //Throw this check in for the overhang
	        				int l = 0, r = regs.size()-1;
	        	            while (r - l > 1) {
	        	                int c = (l + r) / 2;
	        	                if (hpoint >= regs.get(c).getStart()) {
	        	                    l = c;
	        	                } else {
	        	                    r = c;
	        	                }
	        	            }
	        	            boolean inPot = false;
	        	            for(int x=l; x<=r; x++){
	        	            	if(hpoint >= regs.get(x).getStart() && hpoint <= regs.get(x).getEnd()){
	        	            		currPotWeightSig+=hit.getCount(); inPot=true; break;
	        	            	}
	        	            }
	        	            if(!inPot)
	        	            	currNonPotWeightSig+=hit.getCount();
	    				}
					}
				}
				
				synchronized(potRegCountsSigChannelByRep){
					potRegCountsSigChannelByRep.put(rep, potRegCountsSigChannelByRep.get(rep)+currPotWeightSig);
				}
				synchronized(nonPotRegCountsSigChannelByRep){
					nonPotRegCountsSigChannelByRep.put(rep, nonPotRegCountsSigChannelByRep.get(rep)+currNonPotWeightSig);
				}
			}
		}
    }

    
    /**
	 * This main method is only for testing the PotentialRegionFilter
	 * @param args
	 */
	public static void main(String[] args){
		
		System.setProperty("java.awt.headless", "true");
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		econfig.setPerBaseReadFiltering(false);	
		EventsConfig evconfig = new EventsConfig(gconfig, args);
		ChExMixConfig config = new ChExMixConfig(gconfig, args);
		if(config.helpWanted()){
			System.err.println("PotentialRegionFilter:");
			System.err.println(config.getArgsList());
		}else{
			ExperimentManager manager = new ExperimentManager(econfig);
			BindingManager bman = new BindingManager(evconfig, manager, config.isLenientMode());
			//Initialize binding models & binding model record
			Map<ControlledExperiment, List<BindingModel>> repBindingModels = new HashMap<ControlledExperiment, List<BindingModel>>();
			for(ControlledExperiment rep : manager.getReplicates()){
				bman.setUnstrandedBindingModel(rep, new BindingModel(BindingModel.defaultChipExoEmpiricalDistribution));
				repBindingModels.put(rep, new ArrayList<BindingModel>());
				repBindingModels.get(rep).add(bman.getUnstrandedBindingModel(rep));
			}		
						
			System.err.println("Conditions:\t"+manager.getConditions().size());
			for(ExperimentCondition c : manager.getConditions()){
				System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
			}
			for(ExperimentCondition c : manager.getConditions()){
				for(ControlledExperiment r : c.getReplicates()){
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
					if(r.getControl()==null)
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
					else
						System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
				}
			}
			
			// Either execute potential region filter or read from gene track file.
			PotentialRegionFilter filter = new PotentialRegionFilter(evconfig, config, econfig, manager, bman);
			List<Region> potentials= null;
			int stop=10000;
			if (config.getInitialPos()!=null || config.getInitialPos().size()>0){
				stop=1000;
				potentials = filter.executeInitialPositionFiler();
			}else{
				potentials = filter.execute();
			}
			
			RealValuedHistogram histo = new RealValuedHistogram(0, stop, 20);
			int min = Integer.MAX_VALUE;
			int max = -Integer.MAX_VALUE;
			for(Region r : potentials){
				System.out.println(r.getLocationString()+"\t"+r.getWidth());
				histo.addValue(r.getWidth());
				if(r.getWidth()<min)
					min = r.getWidth();
				if(r.getWidth()>max)
					max = r.getWidth();
			}
			System.err.println("Number of Potential Regions: "+potentials.size() +"\tLength: "+filter.getPotRegionLengthTotal());
			System.err.println("Min width: "+min+"\tMax width: "+max);
			histo.printContents(true);
			
			System.err.println("number of reads in potential regions");
			for(ExperimentCondition c : manager.getConditions()){
				System.err.println("Condition "+c.getName());
				System.err.println("PotRegCountsSigChannel: "+ filter.getPotRegCountsSigChannel(c));
				System.err.println("NonPotRegCountsSigChannel: "+ filter.getNonPotRegCountsSigChannel(c));
				System.err.println("PotRegCountsCtrlChannel: "+ filter.getPotRegCountsCtrlChannel(c));
				System.err.println("NonPotRegCountsCtrlChannel: "+ filter.getNonPotRegCountsCtrlChannel(c));
				for (ControlledExperiment r : c.getReplicates()){
					System.err.println("Rep "+r.getName() );
					System.err.println("PotRegCountsSigChannelByRep: "+ filter.getPotRegCountsSigChannelByRep(r));
					System.err.println("NonPotRegCountsSigChannelByRep: "+ filter.getNonPotRegCountsSigChannelByRep(r));
				}
			}
			
			System.err.println("Proportions of tags in potential regions");
			for(ExperimentCondition c : manager.getConditions())
				for(ControlledExperiment r : c.getReplicates())
					System.err.println("Condition "+c.getName()+":\tRep "+r.getName()+"\t"+r.getSignalVsNoiseFraction());
					
			manager.close();
		}
	}
}
