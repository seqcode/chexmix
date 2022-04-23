package org.seqcode.projects.chexmix.motifs;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.data.motifdb.WeightMatrixImport;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.sequence.RandomSequenceGenerator;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.motifs.DrawMotifs;
import org.seqcode.motifs.MarkovMotifThresholdFinder;
import org.seqcode.projects.chexmix.ChExMix;
import org.seqcode.projects.chexmix.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.composite.XLAnalysisConfig;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingSubtype;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.mixturemodel.BindingSubComponents;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.stats.InformationContent;


public class MotifPlatform {
	protected GenomeConfig gconfig;
	protected Genome genome;
	protected ChExMixConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected SequenceGenerator<Region> seqgen;
	protected List<Region> randomRegions = new ArrayList<Region>(); //Randomly chosen regions for motif significance tests
	protected String[] randomSequences; //Randomly chosen sequences for motif significance tests
	protected ArrayList<String> MarkovRandSeq = new ArrayList<String>(); // Random sequences for Markov Threshold
	protected MEMERunner meme;
	protected MarkovBackgroundModel backMod;
	protected SimpleMotifAligner aligner;
	protected List<Region> cachedRegions;
		
	/**
	 * Constructor for motif platform
	 * @param c
	 * @param man
	 * @param regionsOfInterest: this list contains all possible regions that motif-finding/scanning may be run on. 
	 */
	public MotifPlatform(GenomeConfig g, ChExMixConfig c, ExperimentManager man, BindingManager bman, List<Region> regionsOfInterest ){
		gconfig = g;
		genome = gconfig.getGenome();
		config = c;
		manager=man;
		bindingManager = bman;
		seqgen = gconfig.getSequenceGenerator();
		seqgen.useCache(true);
		if (!seqgen.isRegionCached()){
			System.err.println("Caching sequences");
			randomRegions = randomRegionPick(null, config.MOTIF_FINDING_NEGSEQ, config.MOTIF_FINDING_SEQWINDOW);
			randomSequences = seqgen.setupRegionCache(regionsOfInterest, randomRegions);
			System.err.println("Caching completed");
		}
		cachedRegions = regionsOfInterest;
		meme = new MEMERunner(config, man);
		
		backMod = c.getBackModel();
		
        //Pre-load the random sequences
		RandomSequenceGenerator randGen = new RandomSequenceGenerator(backMod, config.RANDOMSEED);
		for(int i=0; i<config.MARKOV_NUM_TEST; i++){
			MarkovRandSeq.add(randGen.execute(config.MOTIF_FINDING_SEQWINDOW));
		}
		
		aligner = new SimpleMotifAligner(config.getMinMotifLength());
	}
	
	/**
	 * Extract sequences around clustered points and call the MEME runner 
	 * @param cond
	 * @param activeComponents
	 * @param trainingRound
	 * @throws Exception 
	 */
	public List<List<List<StrandedPoint>>> findClusterMotifs(List<List<List<StrandedPoint>>> points, int trainingRound) throws Exception{
		
		List<List<List<StrandedPoint>>> adjPoints=new ArrayList<List<List<StrandedPoint>>>();				
		for (ExperimentCondition cond : manager.getConditions()){	
			List<List<StrandedPoint>> modelRefs = new ArrayList<List<StrandedPoint>>();
			int counter=0;
			for (List<StrandedPoint> clusterPoints : points.get(cond.getIndex())){
				List<StrandedPoint> newCenterPos = new ArrayList<StrandedPoint>();
				List<String> seqs = new ArrayList<String>();
				boolean motifFound =false;
				for (StrandedPoint p : clusterPoints){
					int start = p.getLocation()-config.MOTIF_FINDING_SEQWINDOW;
					int end = p.getLocation()+config.MOTIF_FINDING_SEQWINDOW;
					if (start >=1 && end <= p.getGenome().getChromLength(p.getChrom())){
						Region peakReg = new Region(p.getGenome(), p.getChrom(), start, end); 
						boolean containing=false;	// Filter out regions that are edge of cached sequences
						for (Region reg : cachedRegions){
							if (reg.contains(peakReg.expand(config.MOTIF_FINDING_LOCAL_SEQWINDOW/2, config.MOTIF_FINDING_LOCAL_SEQWINDOW/2))){
								containing=true; break;
							}
						} 
						if (containing)
							seqs.add(seqgen.execute(peakReg));
					}
				}
				
				WeightMatrix m = WeightMatrixImport.buildAlignedSequenceMatrix(seqs);	
				InformationContent ic =  new InformationContent(m, backMod);
				// Smooth profiles and get max position
				double[] sic = StatUtil.gaussianSmoother(ic.getMotifIC(),config.getGaussSmoothParam());
				double maxScore = 0.0;
				int maxPos = 0;
				for (int i=0; i < sic.length; i++){
					if (sic[i] > maxScore){
						maxScore=sic[i]; maxPos=i;
					}
				}
				
				// Get sequences for local MEME search at maximum information content
				List<String> localSeqs = new ArrayList<String>();
				for (String currSeq : seqs){
					String seq = currSeq.substring(Math.max(0, maxPos-config.MOTIF_FINDING_LOCAL_SEQWINDOW/2), Math.min(currSeq.length()-1, maxPos+config.MOTIF_FINDING_LOCAL_SEQWINDOW/2));
					if(lowercaseFraction(currSeq)<=config.MOTIF_FINDING_ALLOWED_REPETITIVE){localSeqs.add(seq);}	
				}
					
				if (seqs.size() < config.getMinRefsForBMUpdate()){
					System.err.println("cannot do motif finding due to too few regions ("+localSeqs.size()+"<"+config.getMinRefsForBMUpdate()+").");
				}else{						
					//Execute MEME
					Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(localSeqs, new String("motif_"+cond.getName()+"_t"+trainingRound+"_c"+counter), false);
					List<WeightMatrix> wm = matrices.car();
					List<WeightMatrix> fm = matrices.cdr();
					if(wm.size()>0){
						//Evaluate the significance of the discovered motifs
						int bestMotif=0;
						double rocScores[] = motifROCScores(wm, localSeqs, randomSequences);
						double maxRoc=0;
						for(int i=0; i<rocScores.length; i++)
							if(rocScores[i]>maxRoc){
								maxRoc = rocScores[i]; 
								bestMotif=i;
							}
						//Results summary
						if(config.isVerbose()){
							System.err.println("MEME results for: "+cond.getName());
							for(int w=0; w<fm.size(); w++){
								if(fm.get(w)!=null){
									System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
								}}}		
						
						if(maxRoc >= config.getMotifMinROC()){
							motifFound = true;
							WeightMatrix currMotif = wm.get(bestMotif);
							MarkovMotifThresholdFinder finder = new MarkovMotifThresholdFinder(currMotif, backMod, config.MARKOV_NUM_TEST);
							finder.setRandomSeq(MarkovRandSeq);
							double motifThres = finder.execute(config.MARKOV_BACK_MODEL_THRES);
							for (StrandedPoint p : clusterPoints){
								int center = p.getLocation()-config.MOTIF_FINDING_SEQWINDOW+maxPos;
								Region peakReg = new Region(p.getGenome(), p.getChrom(), center-config.MOTIF_FINDING_LOCAL_SEQWINDOW/2, center+config.MOTIF_FINDING_LOCAL_SEQWINDOW/2);
								StrandedPoint mPoint = getMotifPosition(currMotif, motifThres, peakReg);
								if (mPoint!=null)
									if (Math.abs(center - mPoint.getLocation())<=5)
										newCenterPos.add(mPoint);
							}}}}			
				if(motifFound)
					modelRefs.add(newCenterPos);
				else
					modelRefs.add(clusterPoints);	
				
				counter++;
				
				//testing only
				if(motifFound)
					printMotifsTest(cond, m, clusterPoints, newCenterPos, counter);				
			}	
			adjPoints.add(modelRefs);			
		}	
		
		return adjPoints;
	}
	
	public void printMotifsTest(ExperimentCondition cond, WeightMatrix m, List<StrandedPoint> clusterPoints, List<StrandedPoint> newCenterPos, int counter) throws Exception{
//		String motifLabel = WeightMatrix.getConsensus(m)+"_"+counter+", MEME";
		DrawMotifs.printMotifLogo(m, new File("before_alignment_motif_"+counter+".png"), 75, counter+", MEME"); 
		
		CompositeTagDistribution signalComposite = new CompositeTagDistribution(clusterPoints, cond, config.MAX_BINDINGMODEL_WIDTH,true);							
		TagProbabilityDensity model = new TagProbabilityDensity(signalComposite.getWinSize()-1);
		model.loadData(signalComposite.getCompositeWatson(), signalComposite.getCompositeCrick());
		model.printDensityToFile("before_alignment_motif_"+cond.getName()+"_"+counter);
		
		List<String> newSeqs=new ArrayList<String>();
		for (StrandedPoint p : newCenterPos){
			Region peakReg = new Region(p.getGenome(), p.getChrom(), p.getLocation()-config.MOTIF_FINDING_SEQWINDOW, p.getLocation()+config.MOTIF_FINDING_SEQWINDOW); 
			boolean containing=false;	// Filter out regions that are edge of cached sequences
			for (Region reg : cachedRegions){
				if (reg.contains(peakReg)){
					containing=true; break;
				}
			} 
			if (containing)
				newSeqs.add(seqgen.execute(peakReg));
		}
		WeightMatrix wm = WeightMatrixImport.buildAlignedSequenceMatrix(newSeqs);
//		motifLabel = WeightMatrix.getConsensus(wm)+"_"+counter+", MEME";
		DrawMotifs.printMotifLogo(wm, new File("after_alignment_motif_"+counter+".png"), 75, counter+", MEME"); 
		
		signalComposite = new CompositeTagDistribution(newCenterPos, cond, config.MAX_BINDINGMODEL_WIDTH,true);							
		model = new TagProbabilityDensity(signalComposite.getWinSize()-1);
		model.loadData(signalComposite.getCompositeWatson(), signalComposite.getCompositeCrick());
		model.printDensityToFile("after_alignment_motif_"+cond.getName()+"_"+counter);
	}	
	
	/**
	 * Extract sequences around top BindingComponents and call the MEME runner 
	 * Remove sequences with motifs and recursively do motif discovery
	 * @param cond
	 * @param activeComponents
	 * @param trainingRound
	 * @throws Exception 
	 */
	public void updateSubtypesUsingMotifs(HashMap<Region, List<List<BindingSubComponents>>> activeComponents, int trainingRound) throws Exception{
		Map<ExperimentCondition, List<BindingSubComponents>> allpeaks = new HashMap<ExperimentCondition, List<BindingSubComponents>>();
		for (ExperimentCondition cond : manager.getConditions())
			allpeaks.put(cond, new ArrayList<BindingSubComponents>());
		//Choose which components to include
		for(Region r : activeComponents.keySet()){
			for (ExperimentCondition cond : manager.getConditions()){
				float edge = Math.max(config.getModelRange()/2, config.MOTIF_FINDING_SEQWINDOW/2);
				List<BindingSubComponents> peaks = new ArrayList<BindingSubComponents>();
				for(BindingSubComponents bc : activeComponents.get(r).get(cond.getIndex())){
					//Component must not be at the edge of the region for sequence region cache
					if((bc.getPosition()-r.getStart()>edge) && (r.getEnd()-bc.getPosition()>edge)){
						// Arbitrary minimum read support for distribution updates
						if(bc.getSumResponsibility()>(config.getMinComponentReadFactorForBM()*bindingManager.getAlpha(cond)))
							peaks.add(bc);					
				}}
				allpeaks.get(cond).addAll(peaks);
		}}
		
		for (ExperimentCondition cond : manager.getConditions())
			System.out.println("Condition "+cond.getName()+" : Number of active components used for motif finding "+allpeaks.get(cond).size());
		
		// Group peaks either by overlap with other conditions or subtype probabilities
		Map<ExperimentCondition, List<List<BindingSubComponents>>> groupPeaks;
		if (trainingRound<=2)
			groupPeaks = groupPeaksByOverlap(allpeaks);
		else
			groupPeaks = groupPeaksBySubtypes(allpeaks);		
		
		for (ExperimentCondition cond : manager.getConditions()){			
			// List of motifs and motif references to return for each condition
			List<WeightMatrix> selecMotifs = new ArrayList<WeightMatrix>();
			List<WeightMatrix> selecFreqMatrix = new ArrayList<WeightMatrix>();
			List<Set<StrandedPoint>> motifRefs = new ArrayList<Set<StrandedPoint>>();
			List<List<BindingSubComponents>> groupComps = new ArrayList<List<BindingSubComponents>>();
			int setCounter=0;
			
			for (List<BindingSubComponents> peaks : groupPeaks.get(cond)){	
				boolean groupFoundMotif = false;
				// Region list that gets removed as it finds motifs
				List<Region> sortedRegions = new ArrayList<Region>();
				// Iterate all regions to create hashMap for region to sequences
				Map<Region, String> reg2seq = new HashMap<Region, String>();
				for (BindingSubComponents bc : peaks){
					Region peakReg = new Region(bc.getCoord().getGenome(), bc.getCoord().getChrom(), bc.getPosition()-config.MOTIF_FINDING_SEQWINDOW/2, bc.getPosition()+config.MOTIF_FINDING_SEQWINDOW/2); 
					sortedRegions.add(peakReg);
					String seq = seqgen.execute(peakReg);
					reg2seq.put(peakReg, seq);
				}
		
				// find motifs recursively until no motifs pass ROC threshold scores
				boolean findMotif = true;	
				if (config.getInitialMotifs()!=null && (trainingRound<2 || !config.updateInitialMotifs())){
					findMotif=false;	// No motif finding if initial motifs are defined (first round of training only)
					System.out.println("Skipping motif-finding for this round of training, as initial motifs have been defined.");
				}else if (sortedRegions.size() < config.getMinRefsForBMUpdate()){
					findMotif = false;
					System.err.println("cannot do motif finding due to too few regions ("+sortedRegions.size()+"<"+config.getMinRefsForBMUpdate()+").");
				}
				
				int counter = 0;
				while (findMotif) {		
					//check size of region 
					System.out.println("number of regions "+sortedRegions.size());
        		
					//Get the top sequences that don't have too many lowercase letters (repeats)
					List<String> seqs = new ArrayList<String>();
					int addedSeqs=0;
					for (int r=0; r < sortedRegions.size() && addedSeqs< config.MOTIF_FINDING_TOPSEQS; r++){
						String currSeq = reg2seq.get(sortedRegions.get(r));
						if(lowercaseFraction(currSeq)<=config.MOTIF_FINDING_ALLOWED_REPETITIVE){
							seqs.add(currSeq);
							addedSeqs++;
						}
					}        		
					System.out.println("number of seq added "+addedSeqs);
			
					//Execute MEME
					Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqs, new String("motif_"+cond.getName()+"_t"+trainingRound+"_s"+setCounter+"_c"+counter), false);   	
					List<WeightMatrix> wm = matrices.car();
					List<WeightMatrix> fm = matrices.cdr();
		
					if(wm.size()>0){
						//Evaluate the significance of the discovered motifs
						double rocScores[] = motifROCScores(wm, seqs, randomSequences);
        			
						//Results summary
						if(config.isVerbose()){
							System.err.println("MEME results for: "+cond.getName());
							for(int w=0; w<fm.size(); w++){
								if(fm.get(w)!=null){
									System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
								}}}
			
						//Pick the condition's motif if the ROC is above threshold
						List<WeightMatrix> sigMotifs = new ArrayList<WeightMatrix>();
						List<WeightMatrix> sigFreqMatrix = new ArrayList<WeightMatrix>();						
						for(int i=0; i<rocScores.length; i++){
							if(rocScores[i]>= config.getMotifMinROC()){
								sigMotifs.add(wm.get(i));
								sigFreqMatrix.add(fm.get(i));
								System.out.println("adding motif "+wm.get(i).getName()+"\tROC:"+String.format("%.2f",rocScores[i])+" above "+config.getMotifMinROC());
							}}
						
						for (int m=0; m< sigMotifs.size(); m++){
							WeightMatrix currMotif = sigMotifs.get(m);
							WeightMatrix currFreqMatrix = sigFreqMatrix.get(m);	
							MarkovMotifThresholdFinder finder = new MarkovMotifThresholdFinder(currMotif, backMod, config.MARKOV_NUM_TEST);
							finder.setRandomSeq(MarkovRandSeq);
							Set<StrandedPoint> refs = new HashSet<StrandedPoint>();	// motif references
							List<Region> found = new ArrayList<Region>();	//region with motifs to be removed from next motif finding iteration
							double motifThres = finder.execute(config.MARKOV_BACK_MODEL_THRES); 
							double seqRmThres = finder.execute(config.getMarkovBackSeqRmThres()); 
							for (Region reg : sortedRegions){
								if (getMotifPosition(currMotif, motifThres, reg)!=null)
									refs.add(getMotifPosition(currMotif, motifThres, reg));
								if (getMotifPosition(currMotif, seqRmThres, reg)!=null)
									found.add(reg);
							}
							
							System.out.println("removing "+found.size()+" sequences with "+WeightMatrix.getConsensus(currFreqMatrix));
							
							if (refs.size() > config.getMinRefsForBMUpdate()){
								//remove all regions that has motif hits
								sortedRegions.removeAll(found);
								//Add motif and frequency matrix to list
								selecMotifs.add(currMotif);
								selecFreqMatrix.add(currFreqMatrix);
								motifRefs.add(refs);
								groupFoundMotif=true;	//Significant motif found in current group
							}else{
								findMotif=false;
							}
						} // End of single MEME output evaluation
						if (sigMotifs.size() == 0 || sortedRegions.size() < config.getMinRefsForBMUpdate())
							findMotif=false;
						counter++; //count number of MEME iteration
					}else{
						findMotif = false;
					}
				} // End of while loop to recursively find motifs
				
				// Define subtypes and read distributions using the provided motifs
				if (config.getInitialMotifs()!=null && (trainingRound<2 || !config.updateInitialMotifs())){
					//Define subtypes using initially provided motifs
					List<WeightMatrix> sigMotifs = config.getInitialMotifs();
					List<WeightMatrix> sigFreqMatrix = new ArrayList<WeightMatrix>();
					for (WeightMatrix m : sigMotifs){
						float[][] matrix = new float[m.matrix.length][WeightMatrix.MAXLETTERVAL];
						for (int i=0; i < m.matrix.length; i++)
							for (int j=0; j < WeightMatrix.MAXLETTERVAL; j++)
								matrix[i][j]=m.matrix[i][j];
						WeightMatrix wm = new WeightMatrix(matrix);
						wm.toFrequency();
						sigFreqMatrix.add(wm);
					}
				
					for (int m=0; m< sigMotifs.size(); m++){
						WeightMatrix currMotif = sigMotifs.get(m);
						WeightMatrix currFreqMatrix = sigFreqMatrix.get(m);	
						MarkovMotifThresholdFinder finder = new MarkovMotifThresholdFinder(currMotif, backMod, config.MARKOV_NUM_TEST);
						finder.setRandomSeq(MarkovRandSeq);
						Set<StrandedPoint> refs = new HashSet<StrandedPoint>();	// motif references
						double motifThres = finder.execute(config.MARKOV_BACK_MODEL_THRES); 
						for (Region reg : sortedRegions)
							if (getMotifPosition(currMotif, motifThres, reg)!=null)
								refs.add(getMotifPosition(currMotif, motifThres, reg));
					
						if (refs.size() > config.getMinRefsForBMUpdate()){
							//Add motif and frequency matrix to list
							selecMotifs.add(currMotif);
							selecFreqMatrix.add(currFreqMatrix);
							motifRefs.add(refs);
						}
					} 
				} // End of defining the motif positions using the known motifs
			
				
				// Keep list of peaks that can be used for BM updates later
				if (!groupFoundMotif){
					if (peaks.size() < config.getMinComponentsForBMUpdate())
						System.err.println("The "+cond.getName()+" read distributions cannot be updated due to too few binding components ("+peaks.size()+"<"+config.getMinComponentsForBMUpdate()+").");
					else
						groupComps.add(peaks);	
				}
				setCounter++; // Next set of sequences
			} // End of iteration grouped sequence from a condition
			
			// Make binding subtypes from motif references
			List<BindingSubtype> subtypes = new ArrayList<BindingSubtype>();
			for (int m=0; m< selecMotifs.size(); m++){
				BindingSubtype subtype = new BindingSubtype(cond, new ArrayList<StrandedPoint>(motifRefs.get(m)), config.MAX_BINDINGMODEL_WIDTH);
				subtype.setMotif(selecMotifs.get(m), selecFreqMatrix.get(m));
				subtypes.add(subtype);
			}
			
			// Make binding subtypes from binding events
			for (List<BindingSubComponents> comps : groupComps ){
				List<StrandedPoint> points = new ArrayList<StrandedPoint>();
				for (BindingSubComponents bc : comps){
					StrandedPoint p=null;
					if (bc.isSubtype())
						p = new StrandedPoint(bc.getCoord().getGenome(),bc.getCoord().getChrom(),bc.getPosition(),bc.getMaxTypeStrand());
					else
						p = new StrandedPoint(bc.getCoord().getGenome(),bc.getCoord().getChrom(),bc.getPosition(),'+');
					points.add(p);
				}	
				subtypes.add(new BindingSubtype(cond,points, config.MAX_BINDINGMODEL_WIDTH));
				System.err.println("Made new read distribution from " + points.size() +" binding events.");
			}
			
			// Print bindingModel
			int subtypeC=0;
			for (BindingSubtype subtype : subtypes){
				String distribFile = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_t"+trainingRound+"_ReadDistrib_"+cond.getName()+"_"+subtypeC+".txt";
				subtype.getBindingModel(0).printDensityToFile(distribFile);
				subtypeC++;
			}
			
			if (subtypes.size() > 0)
				bindingManager.addPotentialBindingSubtypes(cond, subtypes);
		}//End of condition
	}
	
	/**
	 * Group peaks by overlap in other conditions and sort based on responsibilities
	 * @param allpeaks
	 * @return sorted peaks
	 */
	public Map<ExperimentCondition, List<List<BindingSubComponents>>> groupPeaksByOverlap(Map<ExperimentCondition, List<BindingSubComponents>> allpeaks){
		Map<ExperimentCondition, List<List<BindingSubComponents>>> groupPeaks = new HashMap<ExperimentCondition, List<List<BindingSubComponents>>>();
		for (ExperimentCondition cond : manager.getConditions())
			groupPeaks.put(cond, new ArrayList<List<BindingSubComponents>>());		
		//Group peaks 
		if (manager.getNumConditions()>1){
			for (int i=0; i < manager.getNumConditions(); i++){
				for (int j=i+1; j< manager.getNumConditions(); j++){
					List<BindingSubComponents> peaksA = allpeaks.get(manager.getIndexedCondition(i));
					List<BindingSubComponents> peaksB = allpeaks.get(manager.getIndexedCondition(j));
					List<BindingSubComponents> overlapA = new ArrayList<BindingSubComponents>();
					List<BindingSubComponents> overlapB = new ArrayList<BindingSubComponents>();
					for (int bcA=0; bcA < peaksA.size(); bcA++){
						for (int bcB=0; bcB< peaksB.size(); bcB++){
							if (peaksA.get(bcA).getCoord().expand(config.EM_MU_UPDATE_WIN/2).overlaps(peaksB.get(bcB).getCoord().expand(config.EM_MU_UPDATE_WIN/2))){
								if (!overlapA.contains(peaksA.get(bcA)))
									overlapA.add(peaksA.get(bcA));	
								if (!overlapB.contains(peaksB.get(bcB)))
									overlapB.add(peaksB.get(bcB));
							}}}	
					if (overlapA.size() > config.getMinRefsForBMUpdate())
						groupPeaks.get(manager.getIndexedCondition(i)).add(overlapA);
					if (overlapB.size() > config.getMinRefsForBMUpdate())
						groupPeaks.get(manager.getIndexedCondition(j)).add(overlapB);			
				}
			}
			for (ExperimentCondition cond : manager.getConditions()){
				List<BindingSubComponents> peaks = allpeaks.get(cond);
				for (List<BindingSubComponents> bc : groupPeaks.get(cond)){
					System.out.println("size of group peaks "+bc.size());
					peaks.removeAll(bc);
				}
				groupPeaks.get(cond).add(peaks);
			}
		}else{
			for (ExperimentCondition cond : manager.getConditions())
				groupPeaks.get(cond).add(allpeaks.get(cond));
		}
		
		for (ExperimentCondition cond : manager.getConditions()){
			for (List<BindingSubComponents> peaks : groupPeaks.get(cond)){
				//Sort by responsibilities
				Collections.sort(peaks, new Comparator<BindingSubComponents>(){
					public int compare(BindingSubComponents o1, BindingSubComponents o2) {return o1.compareByResp(o2);}
				});
				Collections.reverse(peaks);
			}
		}
		
		return groupPeaks;	
	}
	
	public Map<ExperimentCondition, List<List<BindingSubComponents>>> groupPeaksBySubtypes(Map<ExperimentCondition, List<BindingSubComponents>> allpeaks){
		Map<ExperimentCondition, List<List<BindingSubComponents>>> groupPeaks = new HashMap<ExperimentCondition, List<List<BindingSubComponents>>>();
		for (ExperimentCondition cond : manager.getConditions()){
			groupPeaks.put(cond, new ArrayList<List<BindingSubComponents>>());	
			// first group profiles based on binding type association, check to see if there are more than two binding types first
            List<List<BindingSubComponents>> typeComps = new ArrayList<List<BindingSubComponents>>();
            // Initialize
            for (int bt=0; bt< bindingManager.getNumBindingType(cond); bt++) 
                typeComps.add(new ArrayList<BindingSubComponents>());
            // Group binding components based on subtype association
            for(BindingSubComponents comp : allpeaks.get(cond))      
                typeComps.get(comp.getMaxType()).add(comp);    
            
            List<List<BindingSubComponents>> comps2Remove = new ArrayList<List<BindingSubComponents>>(); 
            for (List<BindingSubComponents> btypeComps : typeComps)
            	if (btypeComps.size() < Math.max(config.getMinSubtypeFraction()*allpeaks.get(cond).size(), config.getMinRefsForBMUpdate()))
            		comps2Remove.add(btypeComps);
            typeComps.removeAll(comps2Remove);
            
            for (List<BindingSubComponents> peaks : typeComps){
            	Collections.sort(peaks, new Comparator<BindingSubComponents>(){
//					public int compare(BindingSubComponents o1, BindingSubComponents o2) {return o1.compareByTauProb(o2);}
            		public int compare(BindingSubComponents o1, BindingSubComponents o2) {return o1.compareByResp(o2);}
				});
				Collections.reverse(peaks);				
            }
            //
            groupPeaks.put(cond, typeComps);
		}
		return groupPeaks;		
	}
		
	/**
	 * Align motifs to get relative offsets
	 */
	public void alignMotifs(){		
		// Index 0 motif in the condition is the reference
		for (ExperimentCondition cond : manager.getConditions()){
			List<WeightMatrix> condMotifs = new ArrayList<WeightMatrix>();
			List<WeightMatrix> condFreqMatrix = new ArrayList<WeightMatrix>();
			List<Integer> subIndexes = new ArrayList<Integer>();
			List<BindingSubtype> condSubtypes = new ArrayList<BindingSubtype>();
			for (int i=0; i < bindingManager.getNumBindingType(cond); i++){
				BindingSubtype subtype = bindingManager.getBindingSubtype(cond).get(i);
				condSubtypes.add(subtype);
				if (subtype.hasMotif()){
					condMotifs.add(subtype.getMotif());
					condFreqMatrix.add(subtype.getFreqMatrix());
					subIndexes.add(i);					
				}
			}
			if (condFreqMatrix.size()>1){
				WeightMatrix refMotif = condFreqMatrix.get(0);
				BindingSubtype refSubtype = condSubtypes.get(subIndexes.get(0));
				refSubtype.setMotifOffset(0);
				refSubtype.setReverseMotif(false);
				for (int index=1; index<condFreqMatrix.size(); index++){
					BindingSubtype currSubtype = condSubtypes.get(subIndexes.get(index));
					Pair<Integer,Double> forAlignment = aligner.align(refMotif, condFreqMatrix.get(index));
					Pair<Integer,Double> revAlignment = aligner.align(refMotif, WeightMatrix.reverseComplement(condFreqMatrix.get(index)));
					int refOffset=(int)((refMotif.length()-condFreqMatrix.get(index).length())/2);
					
					if(revAlignment.cdr()>forAlignment.cdr() && refMotif.length()%2 ==0 && condFreqMatrix.get(index).length()%2 ==0)
						refOffset = (int)((refMotif.length()-condFreqMatrix.get(index).length())/2-1);
					if(revAlignment.cdr()>forAlignment.cdr()){
						currSubtype.setMotif(WeightMatrix.reverseComplement(condMotifs.get(index)), WeightMatrix.reverseComplement(condFreqMatrix.get(index)));
						currSubtype.setMotifOffset(refOffset+revAlignment.car()); // is this correct ?
						currSubtype.setReverseMotif(true);
					}else{
						currSubtype.setMotif(condMotifs.get(index), condFreqMatrix.get(index));
						currSubtype.setMotifOffset(refOffset+forAlignment.car()); // is this correct ?
						currSubtype.setReverseMotif(false);
					}
				}
				bindingManager.setBindingSubtypes(cond, condSubtypes);
			}
		}
	}
	
	public double motifAlignMaxScore(WeightMatrix wma, WeightMatrix wmb){
		Pair<Integer,Double> forAlignment = aligner.align(wma,wmb);
		Pair<Integer,Double> revAlignment = aligner.align(wma,WeightMatrix.reverseComplement(wmb));
		double maxscore = 0;
		if (forAlignment.cdr() > maxscore) { maxscore=forAlignment.cdr();}
		if (revAlignment.cdr() > maxscore) { maxscore=revAlignment.cdr();}
		return maxscore;
	}

	/**
	 * Scan the region with the motifs from each condition. Suitable for making a motif prior.  Returns score for either strand.
	 * @param reg
	 * @param boolean indicating forward or reverse strand
	 * @return array of array of motif scores. Indexed by condition. 
	 */
	public double[][][] scanStrandedRegionWithMotifs(Region reg, String regSeq, boolean forwardScores){		
		double[][][] scanScores = new double[manager.getNumConditions()][][];
				
		for(ExperimentCondition cond : manager.getConditions()){
			List<BindingSubtype> subtypes = bindingManager.getBindingSubtype(cond);
			int e = cond.getIndex();
			scanScores[e] = new double[subtypes.size()][reg.getWidth()];
			for (int i=0; i<subtypes.size(); i++)
				for(int z=0; z<reg.getWidth(); z++)
					scanScores[e][i][z]=0;			
			int i=0;	//Subtype counter
			for (BindingSubtype sub : subtypes){
				if (sub.hasMotif()){
					int motifLen = sub.getMotif().length();
					int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;					
					WeightMatrixScorer scorer = new WeightMatrixScorer(sub.getMotif());
					WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
					// profile stores either forward or reverse strand scores
					double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
					for(int z=0; z<reg.getWidth()-motifLen+1; z++){
						double currScore= scores[z];
						int zOff = z+sub.getMotifOffset()+halfMotifWidth;
						if(currScore>0)
							if (zOff>=0 && zOff<reg.getWidth())
								scanScores[e][i][zOff] = currScore;						
					}			
				}		
			i++;	
			}
		}
		return scanScores;
	}
	
	/**
	 * Scan the region with the motifs from each condition. Suitable for making a motif prior.  Returns score for either strand.
	 * @param reg
	 * @param boolean indicating forward or reverse strand
	 * @return array of array of motif scores. Indexed by motif and condition. 
	 */
	public double[][] scanStrandedRegionWithMotifsByCondition(ExperimentCondition cond, Region reg, String regSeq, boolean forwardScores){
		List<BindingSubtype> subtypes = bindingManager.getBindingSubtype(cond);
		double[][]scanScores = new double[subtypes.size()][reg.getWidth()];
		for (int i=0; i<subtypes.size(); i++)
			for(int z=0; z<reg.getWidth(); z++)
				scanScores[i][z]=0;			
		int i=0;	//Subtype counter
		for (BindingSubtype sub : subtypes){
			if (sub.hasMotif()){
				int motifLen = sub.getMotif().length();
				int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;					
				WeightMatrixScorer scorer = new WeightMatrixScorer(sub.getMotif());
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
				// profile stores either forward or reverse strand scores
				double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
				for(int z=0; z<reg.getWidth()-motifLen+1; z++){
					double currScore= scores[z];
					int zOff = z+sub.getMotifOffset()+halfMotifWidth;
					if(currScore>0)
						if (zOff>=0 && zOff<reg.getWidth())
							scanScores[i][zOff] = currScore;						
				}			
			}		
		i++;	
		}
		return scanScores;
	}
		
	
	/**
	 * Scan the region with the motifs from each condition.   
	 * @param reg
	 * @return Pair of: array of array of motif max seqs and array of array of motif scores. Indexed by condition. 
	 */
	public Pair<Double[][][],String[][][]> scanRegionWithMotifsGetSeqs(Region reg, String regSeq, boolean forwardScores){
		Double[][][] scanScores = new Double[manager.getNumConditions()][][];
		String[][][] scanSeqs = new String[manager.getNumConditions()][][];	
		for(ExperimentCondition cond : manager.getConditions()){
			List<BindingSubtype> subtypes = bindingManager.getBindingSubtype(cond);
			int e = cond.getIndex();
			scanScores[e] = new Double[subtypes.size()][reg.getWidth()];
			scanSeqs[e] = new String[subtypes.size()][reg.getWidth()];
			for (int i=0; i<subtypes.size(); i++){
				for(int z=0; z<reg.getWidth(); z++){
					scanScores[e][i][z]=0.0;	
					scanSeqs[e][i][z]="";
				}
			}
			int i=0;	//Subtype counter
			for (BindingSubtype sub : subtypes){
				if (sub.hasMotif()){
					int motifLen = sub.getMotif().length();
					int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;					
					WeightMatrixScorer scorer = new WeightMatrixScorer(sub.getMotif());
					WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
					// profile stores either forward or reverse strand scores
					double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
					for(int z=0; z<reg.getWidth()-motifLen+1; z++){
						double currScore= scores[z];
						String currSeq = regSeq.substring(z, z+sub.getMotif().length());
						if (!forwardScores)
							currSeq = SequenceUtils.reverseComplement(currSeq);
						int zOff = z+sub.getMotifOffset()+halfMotifWidth;
						if (zOff>=0 && zOff<reg.getWidth()){
							if(currScore>0)
								scanScores[e][i][zOff] = currScore;
							scanSeqs[e][i][zOff] = currSeq;
						}					
					}			
				}		
			i++;	
			}
		}
		Pair<Double[][][],String[][][]> scoresAndSeqs = new Pair<Double[][][],String[][][]>(scanScores, scanSeqs);
		return scoresAndSeqs;
	}
	
	/**
	 * Get a sequence for a given region
	 */
	public String getSeq(Region reg){
		String regSeq = seqgen.execute(reg);
		return regSeq;
	}
	
	/**
	 * Scan the region with the motif from one condition. Suitable for making a motif prior.  
	 * @param cond
	 * @param reg
	 * @return array of array of motif scores. Indexed by condition. 
	 */
	public double[][] scanRegionWithMotif(Region reg, ExperimentCondition cond){
		List<BindingSubtype> subtypes = bindingManager.getBindingSubtype(cond);
		double[][] scanScores = new double[subtypes.size()][reg.getWidth()];
		String regSeq = seqgen.execute(reg);
		for (int i=0; i<subtypes.size(); i++)
			for(int z=0; z<reg.getWidth(); z++)
				scanScores[i][z]=0;			
		int i=0;	//Subtype counter
		for (BindingSubtype sub : subtypes){
			if (sub.hasMotif()){
				WeightMatrixScorer scorer = new WeightMatrixScorer(sub.getMotif());
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
				for(int z=0; z<reg.getWidth()-sub.getMotif().length()+1; z++){
					double currScore= profiler.getMaxScore(z);
					int zOff = z+sub.getMotifOffset();
					if(currScore>0){
						if(zOff>0 && zOff<reg.getWidth())
							scanScores[i][zOff] = currScore;
					}
				}
			}		
		i++;	
		}
		return scanScores;
	}

	/**
	 * Compute the fraction of letters in the sequence that are lowercase or N
	 * @param seq
	 * @return
	 */
	protected double lowercaseFraction(String seq){
		double count = 0;
		for (char c:seq.toCharArray())
			if (Character.isLowerCase(c) || c=='N')
				count++;
		return count/(double)seq.length();
	}
	/**
	 * Randomly pick a set of Regions
	 * @param gen
	 * @param blackList
	 * @param numSamples
	 * @param sampleSize
	 * @return
	 */
	protected List<Region> randomRegionPick(List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
		Random rand = new Random(config.RANDOMSEED);
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[config.getGenome().getChromList().size()];
		String [] chromoNames = new String[config.getGenome().getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(config.getGenome());
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false, valid=true;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					valid=true;
					
					int pstart = (int)(randPos-total+1);
					int pend = (int)(pstart+sampleSize);
					if(pstart>0 && pend<chromoSize[c]-1){
						potential = new Region(config.getGenome(), chromoNames[c], pstart, pend);
						//is this region in the blacklist? 
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}

				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
	
	/**
	 * Calculate ROC scores for all motifs. 
	 * @param matrices
	 * @param posSeqs
	 * @param negSeqs
	 * @return
	 */
	protected double[] motifROCScores(List<WeightMatrix> matrices, List<String> posSeqs, String[] negSeqs){
		double[] rocScores = new double[matrices.size()];
		int m=0;
		for(WeightMatrix motif : matrices){
			List<Double> posScores = new ArrayList<Double>();
			List<Double> negScores = new ArrayList<Double>();
			if(motif!=null){
				WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
				for(String posSeq : posSeqs){
					WeightMatrixScoreProfile profiler = scorer.execute(posSeq);
					posScores.add(profiler.getMaxScore());
				}
				for(int s=0; s<negSeqs.length; s++){
					WeightMatrixScoreProfile profiler = scorer.execute(negSeqs[s]);
					negScores.add(profiler.getMaxScore());
				}
			}
			rocScores[m] = calcROCAUC(posScores, negScores);
			m++;
		}
		return rocScores;
	}
	
	/**
	 * Calculate the area under a motif-scoring ROC
	 * @param posMaxScores
	 * @param negMaxScores
	 * @param printROC
	 * @param motif
	 * @return
	 */
	private double calcROCAUC(List<Double> posMaxScores, List<Double> negMaxScores) {
		double auc = 0;
		if(posMaxScores.size()==0)
			return 0;
		if(negMaxScores.size()==0)
			return 1;
		ArrayList<LabeledDouble> data = new ArrayList<LabeledDouble>();
		for(Double d : posMaxScores)
			data.add(new LabeledDouble(d, 1));
		for(Double d : negMaxScores)
			data.add(new LabeledDouble(d, 0));
		
		Collections.sort(data);
		double pCount = (double)posMaxScores.size();
		double nCount = (double)negMaxScores.size();
		int x=0;
		double possum=0;
		double lastsn=0;
		double lastfpr=0;
		double lastdval = 10000000;
		
		for(LabeledDouble d : data){
			possum+=d.label;
			if(d.dat!=lastdval){
				double sn = possum/pCount;
				double fp = (x+1)-possum;
				double sp = (nCount-fp)/nCount;
				double fpr=1-sp;
				if(x>0){
						    //Rectangle             //Triangle
					auc += ((fpr-lastfpr)*lastsn) + ((sn-lastsn)*(fpr-lastfpr)/2);
				}
				lastfpr=fpr;
				lastsn = sn;
			}
			lastdval = d.dat;
			x++;
		}
		return auc;
	}
	
	// Return motif hit position
	protected StrandedPoint getMotifPosition(WeightMatrix motif, double motifThres, Region reg){		
		StrandedPoint p = null;
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif, seqgen);       						
		WeightMatrixScoreProfile profiler = scorer.execute(reg);
		boolean goodMotif=false;
		for(int z=0; z<reg.getWidth(); z++){
			double currScore= profiler.getMaxScore(z);
			if(currScore>=motifThres)
				goodMotif=true;	
		}
		if(goodMotif){
			int halfMotifWidth = (motif.length() % 2 == 0 && profiler.getMaxStrand() == '+') ? motif.length()/2-1 : motif.length()/2;
			int maxCoord = profiler.getMaxIndex()+reg.getStart()+halfMotifWidth;
			p= new StrandedPoint(genome, reg.getChrom(), maxCoord, profiler.getMaxStrand());
		}
		return p;
	}
	
	
	/**
	 * Simple class for ROC analysis
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	protected class LabeledDouble implements Comparable<LabeledDouble>{
		public Double dat;
		public Integer label;
		public LabeledDouble(Double d, Integer i){dat=d; label=i;}
		public int compareTo(LabeledDouble ld) {
			if(dat > ld.dat){return(-1);}
			else if(dat < ld.dat){return(1);}
			else{return 0;}
		}
	}
	
	/**
	 * Main to test local sequence alignment
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
			
			ChExMix gps = new ChExMix(gcon, econ, evconfig, config, xlconfig, sc, manager);
			
			MotifPlatform motifFinder = new MotifPlatform(gcon, config, manager, gps.getBindingManager(), gps.getPotentialFilter().getPotentialRegions());
			
			List<List<StrandedPoint>> initClustPoints = config.getInitialClustPoints();			
			
			// Test on clustered points
			List<List<List<StrandedPoint>>> clustPoints = new ArrayList<List<List<StrandedPoint>>>();
			for (ExperimentCondition cond : manager.getConditions())
				clustPoints.add(initClustPoints);
			
			List<List<List<StrandedPoint>>> adjustedPoints = motifFinder.findClusterMotifs(clustPoints,0);
			
			manager.close();
		}
	}
	
}
