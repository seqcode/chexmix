package org.seqcode.projects.chexmix.motifs;

import java.io.IOException;
import java.text.ParseException;
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

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
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
import org.seqcode.gseutils.NotFoundException;
import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.motifs.MarkovMotifThresholdFinder;
import org.seqcode.projects.chexmix.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.framework.ProfileCluster;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.mixturemodel.BindingSubComponents;


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
	protected ProfileCluster cluster;
		
	/**
	 * Constructor for motif platform
	 * @param c
	 * @param man
	 * @param regionsOfInterest: this list contains all possible regions that motif-finding/scanning may be run on. 
	 */
	public MotifPlatform(GenomeConfig g, ChExMixConfig c, ExperimentManager man, BindingManager bman, List<Region> regionsOfInterest, ProfileCluster cluster ){
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
		meme = new MEMERunner(config, man);
		
		//Load the background model or make background model
        try{       
        	if(config.getBackModel() == null)
        		backMod = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(Genome.findGenome(genome.getDBID()))); // this doesn't seem working for sacCer3
        	else
				backMod = BackgroundModelIO.parseMarkovBackgroundModel(config.getBackModel(), Genome.findGenome(genome.getDBID()));
        } catch (IOException | ParseException | NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        //Pre-load the random sequences
		RandomSequenceGenerator randGen = new RandomSequenceGenerator(backMod);
		for(int i=0; i<config.MARKOV_NUM_TEST; i++){
			MarkovRandSeq.add(randGen.execute(config.MOTIF_FINDING_SEQWINDOW));
		}
		
		aligner = new SimpleMotifAligner(config.getMinMotifLength());
		this.cluster = cluster;
	}
	
	/**
	 * Extract sequences around top BindingComponents and call the MEME runner 
	 * Remove sequences with motifs and recursively do motif discovery
	 * @param cond
	 * @param activeComponents
	 * @param trainingRound
	 */
	public void recursivelyFindMotifs(HashMap<Region, List<List<BindingSubComponents>>> activeComponents, int trainingRound){
		Map<ExperimentCondition, List<BindingSubComponents>> allpeaks = new HashMap<ExperimentCondition, List<BindingSubComponents>>();
		for (ExperimentCondition cond : manager.getConditions())
			allpeaks.put(cond, new ArrayList<BindingSubComponents>());
		//Choose which components to include
		for(Region r : activeComponents.keySet()){
			for (ExperimentCondition cond : manager.getConditions()){
				float edge = Math.max(bindingManager.getMaxInfluenceRange(cond)/2, config.MOTIF_FINDING_SEQWINDOW/2);
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
		if (trainingRound ==1 || trainingRound ==2)
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
				int counter = 0;   	
				if (sortedRegions.size() < config.getMinRefsForBMUpdate()){
					findMotif = false;
					System.err.println("cannot do motif finding due to too few regions ("+sortedRegions.size()+"<"+config.getMinRefsForBMUpdate()+").");
				}
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
							Set<StrandedPoint> refs = getRegionWithMotifs(currMotif,config.MARKOV_BACK_MODEL_THRES, sortedRegions).cdr(); // motif references
							List<Region> found = getRegionWithMotifs(currMotif,config.MARKOV_BACK_SEQ_RM_THRES, sortedRegions).car();	//region with motif to be removed
							
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
				
				// Keep list of peaks that can be used for BM updates later
				if (!groupFoundMotif){
					if (peaks.size() < config.getMinComponentsForBMUpdate())
						System.err.println("The "+cond.getName()+" read distributions cannot be updated due to too few binding components ("+peaks.size()+"<"+config.getMinComponentsForBMUpdate()+").");
					else
						groupComps.add(peaks);	
				}
				setCounter++; // Next set of sequences
			} // End of iteration grouped sequence from a condition
			
			// Calculate scores for all motifs
			int motifSize = selecFreqMatrix.size();
			if (motifSize > 1){ // access motif similarities if more than one motif passes ROC test	
				double[][] freqMatrixScores = new double[motifSize][motifSize];
				for (int i=0; i < motifSize; i++)
					for (int j=i+1; j < motifSize; j++)
						freqMatrixScores[i][j]=0;
				for (int i=0; i < motifSize; i++){
					for (int j=i+1; j < motifSize; j++){
						Pair<Integer,Double> forAlignment = aligner.align(selecFreqMatrix.get(i),selecFreqMatrix.get(j));
						Pair<Integer,Double> revAlignment = aligner.align(selecFreqMatrix.get(i),WeightMatrix.reverseComplement(selecFreqMatrix.get(j)));
						double maxscore = 0;
						if (forAlignment.cdr() > maxscore) { maxscore=forAlignment.cdr();}
						if (revAlignment.cdr() > maxscore) { maxscore=revAlignment.cdr();}
						System.out.println("max score is "+maxscore/config.getMinMotifLength());
						freqMatrixScores[i][j]=maxscore/config.getMinMotifLength();	
					}
				}
				// Identify motifs with similar frequency matrix
				List<Integer> motif2remove = new ArrayList<Integer>();		
				boolean assessSimilarity=true;
				do{			
					double maxscore=0;
					int indexA=0; int indexB=0;
					for (int i=0; i < motifSize; i++){
						for (int j=i+1; j < motifSize; j++){
							if (!motif2remove.contains(i) && !motif2remove.contains(j)){
								if (freqMatrixScores[i][j] > Math.max(maxscore,config.MOTIF_PCC_THRES)){
									maxscore =freqMatrixScores[i][j];
									indexA=i; indexB=j;								
							}}}}	
					if (maxscore > 0){
						/**
						double sumScoreA=0; double sumScoreB=0;
						for (int i=0; i < motifSize; i++){
							for (int j=i+1; j < motifSize; j++){
								if (!motif2remove.contains(i) && !motif2remove.contains(j)){
									if (i==indexA || j==indexA){ sumScoreA+=freqMatrixScores[i][j];}
									if (i==indexB || j==indexB){ sumScoreB+=freqMatrixScores[i][j];}
								}}}		
						if (sumScoreA > sumScoreB)//remove index A
							motif2remove.add(indexA);
						else //remove index B
							motif2remove.add(indexB);
						**/
						double KLLogScore = DistributionSimilarityScore(cond, motifRefs.get(indexA),motifRefs.get(indexB));
						if (KLLogScore < config.KL_DIVERGENCE_BM_THRES){
							// models are similar
							// Instead of picking the motif that maximizes sum of dissimilarity score, pick the motif that is associated with more reference
							if (motifRefs.get(indexA).size() < motifRefs.get(indexB).size())
								motif2remove.add(indexA);
							else
								motif2remove.add(indexB);
						}else{
							freqMatrixScores[indexA][indexB]=0;
						}
					}else{
						assessSimilarity=false;
					}
				}while(assessSimilarity); // End of while loop			
				
				System.out.println("motifs to be removed "+motif2remove.toString());
				System.out.println("number of motifs before removing "+selecMotifs.size());
				removeWeightMatrix(selecMotifs, motif2remove);
				removeWeightMatrix(selecFreqMatrix, motif2remove);
				removeMotifReferences(motifRefs,motif2remove);
				System.out.println("number of motifs after removing "+selecMotifs.size());				
			}			

			if (selecMotifs.size()==0){
				System.err.println("\tNo motif passes minimum ROC score threshold.");
				bindingManager.setMotif(cond, null);
				bindingManager.setFreqMatrix(cond, null);
				bindingManager.setMotifReferece(cond, null);
				bindingManager.setFreqMatrix(cond, null);
			}else{
				//Set the condition's motif if the ROC is above threshold
				ArrayList<Integer> zeros = new ArrayList<Integer>();
				for (WeightMatrix m : selecMotifs){
					System.err.println(WeightMatrix.getConsensus(m) + " chosen as motif above threshould.");
					zeros.add(0);
				}
				bindingManager.setMotif(cond, selecMotifs);
				bindingManager.setMotifOffset(cond, zeros);
				bindingManager.setMotifReferece(cond,motifRefs);
				bindingManager.setFreqMatrix(cond, selecFreqMatrix);
			}
			bindingManager.setComponentsForBMUpdates(cond,groupComps);
			
		} //End of Motif finding for condition
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
					if (overlapA.size() > config.getMinRefsForBMUpdate())
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
            	// Sort by tau probabilities
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
	
	public double DistributionSimilarityScore(ExperimentCondition cond, Set<StrandedPoint> refA, Set<StrandedPoint> refB){
		List<StrandedPoint> compositePointsA = new ArrayList<StrandedPoint>();
		List<StrandedPoint> compositePointsB = new ArrayList<StrandedPoint>();
        compositePointsA.addAll(refA);
        compositePointsB.addAll(refB);
		CompositeTagDistribution signalCompositeA = new CompositeTagDistribution(compositePointsA, cond, config.MAX_BINDINGMODEL_WIDTH,true);
		CompositeTagDistribution signalCompositeB = new CompositeTagDistribution(compositePointsB, cond, config.MAX_BINDINGMODEL_WIDTH,true);
		// Calculate KL divergence score in sliding window and take the lowest
		double[] watsonA = signalCompositeA.getCompositeWatson();
		double[] crickA = signalCompositeA.getCompositeCrick();
		double[] watsonB = signalCompositeB.getCompositeWatson();
		double[] crickB = signalCompositeB.getCompositeCrick();
		double[] rwatsonB = new double[signalCompositeB.getWinSize()];
		double[] rcrickB = new double[signalCompositeB.getWinSize()];
		for (int i=0; i < signalCompositeB.getWinSize(); i++){
			rwatsonB[i] = crickB[signalCompositeB.getWinSize()-i-1];
			rcrickB[i] = watsonB[signalCompositeB.getWinSize()-i-1];	
		}
		
		double minLogKL = Double.MAX_VALUE;
		for (int offset=-config.PCC_SLIDING_WINDOW/2; offset<=config.PCC_SLIDING_WINDOW/2; offset++){
            double currLogKL=0;
            //copy current window
            double[] currW = new double[config.MAX_BINDINGMODEL_WIDTH];
            double[] currC = new double[config.MAX_BINDINGMODEL_WIDTH];
            for (int w=0; w< config.MAX_BINDINGMODEL_WIDTH; w++){
            	currW[w]=0; currC[w]=0;
            }
            for (int w=0; w< config.MAX_BINDINGMODEL_WIDTH; w++){
            	if ((offset+w) >=0 && (offset+w) < config.MAX_BINDINGMODEL_WIDTH){
            		currW[w]=watsonB[offset+w];
            		currC[w]=crickB[offset+w];
            	}
            }
            //Calc KL                                       
            currLogKL += StatUtil.log_KL_Divergence(watsonA, currW) + StatUtil.log_KL_Divergence(currW, watsonA);
            currLogKL += StatUtil.log_KL_Divergence(crickA, currC) + StatUtil.log_KL_Divergence(currC, crickA); 
   
            if (currLogKL < minLogKL) { minLogKL=currLogKL;}
                    
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
            currLogKL += StatUtil.log_KL_Divergence(watsonA, rcurrW) + StatUtil.log_KL_Divergence(rcurrW, watsonA);
            currLogKL += StatUtil.log_KL_Divergence(crickA, rcurrC) + StatUtil.log_KL_Divergence(rcurrC, crickA);  
            
            if (currLogKL < minLogKL) { minLogKL=currLogKL;}
        }
		return minLogKL;
	}
	
	
	
	/**
	 * Align motifs to get relative offsets
	 */
	public void alignMotifs(){		
		// Index 0 motif in the condition is the reference
		WeightMatrix refMotif=null;
		for (ExperimentCondition cond : manager.getConditions()){
			List<WeightMatrix> fm = bindingManager.getFreqMatrices(cond);
			if (fm!=null){
				refMotif = fm.get(0);
				if (refMotif!=null){
					List<WeightMatrix> alignMotifs = new ArrayList<WeightMatrix>();
					List<WeightMatrix> alignFreqMatrix = new ArrayList<WeightMatrix>();
					List<Integer> offsets = new ArrayList<Integer>();
					List<Boolean> reverseStrands = new ArrayList<Boolean>();
					//Reference offset is the center of the motif
					offsets.add(0);
					reverseStrands.add(false);
					alignMotifs.add(bindingManager.getMotifs(cond).get(0));
					alignFreqMatrix.add(fm.get(0));
					for (int index=1; index<fm.size(); index++){
						Pair<Integer,Double> forAlignment = aligner.align(refMotif, fm.get(index));
						Pair<Integer,Double> revAlignment = aligner.align(refMotif, WeightMatrix.reverseComplement(fm.get(index)));
						int refOffset=(int)((refMotif.length()-fm.get(index).length())/2);
						if(revAlignment.cdr()>forAlignment.cdr() && refMotif.length()-fm.get(index).length()%2 !=0)
							refOffset = (int)((refMotif.length()-fm.get(index).length())/2-1);
						System.out.println("fscore "+forAlignment.cdr()+" alignment offset "+forAlignment.car()+" total offset "+(refOffset+forAlignment.car()));
						System.out.println("rscore "+revAlignment.cdr()+" alignment offset "+revAlignment.car()+" total offset "+(refOffset+revAlignment.car()));
						if(revAlignment.cdr()>forAlignment.cdr()){
							alignFreqMatrix.add(WeightMatrix.reverseComplement(fm.get(index)));
							alignMotifs.add(WeightMatrix.reverseComplement(bindingManager.getMotifs(cond).get(index)));
							offsets.add(refOffset+revAlignment.car());	// is this correct ?
							reverseStrands.add(true);
						}else{
							alignFreqMatrix.add(fm.get(index));
							alignMotifs.add(bindingManager.getMotifs(cond).get(index));
							offsets.add(refOffset+forAlignment.car());
							reverseStrands.add(false);
						}
					}
					bindingManager.setMotif(cond, alignMotifs);
					bindingManager.setFreqMatrix(cond, alignFreqMatrix);
					bindingManager.setMotifOffset(cond, offsets);
					bindingManager.setReverseStrand(cond, reverseStrands);
					System.out.println("offset "+offsets.toString()+" reverse "+reverseStrands.toString());
					
				}
			}
		}
	}

	
	/**
	 * Scan the region with the motifs from each condition. Suitable for making a motif prior.  
	 * @param reg
	 * @return array of array of motif scores. Indexed by condition. 
	 */
	public double[][][] scanRegionWithMotifs(Region reg, String regSeq){			
		double[][][] scanScores = new double[manager.getNumConditions()][][];		
		for(ExperimentCondition cond : manager.getConditions()){
			List<WeightMatrix> motifs = bindingManager.getMotifs(cond);
			int e = cond.getIndex();
			scanScores[e] = new double[motifs.size()][reg.getWidth()];
			for(int z=0; z<reg.getWidth(); z++)
				for (int i=0; i<motifs.size(); i++)
					scanScores[e][i][z]=0;
			
			if(motifs!=null){
				List<Integer> motifOffsets = bindingManager.getMotifOffsets(cond);
				for (int i=0; i < motifs.size(); i++){
					WeightMatrixScorer scorer = new WeightMatrixScorer(motifs.get(i));
					WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
					Integer halfMotifWidth = motifs.get(i).length()/2;
					for(int z=0; z<reg.getWidth()-motifs.get(i).length()+1; z++){
						double currScore= profiler.getMaxScore(z);
						int zOff = z+motifOffsets.get(i)+halfMotifWidth;
						if(currScore>0){
							if(zOff > 0 && zOff < reg.getWidth())
								scanScores[e][i][zOff] = currScore;
						}
					}
				}
			}
		}
		return scanScores;
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
			int e = cond.getIndex();
			List<WeightMatrix> motifs = bindingManager.getMotifs(cond);
			
			if(motifs!=null){
				scanScores[e] = new double [motifs.size()][reg.getWidth()];
				for(int z=0; z<reg.getWidth(); z++)
					for (int i=0; i <motifs.size(); i++)
						scanScores[e][i][z]=0;
				
				List<Integer> motifOffsets = bindingManager.getMotifOffsets(cond);
				for (int i=0; i < motifs.size(); i++){	
					int motifLen = bindingManager.getMotifs(cond).get(i).length();
					int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;
					WeightMatrixScorer scorer = new WeightMatrixScorer(motifs.get(i));
					WeightMatrixScoreProfile profiler = scorer.execute(regSeq);				
					// profile stores either forward or reverse strand scores
					double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
					for(int z=0; z<reg.getWidth()-motifs.get(i).length()+1; z++){
						double currScore= scores[z];
						int zOff = z+motifOffsets.get(i)+halfMotifWidth;
						if(currScore>0){
							if (zOff>=0 && zOff<reg.getWidth())
								scanScores[e][i][zOff] = currScore;
						}
					}
				}
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
		double[][] scanScores = null;
		List<WeightMatrix> motifs = bindingManager.getMotifs(cond);
		if(motifs!=null){
			scanScores = new double [motifs.size()][reg.getWidth()];
			List<Integer> motifOffsets = bindingManager.getMotifOffsets(cond);
			for (int i=0; i < motifs.size(); i++){
				for(int z=0; z<reg.getWidth(); z++)
					scanScores[i][z]=0;
				
				int motifLen = motifs.get(i).length();
				int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;
				WeightMatrixScorer scorer = new WeightMatrixScorer(motifs.get(i));
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);				
				// profile stores either forward or reverse strand scores
				double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
				for(int z=0; z<reg.getWidth()-motifs.get(i).length()+1; z++){
					double currScore= scores[z];
					int zOff = z+motifOffsets.get(i)+halfMotifWidth;
					if(currScore>0){
						if (zOff>=0 && zOff<reg.getWidth())
							scanScores[i][zOff] = currScore;
					}
				}
			}
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
			int e = cond.getIndex();
			List<WeightMatrix> motifs = bindingManager.getMotifs(cond);
			if(motifs!=null){			
				scanScores[e] = new Double[motifs.size()][reg.getWidth()];
				scanSeqs[e] = new String[motifs.size()][reg.getWidth()];
				for(int z=0; z<reg.getWidth(); z++){
					for (int i=0; i< motifs.size(); i++){
						scanScores[e][i][z]=0.0;
						scanSeqs[e][i][z]="";
					}
				}		
				List<Integer> motifOffsets = bindingManager.getMotifOffsets(cond);
				for (int i=0; i < motifs.size(); i++){	
					int motifLen = bindingManager.getMotifs(cond).get(i).length();
					int halfMotifWidth = (motifLen % 2 == 0 && forwardScores) ? motifLen/2-1 : motifLen/2;
					WeightMatrixScorer scorer = new WeightMatrixScorer(motifs.get(i));
					WeightMatrixScoreProfile profiler = scorer.execute(regSeq);				
					// profile stores either forward or reverse strand scores
					double[] scores = forwardScores ? profiler.getForwardScores() : profiler.getReverseScores();
					for(int z=0; z<reg.getWidth()-motifs.get(i).length()+1; z++){
						double currScore= scores[z];
						String currSeq = regSeq.substring(z, z+motifs.get(i).length());
						if (!forwardScores)
							currSeq = SequenceUtils.reverseComplement(currSeq);
						int zOff = z+motifOffsets.get(i)+halfMotifWidth;
						if (zOff>=0 && zOff<reg.getWidth()){
							if(currScore>0)
								scanScores[e][i][zOff] = currScore;
							scanSeqs[e][i][zOff] = currSeq;
						}
					}
				}
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
		List<WeightMatrix> motifs = bindingManager.getMotifs(cond);
		double[][] scanScores = new double[motifs.size()][reg.getWidth()];
		String regSeq = seqgen.execute(reg);
		for (int i=0; i < motifs.size(); i++)
			for(int z=0; z<reg.getWidth(); z++)
				scanScores[i][z]=0;
		
		if(motifs!=null){
			List<Integer> motifOffsets = bindingManager.getMotifOffsets(cond);			
			for (int i=0 ; i < motifs.size(); i++){			
				WeightMatrixScorer scorer = new WeightMatrixScorer(motifs.get(i));
				WeightMatrixScoreProfile profiler = scorer.execute(regSeq);
				for(int z=0; z<reg.getWidth()-motifs.get(i).length()+1; z++){
					double currScore= profiler.getMaxScore(z);
					int zOff = z+motifOffsets.get(i);
					if(currScore>0){
						if(zOff>0 && zOff<reg.getWidth())
							scanScores[i][zOff] = currScore;
					}
				}
			}
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
		Random rand = new Random();
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
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]-1){
						int pstart = (int)(randPos-total);
						int pend = (int)(randPos+sampleSize-total);
						potential = new Region(config.getGenome(), chromoNames[c], pstart, pend);
						
						//is this region in the blacklist? 
						boolean valid=true;
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
	
	protected void removeWeightMatrix(List<WeightMatrix> wm, List<Integer> indexes){
		// Iterate and remove similar motifs
		int counter=0;
		for (Iterator<WeightMatrix> iter = wm.listIterator(); iter.hasNext(); ) {
			WeightMatrix motif = iter.next();
		    if (indexes.contains(counter)) { iter.remove();}
		    counter++;
		}
	}
	
	protected void removeMotifReferences(List<Set<StrandedPoint>> mrefs, List<Integer> indexes){
		// Iterate and remove similar motifs
		int counter=0;
		for (Iterator<Set<StrandedPoint>> iter = mrefs.listIterator(); iter.hasNext(); ) {
			Set<StrandedPoint> refs = iter.next();
		    if (indexes.contains(counter)) { iter.remove();}
		    counter++;
		}
	}
	
	protected Pair<List<Region>, Set<StrandedPoint>> getRegionWithMotifs(WeightMatrix motif, double markov_thres, List<Region> regs){
		List<Region> regWithMotif=new ArrayList<Region>();
		Set<StrandedPoint> motifPos = new HashSet<StrandedPoint>();
		MarkovMotifThresholdFinder finder = new MarkovMotifThresholdFinder(motif, backMod, config.MARKOV_NUM_TEST);
		finder.setRandomSeq(MarkovRandSeq);
		double motifThres = finder.execute(markov_thres);       				        				
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif, seqgen);       						
		for (Region reg : regs){
			WeightMatrixScoreProfile profiler = scorer.execute(reg);
			boolean goodMotif=false;
			for(int z=0; z<reg.getWidth(); z++){
				double currScore= profiler.getMaxScore(z);
				if(currScore>=motifThres)
					goodMotif=true;	
			}
			if(goodMotif){
				regWithMotif.add(reg);
				int halfMotifWidth = (motif.length() % 2 == 0 && profiler.getMaxStrand() == '+') ? motif.length()/2-1 : motif.length()/2;
				int maxCoord = profiler.getMaxIndex()+reg.getStart()+halfMotifWidth;
				motifPos.add(new StrandedPoint(genome, reg.getChrom(), maxCoord, profiler.getMaxStrand()));
			}
		}
		return new Pair<List<Region>, Set<StrandedPoint>>(regWithMotif,motifPos);
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
}
