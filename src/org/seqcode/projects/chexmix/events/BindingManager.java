package org.seqcode.projects.chexmix.events;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.chexmix.composite.ProteinDNAInteractionModel;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.mixturemodel.BindingSubComponents;


/**
 * BindingManager stores lists of binding events and motifs associated with experiment conditions, 
 * and binding models associated with replicates. 
 * These data structures used to be under the relevant experiment components, but we moved them here to 
 * allow experiment loading to be independent.  
 * 
 * @author mahony
 *
 */
public class BindingManager {

	protected EventsConfig config;
	protected ExperimentManager manager;
	protected List<BindingEvent> events;
	protected Map<ExperimentCondition, List <BindingEvent>> conditionEvents;
	protected Map<ExperimentCondition, List<WeightMatrix>> motifs;
	protected Map<ExperimentCondition, List<Integer>> motifIndexes;
	protected Map<ExperimentCondition, List<WeightMatrix>> freqMatrices;
	protected Map<ExperimentCondition, List<Integer>> motifOffsets;
	protected Map<ExperimentCondition, List<Set<StrandedPoint>>> motifReferencePoints;
	protected Map<ExperimentCondition, List<Boolean>> motifReverseStrands;
	protected Map<ExperimentCondition, Integer> numBindingType;
	protected Map<ControlledExperiment, List<TagProbabilityDensity>> models;
	protected Map<ControlledExperiment, BindingModel> unstrandedModel;
	protected Map<ExperimentCondition, Integer> maxInfluenceRange;
	protected Map<ExperimentCondition, Double> alpha;
	protected Map<ControlledExperiment, List<ProteinDNAInteractionModel>> ProteinDNAInteractionModels;
	protected Map<ExperimentCondition, List<List<BindingSubComponents>>> ComponentsForBMUpdates;
	protected Map<ExperimentCondition, List<StrandedPoint>> alignedEventPoints;
	
	public BindingManager(EventsConfig con, ExperimentManager exptman){
		config = con;
		manager = exptman;
		events  = new ArrayList<BindingEvent>();
		conditionEvents = new HashMap<ExperimentCondition, List <BindingEvent>>();
		motifs = new HashMap<ExperimentCondition, List<WeightMatrix>>();
		motifIndexes = new HashMap<ExperimentCondition, List<Integer>>();
		freqMatrices = new HashMap<ExperimentCondition, List<WeightMatrix>>();
		motifOffsets = new HashMap<ExperimentCondition, List<Integer>>();
		motifReverseStrands = new HashMap<ExperimentCondition, List<Boolean>>();
		maxInfluenceRange = new HashMap<ExperimentCondition, Integer>();
		alpha = new HashMap<ExperimentCondition, Double>();
		motifReferencePoints = new HashMap<ExperimentCondition,List<Set<StrandedPoint>>>();
		numBindingType = new HashMap<ExperimentCondition, Integer>();
		models = new HashMap<ControlledExperiment, List<TagProbabilityDensity>>();
		unstrandedModel = new HashMap<ControlledExperiment, BindingModel>();
		ProteinDNAInteractionModels = new HashMap<ControlledExperiment, List<ProteinDNAInteractionModel>>();
		ComponentsForBMUpdates = new HashMap<ExperimentCondition, List<List<BindingSubComponents>>>();
		alignedEventPoints = new HashMap<ExperimentCondition, List<StrandedPoint>>();
		for(ExperimentCondition cond : manager.getConditions()){
			conditionEvents.put(cond, new ArrayList<BindingEvent>());
			motifOffsets.put(cond,new ArrayList<Integer>());
			motifIndexes.put(cond, new ArrayList<Integer>());
			numBindingType.put(cond, 1);
			maxInfluenceRange.put(cond,0);
			alpha.put(cond, 0.0);
			motifReferencePoints.put(cond, new ArrayList<Set<StrandedPoint>>());	
			ComponentsForBMUpdates.put(cond, new ArrayList<List<BindingSubComponents>>());
			alignedEventPoints.put(cond, new ArrayList<StrandedPoint>());
		}
	}
	
	public List<BindingEvent> getBindingEvents(){return events;}
	public List<BindingEvent> getConditionBindingEvents(ExperimentCondition ec){return conditionEvents.get(ec);}
	public List<WeightMatrix> getMotifs(ExperimentCondition ec){return motifs.get(ec);}
	public List<Integer> getMotifIndexes(ExperimentCondition ec){return motifIndexes.get(ec);}
	public List<WeightMatrix> getFreqMatrices(ExperimentCondition ec){return freqMatrices.get(ec);}
	public List<Integer> getMotifOffsets(ExperimentCondition ec){return motifOffsets.get(ec);}
	public List<Set<StrandedPoint>> getMotifReferences(ExperimentCondition ec){return motifReferencePoints.get(ec);}
	public List<Boolean> getReverseMotifs(ExperimentCondition ec){return motifReverseStrands.get(ec);}
	public Integer getNumBindingType(ExperimentCondition ec){return numBindingType.get(ec);}
	public List<TagProbabilityDensity> getBindingModel(ControlledExperiment ce){return models.get(ce);}
	public BindingModel getUnstrandedBindingModel(ControlledExperiment ce){return unstrandedModel.get(ce);}
	public Integer getMaxInfluenceRange(ExperimentCondition ec){return maxInfluenceRange.get(ec);}
	public Double getAlpha(ExperimentCondition ec){return alpha.get(ec);}
	public List<ProteinDNAInteractionModel> getProteinDNAInteractionModel(ControlledExperiment ce){return ProteinDNAInteractionModels.get(ce);}
	public List<List<BindingSubComponents>> getComponentsForBMUpdates(ExperimentCondition ec){return ComponentsForBMUpdates.get(ec);}
	public List<StrandedPoint> getAlignedEventPoints(ExperimentCondition ec){return alignedEventPoints.get(ec);}

	public void setBindingEvents(List<BindingEvent> e){events =e;}
	public void setConditionBindingEvents(ExperimentCondition ec, List<BindingEvent> e){conditionEvents.put(ec, e);}
	public void setMotif(ExperimentCondition ec, List<WeightMatrix> m){motifs.put(ec, m);}
	public void setFreqMatrix(ExperimentCondition ec, List<WeightMatrix> m){freqMatrices.put(ec, m);}
	public void setMotifOffset(ExperimentCondition ec, List<Integer> offsets){motifOffsets.put(ec, offsets);}
	public void setReverseStrand(ExperimentCondition ec, List<Boolean> reverse){motifReverseStrands.put(ec, reverse);}
	public void setMotifReferece(ExperimentCondition ec, List<Set<StrandedPoint>> motifRefs){motifReferencePoints.put(ec,motifRefs);}
	public void setMotifIndexes(ExperimentCondition ec, List<Integer> index){motifIndexes.put(ec, index);}
	public void setAlpha(ExperimentCondition ec, Double a){alpha.put(ec,a);}
	public void setBindingModel(ControlledExperiment ce, List<TagProbabilityDensity> mod){models.put(ce, mod); numBindingType.put(ce.getCondition(), mod.size());}
	public void setUnstrandedBindingModel(ControlledExperiment ce, BindingModel mod){unstrandedModel.put(ce, mod);}
	public void setProteinDNAInteractionModel(ControlledExperiment ce, List<ProteinDNAInteractionModel> model) {ProteinDNAInteractionModels.put(ce, model);}
	public void setComponentsForBMUpdates(ExperimentCondition ec, List<List<BindingSubComponents>> comps){ ComponentsForBMUpdates.put(ec, comps);}
	public void setAlignedEventPoints(ExperimentCondition ec, List<StrandedPoint> points){alignedEventPoints.put(ec, points);}

	public void updateMaxInfluenceRange(ExperimentCondition ec, boolean firstround){
		int max=0; 
		for(ControlledExperiment rep : ec.getReplicates()){
			if (firstround && config.getDefaultBindingModel()==null){
				if (getUnstrandedBindingModel(rep).getInfluenceRange()>max)
					max=getUnstrandedBindingModel(rep).getInfluenceRange();
			}else{
				int min=Integer.MAX_VALUE;
				for (TagProbabilityDensity density : getBindingModel(rep)){
					if (density.getInfluenceRange()< min)
						min=density.getInfluenceRange(); 
				}
				max=min;
			}
		}maxInfluenceRange.put(ec, max);
	}
	
	public void updateNumBindingTypes(){
		int[] numBindingTypes = new int[manager.getNumConditions()];
		for (ExperimentCondition cond : manager.getConditions())
			numBindingTypes[cond.getIndex()]=getNumBindingType(cond);
		BindingEvent.setNumBindingTypes(numBindingTypes);	
	}
		
	/**
	 * Print the events
	 * @param outRoot
	 */
	public void printConditionEvents(ExperimentCondition ec, String outRoot){
		try{
			String outName = outRoot+"."+ec.getName()+".events";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.fullHeadString()+"\n");
			for(BindingEvent e : getConditionBindingEvents(ec)){
				fw.write(e.toString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Print the replicate counts at events
	 * @param outRoot
	 */
	public void printReplicateCounts(ExperimentCondition ec, String outRoot){
		try{
			String outName = outRoot+"."+ec.getName()+".repcounts";
			FileWriter fw = new FileWriter(outName);
			
			fw.write(BindingEvent.repCountHeadString()+"\n");
			for(BindingEvent e : getConditionBindingEvents(ec)){
				fw.write(e.getRepCountString()+"\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * For each controlled experiment, simply calculate the proportion of reads in the provided 
	 * list of binding events to everything else. 
	 * @param regs
	 */
	public void estimateSignalVsNoiseFractions(List<BindingEvent> signalEvents){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				double repSigCount =0, repNoiseCount=0;
				for(BindingEvent event : signalEvents){
					if(event.isFoundInCondition(c))
						repSigCount += event.getRepSigHits(r);
				}
				repNoiseCount = r.getSignal().getHitCount() - repSigCount;
				r.setSignalVsNoiseFraction(repSigCount/r.getSignal().getHitCount());
				System.err.println(r.getName()+"\t"+r.getIndex()+"\tsignal-noise ratio:\t"+String.format("%.4f",r.getSignalVsNoiseFraction()));
			}
		}
	}
	/**
	 * Count the binding events present in a given condition
	 * @param cond
	 * @return
	 */
	public int countEventsInCondition(ExperimentCondition cond, double qMinThres){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=qMinThres)
				count++;
		}
		return count;
	}
	/**
	 * Count the differential binding events present in a given pair of conditions
	 * @param cond
	 * @return
	 */
	public int countDiffEventsBetweenConditions(ExperimentCondition cond, ExperimentCondition othercond, double qMinThres, double diffPMinThres){
		int count=0;
		for(BindingEvent e : events){
			if(e.isFoundInCondition(cond) && e.getCondSigVCtrlP(cond) <=qMinThres)
    			if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0)
    				count++;
		}
		return count;
	}
    /**
     * Print all binding events to files
     */
    public void writeBindingEventFiles(String filePrefix, double qMinThres, boolean runDiffTests, double diffPMinThres){
    	if(events.size()>0){
    		
	    	try {
	    		//Full output table (all non-zero components)
	    		String filename = filePrefix+".all.events.table";
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.fullHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.toString()+"\n");
				fout.close();
	    		
	    		//Per-condition event files
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			//Sort on the current condition
	    			BindingEvent.setSortingCond(cond);
	    			Collections.sort(events, new Comparator<BindingEvent>(){
	    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	    	        });
	    			//Print events
	    			String condName = cond.getName(); 
	    			condName = condName.replaceAll("/", "-");
	    			filename = filePrefix+"_"+condName+".events";
					fout = new FileWriter(filename);
					fout.write(BindingEvent.conditionHeadString(cond)+"\n");
			    	for(BindingEvent e : events){
			    		double Q = e.getCondSigVCtrlP(cond);
			    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
			    		if(e.isFoundInCondition(cond) && Q <=qMinThres)
			    			fout.write(e.getConditionString(cond)+"\n");
			    	}
					fout.close();
					
					//Sort on the current condition
					BindingEvent.setSortingCond(cond);
					Collections.sort(events, new Comparator<BindingEvent>(){
	    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCounts(o2);}
	    	        });
					//Print ML subtype events
					filename = filePrefix+"_"+condName+".subtype.events";
					fout = new FileWriter(filename);
					for(BindingEvent e : events){
						double Q = e.getCondSigVCtrlP(cond);
						//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
						if(e.isFoundInCondition(cond) && Q <=qMinThres)
							fout.write(e.getSubtypeString(cond)+"\n");
					}	
					fout.close();
					
					if (getMotifs(cond)!=null){
						//Print aligned points
						List<List<StrandedPoint>> subtypePoints = new ArrayList<List<StrandedPoint>>();
						List<List<StrandedPoint>> alignedPoints = new ArrayList<List<StrandedPoint>>();
						List<StrandedPoint> allPoints = new ArrayList<StrandedPoint>();
						for (int bt=0; bt < getNumBindingType(cond); bt++){
							subtypePoints.add(new ArrayList<StrandedPoint>());
							alignedPoints.add(new ArrayList<StrandedPoint>());
						}
						for(BindingEvent e : events){
							double Q = e.getCondSigVCtrlP(cond);
							//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
							if(e.isFoundInCondition(cond) && Q <=qMinThres){
								Pair<Integer, StrandedPoint>p=e.getMaxSubtypePoint(cond);
								subtypePoints.get(p.car()).add(p.cdr());
							}
						}	
						
						// Align points
						for (int bt=0;bt< getNumBindingType(cond); bt++){
							int motifIndex = getMotifIndexes(cond).get(bt);
							List<StrandedPoint> points = subtypePoints.get(bt);
							if (motifIndex != -1){
		    					int offset = getMotifOffsets(cond).get(motifIndex);
		    					boolean reverse = getReverseMotifs(cond).get(motifIndex);
		    					if (reverse){
		    						for (StrandedPoint p : points){
		    							int location = p.getStrand()=='+' ? p.getLocation()-offset : p.getLocation()+offset;
		    							alignedPoints.get(bt).add(new StrandedPoint(p.getGenome(),p.getChrom(),location,p.getStrand() =='+' ? '-' : '+'));
		    						}	    						
		    					}else{
		    						for (StrandedPoint p : points){
		    							int location = p.getStrand()=='+' ? p.getLocation()+offset : p.getLocation()-offset;
		    							alignedPoints.get(bt).add(new StrandedPoint(p.getGenome(),p.getChrom(),location,p.getStrand()));
		    						}	    						
		    					}	    					
							}else{
								for (StrandedPoint p : points)
									alignedPoints.get(bt).add(p);
							}
						}
						filename = filePrefix+"_"+condName+".subtype.aligned.events";
						fout = new FileWriter(filename);
						for (List<StrandedPoint> points : alignedPoints){
							allPoints.addAll(points);
							for (StrandedPoint p : points)
								fout.write(p.toString()+'\n');
						}
						fout.close();
						setAlignedEventPoints(cond,allPoints);
					}
	    		}
	    		
	    		//Differential event files
	    		if(manager.getNumConditions()>1 && runDiffTests){
	    			for(ExperimentCondition cond : manager.getConditions()){
		    			//Sort on the current condition
		    			BindingEvent.setSortingCond(cond);
		    			Collections.sort(events, new Comparator<BindingEvent>(){
		    	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
		    	        });
		    			
		    			for(ExperimentCondition othercond : manager.getConditions()){
		    				if(!cond.equals(othercond)){
				    			//Print diff events
				    			String condName = cond.getName(); 
				    			String othercondName = othercond.getName(); 
				    			condName = condName.replaceAll("/", "-");
				    			filename = filePrefix+"_"+condName+"_gt_"+othercondName+".diff.events";
								fout = new FileWriter(filename);
								fout.write(BindingEvent.conditionShortHeadString(cond)+"\n");
						    	for(BindingEvent e : events){
						    		double Q = e.getCondSigVCtrlP(cond);
						    		//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
						    		if(e.isFoundInCondition(cond) && Q <=qMinThres){
						    			if(e.getInterCondP(cond, othercond)<=diffPMinThres && e.getInterCondFold(cond, othercond)>0){
						    				fout.write(e.getConditionString(cond)+"\n");
						    			}
						    		}
						    	}
								fout.close();
		    				}
		    			}
		    		}
	    		}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
 
    /**
     * Print all binding events to files
     */
    public void writeFullEventFile(String filename){
    	if(events.size()>0){
    		try {
	    		//Full dataset table
		    	FileWriter fout = new FileWriter(filename);
		    	fout.write(BindingEvent.fullHeadString()+"\n");
		    	for(BindingEvent e : events)
		    		fout.write(e.toString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    
    /**
     * Print all motifs to file
     */
    public void writeMotifFile(String filename){
		try {
    		//Full dataset table
	    	FileWriter fout = new FileWriter(filename);
	    	for(ExperimentCondition cond : manager.getConditions()){
	    		if(getFreqMatrices(cond)!=null){
	    			for (WeightMatrix m : getFreqMatrices(cond))
	    				fout.write(WeightMatrix.printTransfacMatrix(m, cond.getName()));
	    		}
	    	}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }
 
    
    /**
     * Print replicate data counts to a file
     */
    public void writeReplicateCounts(String filename){
    	if(events.size()>0){
    		try {
	    		//Full dataset table
	    		FileWriter fout = new FileWriter(filename);
	    		fout.write(BindingEvent.repCountHeadString()+"\n");
	    		for(BindingEvent e : events)
	    			fout.write(e.getRepCountString()+"\n");
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
    }
    /**
     * Print all binding events to screen
     * TESTING ONLY
     */
    public void printBindingEvents(){
    	System.err.println(events.size()+" events found");
    	for(BindingEvent e : events){
    		System.err.println(e.toString());
    	}
    }
    
    /**
     * Convert binding events to components
     */
    public HashMap<Region, List<List<BindingSubComponents>>> getComponentsFromEnrichedEvents(){   	
    	HashMap<Region, List<List<BindingSubComponents>>> enrichedComps = new HashMap<Region,List<List<BindingSubComponents>>>();   
    	List<List<BindingEvent>> enrichedBindingEventbyCond = new ArrayList<List<BindingEvent>>();
    	for(ExperimentCondition cond : manager.getConditions())
    		enrichedBindingEventbyCond.add(new ArrayList<BindingEvent>());
    	for(BindingEvent e : events){
    		for(ExperimentCondition cond : manager.getConditions()){
    			double Q = e.getCondSigVCtrlP(cond);
    			if(e.isFoundInCondition(cond) && Q <=config.getMultiGPSQMinThres()){
    				if (!enrichedComps.containsKey(e.getContainingRegion()))
    					enrichedComps.put(e.getContainingRegion(), new ArrayList<List<BindingSubComponents>>());
    				enrichedBindingEventbyCond.get(cond.getIndex()).add(e);
    			}}}		
    	for (Region reg : enrichedComps.keySet())
    		for(ExperimentCondition cond : manager.getConditions())
    			enrichedComps.get(reg).add(new ArrayList<BindingSubComponents>());	
    	for(ExperimentCondition cond : manager.getConditions()){
    		for (BindingEvent e : enrichedBindingEventbyCond.get(cond.getIndex())){
    			BindingSubComponents currComp = new BindingSubComponents(e.getPoint(),cond.getReplicates().size());
    			currComp.setUnstrandedSumResponsibilities(e.getCondSigHits(cond));
    			enrichedComps.get(e.getContainingRegion()).get(cond.getIndex()).add(currComp);
    		}
    	}
		return enrichedComps;
    }    
}
