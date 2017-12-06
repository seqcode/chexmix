package org.seqcode.projects.chexmix.events;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Pair;
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
	protected Map<ExperimentCondition, List<BindingSubtype>> bindingSubtypes;
	protected Map<ExperimentCondition, List<Boolean>> motifReverseStrands;
	protected Map<ExperimentCondition, Integer> numBindingType;
	protected Map<ControlledExperiment, BindingModel> unstrandedModel;
	protected Map<ExperimentCondition, Integer> maxInfluenceRange;
	protected Map<ExperimentCondition, Double> alpha;
	protected Map<ExperimentCondition, List<List<StrandedPoint>>> alignedEventPoints;
	protected Map<ExperimentCondition, List<BindingSubtype>> potentialBindingSubtypes;
	
	public BindingManager(EventsConfig con, ExperimentManager exptman){
		config = con;
		manager = exptman;
		events  = new ArrayList<BindingEvent>();
		conditionEvents = new HashMap<ExperimentCondition, List <BindingEvent>>();
		bindingSubtypes = new HashMap<ExperimentCondition, List<BindingSubtype>>();
		motifReverseStrands = new HashMap<ExperimentCondition, List<Boolean>>();
		maxInfluenceRange = new HashMap<ExperimentCondition, Integer>();
		alpha = new HashMap<ExperimentCondition, Double>();
		numBindingType = new HashMap<ExperimentCondition, Integer>();
		unstrandedModel = new HashMap<ControlledExperiment, BindingModel>();
		alignedEventPoints = new HashMap<ExperimentCondition, List<List<StrandedPoint>>>();
		for(ExperimentCondition cond : manager.getConditions()){
			conditionEvents.put(cond, new ArrayList<BindingEvent>());
			bindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
			numBindingType.put(cond, 1);
			maxInfluenceRange.put(cond,0);
			alpha.put(cond, 0.0);
			alignedEventPoints.put(cond, new ArrayList<List<StrandedPoint>>());
			potentialBindingSubtypes.put(cond, new ArrayList<BindingSubtype>());
		}
	}
	
	public List<BindingEvent> getBindingEvents(){return events;}
	public List<BindingEvent> getConditionBindingEvents(ExperimentCondition ec){return conditionEvents.get(ec);}
	public List<BindingSubtype> getBindingSubtype(ExperimentCondition ec){return bindingSubtypes.get(ec);}
	public List<Boolean> getReverseMotifs(ExperimentCondition ec){return motifReverseStrands.get(ec);}
	public Integer getNumBindingType(ExperimentCondition ec){return numBindingType.get(ec);}
	public BindingModel getUnstrandedBindingModel(ControlledExperiment ce){return unstrandedModel.get(ce);}
	public Integer getMaxInfluenceRange(ExperimentCondition ec){return maxInfluenceRange.get(ec);}
	public Double getAlpha(ExperimentCondition ec){return alpha.get(ec);}
	public List<List<StrandedPoint>> getAlignedEventPoints(ExperimentCondition ec){return alignedEventPoints.get(ec);}
	public List<BindingSubtype> getPotentialBindingSubtypes(ExperimentCondition ec){return potentialBindingSubtypes.get(ec);}

	public void setBindingEvents(List<BindingEvent> e){events =e;}
	public void setConditionBindingEvents(ExperimentCondition ec, List<BindingEvent> e){conditionEvents.put(ec, e);}
	public void setBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> sub){bindingSubtypes.put(ec, sub); numBindingType.put(ec, sub.size());}
	public void setReverseStrand(ExperimentCondition ec, List<Boolean> reverse){motifReverseStrands.put(ec, reverse);}
	public void setAlpha(ExperimentCondition ec, Double a){alpha.put(ec,a);}
	public void setUnstrandedBindingModel(ControlledExperiment ce, BindingModel mod){unstrandedModel.put(ce, mod);}
	public void setAlignedEventPoints(ExperimentCondition ec, List<List<StrandedPoint>> points){alignedEventPoints.put(ec, points);}
	public void addPotentialBindingSubtypes(ExperimentCondition ec, List<BindingSubtype> subtypes){potentialBindingSubtypes.get(ec).addAll(subtypes);}
	public void clearPotentialBindingSubtypes(ExperimentCondition ec){ potentialBindingSubtypes.put(ec, new ArrayList<BindingSubtype>());}

	public void updateMaxInfluenceRange(ExperimentCondition ec, boolean firstround){
		int max=0; 
		if (firstround && config.getDefaultBindingModel()==null){		
			for(ControlledExperiment rep : ec.getReplicates())
				if (getUnstrandedBindingModel(rep).getInfluenceRange()>max)
					max=getUnstrandedBindingModel(rep).getInfluenceRange();
		}else{
			int min=Integer.MAX_VALUE;
			for ( BindingSubtype subtype : getBindingSubtype(ec)){
				TagProbabilityDensity density =subtype.getBindingModel(0);
				if (density.getInfluenceRange()< min)
					min=density.getInfluenceRange();		
			}
			max=min;
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
					
					//Print aligned points
					List<List<StrandedPoint>> subtypePoints = new ArrayList<List<StrandedPoint>>();
					List<List<Region>> potReg = new ArrayList<List<Region>>();
					List<List<StrandedPoint>> alignedPoints = new ArrayList<List<StrandedPoint>>();
					for (int bt=0; bt < getNumBindingType(cond); bt++){
						subtypePoints.add(new ArrayList<StrandedPoint>());
						potReg.add(new ArrayList<Region>());
						alignedPoints.add(new ArrayList<StrandedPoint>());
					}
					for(BindingEvent e : events){
						double Q = e.getCondSigVCtrlP(cond);
						//Because of the ML step and component sharing, I think that an event could be assigned a significant number of reads without being "present" in the condition's EM model.
						if(e.isFoundInCondition(cond) && Q <=qMinThres){
							Pair<Integer, StrandedPoint>p=e.getMaxSubtypePoint(cond);
							subtypePoints.get(p.car()).add(p.cdr());
							potReg.get(p.car()).add(e.getContainingRegion());
						}
					}	
						
					// Align points
					for (int bt=0;bt< getNumBindingType(cond); bt++){
						List<StrandedPoint> points = subtypePoints.get(bt);
						BindingSubtype currSubtype = bindingSubtypes.get(cond).get(bt);
						int c=0;
						if (currSubtype.hasMotif()){
		    				int offset = currSubtype.getMotifOffset();
		    				if (currSubtype.reverseMotif()){ //reverse complement
		    					for (StrandedPoint p : points){
		    						int location = p.getStrand()=='+' ? p.getLocation()-offset : p.getLocation()+offset;
		    						// Filter out sequences that are on the edge of cashed sequences
									if((location-potReg.get(bt).get(c).getStart()>config.SEQPLOTWIN) && (potReg.get(bt).get(c).getEnd()-location>config.SEQPLOTWIN))
										alignedPoints.get(bt).add(new StrandedPoint(p.getGenome(),p.getChrom(),location,p.getStrand() =='+' ? '-' : '+'));
									c++;
		    					}	    						
		    				}else{
		    					for (StrandedPoint p : points){
		    						int location = p.getStrand()=='+' ? p.getLocation()+offset : p.getLocation()-offset;
		    						if((location-potReg.get(bt).get(c).getStart()>config.SEQPLOTWIN) && (potReg.get(bt).get(c).getEnd()-location>config.SEQPLOTWIN))
		    							alignedPoints.get(bt).add(new StrandedPoint(p.getGenome(),p.getChrom(),location,p.getStrand()));
		    						c++;
		    					}	    						
		    				}	    					
						}else{ // no motif found
							for (StrandedPoint p : points){
								if((p.getLocation()-potReg.get(bt).get(c).getStart()>config.SEQPLOTWIN) && (potReg.get(bt).get(c).getEnd()-p.getLocation()>config.SEQPLOTWIN))
									alignedPoints.get(bt).add(p);
								c++;
							}
						}
					}
					filename = filePrefix+"_"+condName+".subtype.aligned.events";
					fout = new FileWriter(filename);
					int subtypeC=0;
					for (List<StrandedPoint> points : alignedPoints){
						for (StrandedPoint p : points)
							fout.write(p.toString()+"\tSubtype"+subtypeC+'\n');
						subtypeC++;
					}
					fout.close();
					setAlignedEventPoints(cond,alignedPoints);
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
	    	for(ExperimentCondition cond : manager.getConditions())
	    		for (BindingSubtype sub : getBindingSubtype(cond))
	    			if (sub.hasMotif())
	    				fout.write(WeightMatrix.printTransfacMatrix(sub.getFreqMatrix(), cond.getName()));
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
    public HashMap<Region, List<List<BindingSubComponents>>> getComponentsFromEnrichedEvents(List<Region> potRegions){      	
    	HashMap<Region, List<List<BindingSubComponents>>> enrichedComps = new HashMap<Region,List<List<BindingSubComponents>>>();  
    	for (Region reg : potRegions){
    		enrichedComps.put(reg, new ArrayList<List<BindingSubComponents>>());
    		for(ExperimentCondition cond : manager.getConditions())
    			enrichedComps.get(reg).add(new ArrayList<BindingSubComponents>());
    	}
    	for(BindingEvent e : events){
    		for(ExperimentCondition cond : manager.getConditions()){
    			double Q = e.getCondSigVCtrlP(cond);
    			if(e.isFoundInCondition(cond) && Q <=config.getMultiGPSQMinThres()){
    				BindingSubComponents currComp = new BindingSubComponents(e.getPoint(),cond.getReplicates().size());
        			currComp.setUnstrandedSumResponsibilities(e.getCondSigHits(cond));
    				enrichedComps.get(e.getContainingRegion()).get(cond.getIndex()).add(currComp);
    			}}}		
		return enrichedComps;
    }    
}
