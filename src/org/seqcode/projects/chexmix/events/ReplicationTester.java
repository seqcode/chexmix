package org.seqcode.projects.chexmix.events;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExperimentScaler;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

/** 
 * This class assesses reproducibility of each peak across replicates in a given condition.
 * A binding event is called "replicated" if:
 *  1) it passes significance filters (vs. background) in multiple replicates
 *  2) it is not significantly difference according to a pairwise Binomial after regression-based scaling
 *  3) the condition only has one replicate
 *  
 * @author mahony
 *
 */
public class ReplicationTester {

	protected EventsConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected Binomial binomial;
	protected ExperimentScaler scaler;
	
	public ReplicationTester(EventsConfig con, ExperimentManager exptman, BindingManager bman){
		this.config = con;
		this.manager = exptman;
		this.bindingManager = bman;
		binomial = new Binomial(100, .5, new DRand());
		scaler = new ExperimentScaler();
	}
	
	public void execute(double qMinThres){
		List<BindingEvent> events = bindingManager.getBindingEvents();
				
		// 1) calculate scaling factors using regression on the binding event read counts
		double repScaling[][] = new double [manager.getReplicates().size()][manager.getReplicates().size()];
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment r1 : c.getReplicates()){
				for(ControlledExperiment r2 : c.getReplicates()){
					if(r1.getIndex()==r2.getIndex()){repScaling[r1.getIndex()][r2.getIndex()]=1;}
					else if(r1.getIndex()<r2.getIndex()){
						List<Float> r1Counts = new ArrayList<Float>();
						List<Float> r2Counts = new ArrayList<Float>();
						for(BindingEvent e : events){
							if(e.isFoundInCondition(c)){
								r1Counts.add((float)e.getRepSigHits(r1));
								r2Counts.add((float)e.getRepSigHits(r2));
							}
						}
						double s = scaler.scalingRatioByRegression(r1Counts, r2Counts);
						repScaling[r1.getIndex()][r2.getIndex()]=s;
						repScaling[r2.getIndex()][r1.getIndex()]=s;
					}
				}
			}
		}
		
		// 2) implement replication rules
		for(ExperimentCondition c : manager.getConditions()){
			
			// 2.1) if the condition contains only one replicate, all events are "replicated"
			if(c.getReplicates().size()==1){
				for(BindingEvent e : events)
					if(e.isFoundInCondition(c))
						e.setReplicated(c, true);
			}else{
				// 2.2) if the event has passed significance filter in multiple replicates, it is "replicated"
				for(BindingEvent e : events){
					int repC=0;
					for(ControlledExperiment r : c.getReplicates())
						if(e.isFoundInCondition(c) && e.getRepSigVCtrlQ(r)<=qMinThres)
							repC++;
					if(repC>1)
						e.setReplicated(c, true);
				}
				// 2.3) test the difference between replicates using Binomial
				
			}
		}
		
		
		
	}
	
}
