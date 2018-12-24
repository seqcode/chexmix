package org.seqcode.projects.chexmix.events;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExperimentScaler;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

/** 
 * This class assesses replication of each peak across replicates in a given condition.
 * Replication codes (added to BindingEvents):
 * 		-2 : This event was not called in this condition  
 * 		 0 : This condition only contains a single replicate
 * 		-1 : Called in one replicate or condition as a whole, significant difference detected across some replicates
 * 		 1 : Called in one replicate or condition as a whole, no significant difference detected across any replicates
 * 		>1 : Number of replicates in which binding event is called significant
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
	protected double minFoldChange;
	protected double qMinThres;
	
	public ReplicationTester(EventsConfig con, ExperimentManager exptman, BindingManager bman, double minFoldChange, double qMinThres){
		this.config = con;
		this.manager = exptman;
		this.bindingManager = bman;
		binomial = new Binomial(100, .5, new DRand());
		scaler = new ExperimentScaler();
		this.minFoldChange = minFoldChange;
		this.qMinThres = qMinThres;
	}
	
	public void execute(){
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
						repScaling[r2.getIndex()][r1.getIndex()]=1/s;
					}
				}
			}
		}
		
		// 2) assign replication codes
		for(ExperimentCondition c : manager.getConditions()){
			
			// 2.1) if the condition contains only one replicate, code = 0
			if(c.getReplicates().size()==1){
				for(BindingEvent e : events)
					e.setReplicationCode(c, 0);
			}else{

				for(BindingEvent e : events){
					// 2.2) if this event is not called in this condition, code = -2
					if(!e.isFoundInCondition(c))
						e.setReplicationCode(c, -2);
					
					// 2.3) if this event passes per-replicate significant filters in multiple replicates, code = number of passing replicates
					int repC=0;
					for(ControlledExperiment r : c.getReplicates())
						if(e.isFoundInCondition(c) && e.getRepSigVCtrlQ(r)<=qMinThres)
							repC++;
					if(repC>1){
						e.setReplicationCode(c, repC);
					}else if(repC==1 || e.getCondSigVCtrlQ(c)<=qMinThres){
						// 2.4) if the event is called in only one replicate or the condition as a whole, test the difference between replicates using Binomial
						boolean consistent=true;
						for(ControlledExperiment r1 : c.getReplicates()){
							for(ControlledExperiment r2 : c.getReplicates()){
								if(r1.getIndex()!=r2.getIndex()){
									double countA = e.getRepSigHits(r1), countB = e.getRepSigHits(r2)*repScaling[r1.getIndex()][r2.getIndex()];
									int cAB = (int)Math.ceil(countA + countB);
									if(cAB>0){
										binomial.setNandP((int)Math.ceil(countA + countB), 1.0 / (minFoldChange + 1));
										consistent = consistent && (binomial.cdf((int)Math.ceil(countB))>qMinThres);
									}
								}
							}
						}
						// 2.4 cont) if there is any significant difference, code = -1. otherwise, code = 1
						e.setReplicationCode(c, consistent ? 1: -1);
					}
				}
			}
		}
	}
	
	/**
	 * Write a file with the replication codes for every binding event in every experiment
	 */
	public void writeReplicationInfoFile(String filename){
		List<BindingEvent> events = bindingManager.getBindingEvents();
		try {
    		//Full output table (all non-zero components)
    		FileWriter fout = new FileWriter(filename);
    		String head = "### ChExMix replication codes\n"+
    				"# Codes:\n"+
    				"#    -2 : Event was not called in this condition\n" +  
    				"#     0 : This condition only contains a single replicate\n"+
    				"#    -1 : Called in one replicate or condition as a whole, significant difference detected across some replicates\n"+
    				"#     1 : Called in one replicate or condition as a whole, no significant difference detected across any replicates\n"+
    				"#    >1 : Number of replicates in which binding event is called significant\n"+
    				"#\n"+
    				"#BindingEvent";
    		for(ExperimentCondition c : manager.getConditions())
    			head = head +"\t"+c.getName();
    		head = head +"\n";
			fout.write(head);

    		for(BindingEvent e : events){
    			fout.write(e.getPoint().getLocationString());
    			for(ExperimentCondition c : manager.getConditions())
    				fout.write("\t"+e.getReplicationCode(c));
    			fout.write("\n");
    		}
			fout.close();
    		
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
