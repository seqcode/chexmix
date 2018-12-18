package org.seqcode.projects.chexmix.events;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;

import cern.jet.random.Binomial;
import cern.jet.random.ChiSquare;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

/**
 * Test the significance of count enrichment vs control
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class EnrichmentSignificance {

	protected EventsConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected double minFoldChange;
	protected double genomeLength;
	protected Binomial binomial;
	protected Poisson poisson;
	protected ChiSquare chisquare;
	
	public EnrichmentSignificance(EventsConfig con, ExperimentManager exptman, BindingManager bman, double minFoldChange, double genomeLength){
		this.config = con;
		this.manager = exptman;
		this.bindingManager = bman;
		this.minFoldChange = minFoldChange;
		this.genomeLength = genomeLength;
		binomial = new Binomial(100, .5, new DRand());
		poisson = new Poisson(1, new DRand());
		chisquare = new ChiSquare(1, new DRand());
	}

	/**
	 * Evaluate the significance of a set of EnrichedFeatures using Binomial and Poisson tests
	 * Assumes that the counts in the features are not already scaled
	 * @param modelRange
	 */
	public void execute() {this.execute(-1);}
	public void execute(int modelRange){

		//Calculate relative replicate weights using signal vs noise fractions in each signal channel
		double[] repWeights = new double[manager.getReplicates().size()];
		for(ExperimentCondition c : manager.getConditions()){
			double totalSig =0;
			for(ControlledExperiment r : c.getReplicates())
				totalSig += r.getSignalVsNoiseFraction() * r.getSignal().getHitCount();
			for(ControlledExperiment r : c.getReplicates())
				repWeights[r.getIndex()]=(r.getSignalVsNoiseFraction() * r.getSignal().getHitCount())/totalSig;
		}
		
		//Compute fold difference, log-likelihood, and p-values. 
		//To be clear about what each variable is measuring:
		// sigCtrlFold: weighted average of per-replicate fold difference between signal reads assigned to event and scaled control reads assigned to event.
		// sigCtrlP: Outcome of binomial test between sum of per-replicate signal read counts at event and sum of scaled per-replicate control counts at event. If multiple replicates use the same control, the counts are used redundantly since this is equivalent to summing the scaling factors across replicates. 
		// LL: log-likelihood loss if component was eliminated from model.
		// LLp: Chi-square distributed p-value corresponding to LL.  
		for (BindingEvent cf: bindingManager.getBindingEvents()){
			for(ExperimentCondition c1 : manager.getConditions()){
				double c1Sig = cf.getCondSigHitsFromReps(c1);
				double ctrlCountScaled = cf.getCondCtrlHitsScaledFromReps(c1);
				
				//Weighted fold difference, signal vs control
				double sigCtrlFold = 0;
				for(ControlledExperiment r : c1.getReplicates()){
					double repFold = cf.getRepCtrlHits(r)>1 ? cf.getRepSigHits(r)/(cf.getRepCtrlHits(r)*r.getControlScaling()) : cf.getRepSigHits(r);
					sigCtrlFold += repFold * repWeights[r.getIndex()];
				}
				
				//P-value, signal vs control
				double sigCtrlP = c1.getTotalControlCount()>=0 ? 
						evaluateSignificanceBinomial(c1Sig, ctrlCountScaled, cf.getCondTotalSigHitsFromReps(c1)) :
						evaluateSignificancePoisson(c1Sig, cf.getCondTotalSigHitsFromReps(c1), modelRange);
				cf.setCondSigVCtrlFold(c1, sigCtrlFold);
				cf.setCondSigVCtrlP(c1, sigCtrlP);
				
				//P-value, signal vs control for each replicate
				for(ControlledExperiment r : c1.getReplicates()){
					double repFold = cf.getRepCtrlHits(r)>1 ? cf.getRepSigHits(r)/(cf.getRepCtrlHits(r)*r.getControlScaling()) : cf.getRepSigHits(r);
					double repSigCtrlFold = repFold * repWeights[r.getIndex()];
					double repSigCtrlP = r.hasControl() ? 
							evaluateSignificanceBinomial(cf.getRepSigHits(r), cf.getRepCtrlHits(r)*r.getControlScaling(), r.getSignal().getHitCount()) :
							evaluateSignificancePoisson(cf.getRepSigHits(r), r.getSignal().getHitCount(), modelRange);
					cf.setRepSigVCtrlFold(r, repSigCtrlFold);
					cf.setRepSigVCtrlP(r, repSigCtrlP);
				}
				
				//Log-likelihood p-value
				if(config.CALC_EVENTS_LL){
					cf.setLLp(c1, chisquare.cdf(cf.getLLd(c1)));
					System.out.println(String.format("%s\t%s\t%.0f\t%e",cf.getPoint().getLocationString(),c1.getName(),cf.getLLd(c1),cf.getLLp(c1)));
				}
			}
		}
		// calculate q-values, correction for multiple testing
		benjaminiHochbergCorrection(bindingManager.getBindingEvents());
	}//end of evaluateConfidence method

	
	/**
	 * Evaluate the significance using Binomial 
	 */
	private double evaluateSignificanceBinomial(double countA, double countB, double total) {
        double pValue;
		
		if(countA+countB<=0 || (countA/countB)<=minFoldChange){
			return(1);
		}else{
	        try{
	
	            binomial.setNandP((int)Math.ceil(countA + countB), 1.0 / (minFoldChange + 1));
	            pValue = binomial.cdf((int)Math.ceil(countB));
	            
	        } catch(Exception err){
	            err.printStackTrace();
	            throw new RuntimeException(err.toString(), err);
	        }
	        return(pValue);
		}
	}

	/**
	 * Evaluate the significance using Poisson (only used when there are no controls) 
	 */
	private double evaluateSignificancePoisson(double countA, double total, int modelWidth) {
        double pValue;
		
		if(countA<=0){
			return(1);
		}else{
	        try{
	
	            poisson.setMean(minFoldChange * total * (double)modelWidth / (double)genomeLength );
	            int cA = (int)Math.ceil(countA);
	            pValue = 1 - poisson.cdf(cA) + poisson.pdf(cA);
	            
	        } catch(Exception err){
	            err.printStackTrace();
	            throw new RuntimeException(err.toString(), err);
	        }
	        return(pValue);
		}
	}

	/**
	 * Multiple hypothesis testing correction
	 */
	private void benjaminiHochbergCorrection(List<BindingEvent> features){
		double total = features.size();
		
		//Signal-vs-Control corrections by condition
		for(ExperimentCondition c : manager.getConditions()){
			BindingEvent.setSortingCond(c);		
			
			Collections.sort(features, new Comparator<BindingEvent>(){
	            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
	        });
			
			
			double rank =1.0;
			for(BindingEvent cf : features){
				cf.setCondSigVCtrlQ(c, Math.min(1.0, cf.getCondSigVCtrlP(c)*(total/rank)));
				rank++;
			}
			
			for(ControlledExperiment r : c.getReplicates()){
				Collections.sort(features, new Comparator<BindingEvent>(){
		            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareByRepSigCtrlPvalue(o2, r);}
		        });
				double rRank =1.0;
				for(BindingEvent cf : features){
					cf.setRepSigVCtrlQ(r, Math.min(1.0, cf.getRepSigVCtrlP(r)*(total/rRank)));
					rRank++;
				}
			}
		}
		
		//LL p-value corrections by condition
		if(config.CALC_EVENTS_LL){
			for(ExperimentCondition c : manager.getConditions()){
				BindingEvent.setSortingCond(c);
				Collections.sort(features, new Comparator<BindingEvent>(){
		            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareByLLPvalue(o2);}
		        });
				
				
				double rank =1.0;
				for(BindingEvent cf : features){
					cf.setCondSigVCtrlQ(c, Math.min(1.0, cf.getCondSigVCtrlP(c)*(total/rank)));
					for(ControlledExperiment r : c.getReplicates())
						cf.setRepSigVCtrlQ(r, Math.min(1.0, cf.getRepSigVCtrlP(r)*(total/rank)));
					rank++;
				}
			}
		}
		
		//Finally, sort on the first condition
		BindingEvent.setSortingCond(manager.getConditions().get(0));
		Collections.sort(features, new Comparator<BindingEvent>(){
            public int compare(BindingEvent o1, BindingEvent o2) {return o1.compareBySigCtrlPvalue(o2);}
        });
	}//end of benjaminiHochbergCorrection method
}
