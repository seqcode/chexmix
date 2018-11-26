package org.seqcode.projects.chexmix.events;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.location.AnnotationLoader;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.Pair;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;

/**
 * BindingEvent: a class representing a potential binding event and associated read counts and 
 * significance P-values across experimental conditions and replicates. 
 * 
 * BindingEvents are associated with a particular configuration of binding components, since the read counts per event can depend on 
 * other proximal binding events. Therefore, and depending on what configuration the binding event comes from, multiple binding events 
 * may share the same coordinate without having the same read counts associated per experiment. While this may seem strange, defining 
 * events in this way reduces redundant records while representing all discovered configurations.   
 * The conditions that shared the exact configuration that led to the creation of this binding event are represented in foundInCond.
 * 
 * @author Naomi Yamada
 *
 */
public class BindingEvent implements Comparable<BindingEvent>{
	
	protected static ExperimentManager experiments=null;
	protected static EventsConfig config=null;
	protected static ChExMixConfig gpsconfig=null;
	protected static ExperimentCondition sortingCondition = null;
	protected static ExperimentCondition sortingConditionB = null;
	protected static final int numSingleCondCols = 4; //Number of columns in the output file for each single-condition 
	protected static final int numInterCondCols = 3; //Number of columns in the output file for each inter-condition 
	protected static boolean findMotif=true;
	protected static int[] numBindingTypes=null;
	protected Point point;
	protected StrandedPoint[][][] typePoints;
	protected double[][][] typeProbs;
	protected Region containingReg;
	protected boolean [] foundInCond; //The binding event was discovered in these conditions [indexed by condition]
	protected double [] condSigHits;  //Signal counts by condition  (not scaled)      [indexed by condition]
	protected double [] condCtrlHits;  //Control counts by condition  (not scaled)    [indexed by condition]
	protected double [] condSigVCtrlFold;  //Signal vs Control fold by condition      [indexed by condition]
	protected double [] condSigVCtrlP;  //Signal vs Control P-value by condition      [indexed by condition]
	protected double [] repSigHits;  //Signal counts by replicate  (not scaled)     [indexed by replicate]
	protected double [] repCtrlHits;  //Control counts by replicate  (not scaled)   [indexed by replicate]
	protected double [] repSigVCtrlFold;  //Signal vs Control fold by replicate     [indexed by replicate]
	protected double [] repSigVCtrlP;  //Signal vs Control P-value by replicate     [indexed by replicate]
	protected double [][] interCondScMean;   //Signal vs signal scaled mean (logCPM from EdgeR), inter-condition  [indexed by condition & condition]
	protected double [][] interCondFold;   //Signal vs signal fold difference  (logFold from EdgeR), inter-condition  [indexed by condition & condition]
	protected double [][] interCondP;   //Signal vs signal P, inter-condition         [indexed by condition & condition]
	protected double [][] interRepP;   //Signal vs signal P, inter-replicate        [indexed by replicate & replicate]
	protected double []   LLd; 			//Log-likelihood loss test statistic resulting from eliminating component [indexed by condition]
	protected double []   LLp;			//P-value for LL [indexed by condition]
	protected double [][][] motifScores;	//LL score for motif match at binding event
	protected String sequences[][][];   //Sequences at binding event	
	protected Gene nearestGene=null;
	protected int distToGene=0;
	
    
 	public BindingEvent(Point p, Region potentialReg){
		point=p;
		containingReg=potentialReg;
			
		int numC = experiments.getConditions().size();
		int numR = experiments.getReplicates().size();
		typePoints = new StrandedPoint [numC][][];
		typeProbs = new double[numC][][];
		foundInCond = new boolean[numC];
		condSigHits = new double [numC];
		condCtrlHits = new double [numC];
		condSigVCtrlFold = new double [numC];
		condSigVCtrlP = new double [numC];
		interCondScMean = new double [numC][numC];
		interCondFold = new double [numC][numC];
		interCondP = new double [numC][numC];
		repSigHits = new double [numR];
		repCtrlHits = new double [numR];
		repSigVCtrlFold = new double [numR];
		repSigVCtrlP = new double [numR];
		interRepP = new double [numR][numR];
		LLd = new double [numC];
		LLp = new double [numC];
		motifScores = new double[numC][][];
		sequences = new String[numC][][];
		for(int c=0; c<numC; c++){
			foundInCond[c]=false; condSigHits[c]=0; condCtrlHits[c]=0; condSigVCtrlP[c]=1; condSigVCtrlFold[c]=1; 
//			motifScores[c]=0; sequences[c]="--";
			for(int d=0; d<numC; d++){
				interCondScMean[c][d]=0; interCondFold[c][d]=1; interCondP[c][d]=1;
			}
		}
		for(int r=0; r<numR; r++){
			repSigHits[r]=0; repCtrlHits[r]=0; repSigVCtrlP[r]=1; repSigVCtrlFold[r]=1;
		}
		findMotif=gpsconfig.getFindingMotifs();
	}
 	
 	
 	//Getters
 	public StrandedPoint[][] getTypePoints(ExperimentCondition c){return typePoints[experiments.getConditionIndex(c)];}
 	public double[][] getTypeProbs(ExperimentCondition c){return typeProbs[experiments.getConditionIndex(c)];}
 	public boolean isFoundInCondition(ExperimentCondition c){return foundInCond[experiments.getConditionIndex(c)];}
 	public boolean isFoundInCondition(int conditionIndex){return foundInCond[conditionIndex];}
 	public double getCondSigHits(ExperimentCondition c){return(condSigHits[experiments.getConditionIndex(c)]);}
 	public double getCondCtrlHits(ExperimentCondition c){return(condCtrlHits[experiments.getConditionIndex(c)]);}
 	public double getCondSigVCtrlP(ExperimentCondition c){return(condSigVCtrlP[experiments.getConditionIndex(c)]);}
 	public double getCondSigVCtrlFold(ExperimentCondition c){return(condSigVCtrlFold[experiments.getConditionIndex(c)]);}
 	public double getRepSigHits(ControlledExperiment r){return(repSigHits[r.getIndex()]);}
 	public double getRepCtrlHits(ControlledExperiment r){return(repCtrlHits[r.getIndex()]);}
 	public double getRepSigVCtrlP(ControlledExperiment r){return(repSigVCtrlP[r.getIndex()]);}
 	public double getRepSigVCtrlFold(ControlledExperiment r){return(repSigVCtrlFold[r.getIndex()]);}
 	public double getInterCondScMean(ExperimentCondition c1, ExperimentCondition c2){return(interCondScMean[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
 	public double getInterCondFold(ExperimentCondition c1, ExperimentCondition c2){return(interCondFold[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
 	public double getInterCondP(ExperimentCondition c1, ExperimentCondition c2){return(interCondP[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]);}
	public double getInterRepP(ControlledExperiment r1, ControlledExperiment r2){return interRepP[r1.getIndex()][r2.getIndex()];}
	public static int getNumSingleCondCols(){return numSingleCondCols;}
	public static int getNumInterCondCols(){return numInterCondCols;}
	public double getLLd(ExperimentCondition c1){return(LLd[experiments.getConditionIndex(c1)]);}
	public double getLLp(ExperimentCondition c1){return(LLp[experiments.getConditionIndex(c1)]);}
	public double[][] getMotifScore(ExperimentCondition c1){return(motifScores[experiments.getConditionIndex(c1)]);}
	public String[][] getSequence(ExperimentCondition c1){return sequences[experiments.getConditionIndex(c1)];}
	public String getFoundInConditionsString(){
		String foundStr="";
		for(ExperimentCondition c : experiments.getConditions()){
			String found = foundInCond[c.getIndex()] ? "1" : "0";
			if(foundStr.length()>0)
				foundStr=foundStr+","+found;
			else
				foundStr=found;
		}
		return foundStr;
	}
	public double getCondSigHitsFromReps(ExperimentCondition c){
		double hits=0;
		for(ControlledExperiment rep : c.getReplicates())
			hits+=repSigHits[rep.getIndex()];
		return hits;
	}
	public double getCondCtrlHitsScaledFromReps(ExperimentCondition c){
		double hits=0;
		for(ControlledExperiment rep : c.getReplicates())
			hits+=repCtrlHits[rep.getIndex()]*rep.getControlScaling();
		return hits;
	}
	public double getCondTotalSigHitsFromReps(ExperimentCondition c){
		double total=0;
		for(ControlledExperiment rep : c.getReplicates())
			total+=rep.getSignal().getHitCount();
		return total;
	}
	public double getCondTotalCtrlHitsScaledFromReps(ExperimentCondition c){
		double total=0;
		for(ControlledExperiment rep : c.getReplicates())
			total+=rep.getControl().getHitCount()*rep.getControlScaling();
		return total;
	}
	
	//Setters
	public void setTypePoints(ExperimentCondition c, StrandedPoint[][] points){ typePoints[experiments.getConditionIndex(c)]=points;}
	public void setTypeProbs(ExperimentCondition c, double[][] probs) { typeProbs[experiments.getConditionIndex(c)]=probs;}
	public void setIsFoundInCondition(ExperimentCondition c, boolean found){foundInCond[experiments.getConditionIndex(c)] = found;}
	public void setIsFoundInCondition(int c, boolean found){foundInCond[c] = found;}
	public void setCondSigHits(ExperimentCondition c, double x){condSigHits[experiments.getConditionIndex(c)]=x;}
	public void setCondCtrlHits(ExperimentCondition c, double x){condCtrlHits[experiments.getConditionIndex(c)]=x;}
	public void setCondSigVCtrlFold(ExperimentCondition c, double x){condSigVCtrlFold[experiments.getConditionIndex(c)]=x;}
	public void setCondSigVCtrlP(ExperimentCondition c, double x){condSigVCtrlP[experiments.getConditionIndex(c)]=x;}
	public void setLLd(ExperimentCondition c, double l){LLd[experiments.getConditionIndex(c)]=l;}
	public void setLLp(ExperimentCondition c, double p){LLp[experiments.getConditionIndex(c)]=p;}
	public void setRepSigHits(ControlledExperiment r, double x){repSigHits[r.getIndex()]=x;}
	public void setRepCtrlHits(ControlledExperiment r, double x){repCtrlHits[r.getIndex()]=x;}
	public void setRepSigVCtrlFold(ControlledExperiment r, double x){repSigVCtrlFold[r.getIndex()]=x;}
	public void setRepSigVCtrlP(ControlledExperiment r, double x){repSigVCtrlP[r.getIndex()]=x;}
	public void setInterCondScMean(ExperimentCondition c1, ExperimentCondition c2, double x){interCondScMean[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterCondFold(ExperimentCondition c1, ExperimentCondition c2, double x){interCondFold[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterCondP(ExperimentCondition c1, ExperimentCondition c2, double x){interCondP[experiments.getConditionIndex(c1)][experiments.getConditionIndex(c2)]=x;}
	public void setInterRepP(ControlledExperiment r1, ControlledExperiment r2, double x){interRepP[r1.getIndex()][r2.getIndex()]=x;}
	public void setMotifScore(ExperimentCondition c1, double[][] s){motifScores[experiments.getConditionIndex(c1)] = s;}
	public void setSequence(ExperimentCondition c1, String[][] seq){sequences[experiments.getConditionIndex(c1)]=seq;}
	
	/**
	 * Rank according to increasing p-value (score) for the sorting condition, then by signal count for the sorting condition
	 */
	public int compareBySigCtrlPvalue(BindingEvent p) {
		if(this.getCondSigVCtrlP(sortingCondition)<p.getCondSigVCtrlP(sortingCondition)){return(-1);}
		else if(this.getCondSigVCtrlP(sortingCondition)>p.getCondSigVCtrlP(sortingCondition)){return(1);}
		else{
			if(this.getCondSigHits(sortingCondition)>p.getCondSigHits(sortingCondition)){return(-1);}
			else if(this.getCondSigHits(sortingCondition)<p.getCondSigHits(sortingCondition)){return(1);}
		}return(0);	
	}
	
	/**
	 * Rank according to increasing replicate p-value (score) for the sorting condition, then by signal count for the sorting condition
	 */
	public int compareByRepSigCtrlPvalue(BindingEvent p, ControlledExperiment rep) {
		if(this.getRepSigVCtrlP(rep)<p.getRepSigVCtrlP(rep)){return(-1);}
		else if(this.getRepSigVCtrlP(rep)>p.getRepSigVCtrlP(rep)){return(1);}
		else{
			if(this.getRepSigHits(rep)>p.getRepSigHits(rep)){return(-1);}
			else if(this.getRepSigHits(rep)<p.getRepSigHits(rep)){return(1);}
		}return(0);	
	}
	
	/**
	 * Rank according to increasing p-value (score) for the pair of sorting conditions, then by signal count for the sorting condition
	 */
	public int compareByInterCondPvalue(BindingEvent p) {
		if(this.getInterCondP(sortingCondition, sortingConditionB) < p.getInterCondP(sortingCondition, sortingConditionB)){return(-1);}
		else if(this.getInterCondP(sortingCondition, sortingConditionB) > p.getInterCondP(sortingCondition, sortingConditionB)){return(1);}
		else{
			if(this.getCondSigHits(sortingCondition)>p.getCondSigHits(sortingCondition)){return(-1);}
			else if(this.getCondSigHits(sortingCondition)<p.getCondSigHits(sortingCondition)){return(1);}
		}return(0);
	}
	
	/**
	 * Rank according to increasing LL p-value for the sorting condition, then by signal count for the sorting condition
	 */
	public int compareByLLPvalue(BindingEvent p) {
		if(this.getLLp(sortingCondition)<p.getLLp(sortingCondition)){return(-1);}
		else if(this.getLLp(sortingCondition)>p.getLLp(sortingCondition)){return(1);}
		else{
			if(this.getCondSigHits(sortingCondition)>p.getCondSigHits(sortingCondition)){return(-1);}
			else if(this.getCondSigHits(sortingCondition)<p.getCondSigHits(sortingCondition)){return(1);}
		}return(0);	
	}
	
	/**
	 * Rank according to decreasing signal count for the sorting condition
	 */
	public int compareBySigCounts(BindingEvent p){
		if(this.getCondSigHits(sortingCondition)>p.getCondSigHits(sortingCondition)){return(-1);}
		else if(this.getCondSigHits(sortingCondition)<p.getCondSigHits(sortingCondition)){return(1);}
		return(0);
	}
	
	/**
	 * Rank according to location
	 * @param f
	 * @return
	 */
	public int compareByLocation(BindingEvent f) {
		return getPoint().compareTo(f.getPoint());
	}
	
	//Comparable default method
	public int compareTo(BindingEvent f) {
		return compareByLocation(f);
	}
	
	public Point getPoint(){return(point);}
	public Region getContainingRegion(){return(containingReg);}
	 
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of BindingEvent.toString().
	 */
	public static String fullHeadString(){
		String head="### ChExMix output\n";
		
		head = head + "#Condition\tName\tIndex\tTotalSigCount\tSignalFraction\n";
		//Add some basic information on the experiments
		for(ExperimentCondition c : experiments.getConditions()){
			head = head + "#Condition\t"+c.getName()+"\t"+c.getIndex()+"\t"+c.getTotalSignalCount()+"\t"+String.format("%.3f",c.getTotalSignalVsNoiseFrac())+"\n";
		}
		head = head + "#Replicate\tParentCond\tName\tIndex\tSigCount\tCtrlCount\tSigCtrlScaling\tSignalFraction\n";
		for(ExperimentCondition c : experiments.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				if(r.getControl()==null)
					head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t0\t1\t"+String.format("%.2f",r.getSignalVsNoiseFraction())+"\n";
				else
					head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t"+r.getControl().getHitCount()+"\t"+String.format("%.4f",r.getControlScaling())+"\t"+String.format("%.3f",r.getSignalVsNoiseFraction())+"\n";
			}
		}
		
		head = head + "#\n#Point";
		for(ExperimentCondition c : experiments.getConditions())
			head = head +"\t"+c.getName()+"_Sig"+"\t"+c.getName()+"_Ctrl"+"\t"+c.getName()+"_log2Fold"+"\t"+c.getName()+"_log2P";
		for(ExperimentCondition c : experiments.getConditions())
			for(ExperimentCondition c2 : experiments.getConditions()){
				if(c!=c2){
					head = head +"\t"+c.getName()+"_vs_"+c2.getName()+"_log2CPM"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2Fold"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2P";
				}
			}
		head = head+"\tActiveConds";
		if(config.isAddingAnnotations())
			head = head +"\tnearestGene\tdistToGene";
		
		//head = head +"\n";
		
		return head;
	}
	
	/**
	 * Print info on each condition and inter-condition
	 */
	public String toString() {
		//String gene = nearestGene == null ? "NONE" : nearestGene.getName();
		String out = point.getLocationString();
		for(ExperimentCondition c : experiments.getConditions()){
			double logP = Math.log(getCondSigVCtrlP(c))/config.LOG2;
			double logF = Math.log(getCondSigVCtrlFold(c))/config.LOG2;
			out = out+"\t"+String.format("%.1f", getCondSigHits(c))+"\t"+String.format("%.1f", getCondCtrlHits(c))+"\t"+String.format("%.3f", logF)+"\t"+String.format("%.3f", logP); 
		}for(ExperimentCondition c : experiments.getConditions())
			for(ExperimentCondition c2 : experiments.getConditions()){
				if(c!=c2){
					double logP = Math.log(getInterCondP(c, c2))/config.LOG2; 
					out = out+"\t"+String.format("%.2f", getInterCondScMean(c, c2))+"\t"+String.format("%.3f", getInterCondFold(c,c2))+"\t"+String.format("%.3f",logP);
				}
			}
		out = out+"\t"+this.getFoundInConditionsString();
		if(config.isAddingAnnotations())
			if(nearestGene==null)
				out = out + "\tNONE\t"+config.getMaxAnnotDistance();
			else
				out = out + "\t"+nearestGene.getName()+"\t"+distToGene;
		//out = out+"\n";
		return out;
	}
	
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of BindingEvent.getRepCountString().
	 */
	public static String repCountHeadString(){
		/* Use a simple header instead of the more complex format below
		String head="### ReplicateCount output\n";
		head = head + "#Condition\tName\tIndex\tTotalSigCount\tSignalFraction\n";
		//Add some basic information on the experiments
		for(ExperimentCondition c : experiments.getConditions()){
			head = head + "#Condition\t"+c.getName()+"\t"+c.getIndex()+"\t"+c.getTotalSignalCount()+"\t"+c.getSigProp()+"\n";
		}
		head = head + "#Replicate\tParentCond\tName\tIndex\tSigCount\tCtrlCount\tSignalFraction\n";
		for(ExperimentCondition c : experiments.getConditions()){
			for(ControlledExperiment r : c.getReplicates()){
				if(r.getControl()==null)
					head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t0\t1\t"+String.format("%.3f",r.getSigProp())+"\n";
				else
					head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t"+r.getControl().getHitCount()+"\t"+String.format("%.2f",r.getControlScaling())+"\t"+String.format("%.3f",r.getSigProp())+"\n";
			}
		}
		head = head + "#\n#Point";
		*/
		String head = "Point";
		for(ControlledExperiment r : experiments.getReplicates())
			head = head +"\t"+r.getName();
			
		return head;
	}
	
	/**
	 * Print info on each replicate's & condition's unnormalized counts
	 */
	public String getRepCountString() {
		String out = point.getLocationString();
		
		for(ControlledExperiment r : experiments.getReplicates())
			out = out+"\t"+String.format("%.0f", getRepSigHits(r)); 
		
		return out;
	}
	
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of BindingEvent.getConditionString().
	 * (i.e. single condition files)
	 */
	public static String conditionHeadString(ExperimentCondition c){
		String head="### ChExMix output\n";
		
		head = head + "#Condition\tName\tIndex\tTotalSigCount\tSignalFraction\n";
		head = head + "#Condition\t"+c.getName()+"\t"+c.getIndex()+"\t"+c.getTotalSignalCount()+"\t"+String.format("%.3f",c.getTotalSignalVsNoiseFrac())+"\n";
		
		head = head + "#Replicate\tParentCond\tName\tIndex\tSigCount\tCtrlCount\tCtrlScaling\tSignalFraction\n";
		for(ControlledExperiment r : c.getReplicates()){
			if(r.getControl()==null)
				head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t0\t1\t"+String.format("%.3f",r.getSignalVsNoiseFraction())+"\n";
			else
				head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t"+r.getControl().getHitCount()+"\t"+String.format("%.3f",r.getControlScaling())+"\t"+String.format("%.3f",r.getSignalVsNoiseFraction())+"\n";
		}
		
		head = head + "#\n#Point";
		head = head +"\t"+c.getName()+"_Sig"+"\t"+c.getName()+"_Ctrl"+"\t"+c.getName()+"_log2Fold"+"\t"+c.getName()+"_log2P";
		
		for(ExperimentCondition c2 : experiments.getConditions()){
			if(c!=c2){
				head = head +"\t"+c.getName()+"_vs_"+c2.getName()+"_log2CPM"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2Fold"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2P";
			}
		}
		for (int bt=0; bt < numBindingTypes[c.getIndex()]; bt++){
			for (int s=0; s< 2;s++){
				head = head+"\tPoint\tTau";
				if(config.isAddingSequences() && findMotif)
						head = head +"\tSequence\tMotifScore";
			}	
		}

		if(config.isAddingAnnotations())
			head = head +"\tnearestGene\tdistToGene";
		return head;
	}
	
	/**
	 * Returns a string suitable for use as the header of a table whose rows are 
	 * the output of BindingEvent.getConditionString().
	 * (i.e. single condition files)
	 */
	public static String conditionSubtypeHeadString(ExperimentCondition c){
		String head="### ChExMix output\n";
		
		head = head + "#Condition\tName\tIndex\tTotalSigCount\tSignalFraction\n";
		head = head + "#Condition\t"+c.getName()+"\t"+c.getIndex()+"\t"+c.getTotalSignalCount()+"\t"+String.format("%.3f",c.getTotalSignalVsNoiseFrac())+"\n";
		
		head = head + "#Replicate\tParentCond\tName\tIndex\tSigCount\tCtrlCount\tCtrlScaling\tSignalFraction\n";
		for(ControlledExperiment r : c.getReplicates()){
			if(r.getControl()==null)
				head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t0\t1\t"+String.format("%.3f",r.getSignalVsNoiseFraction())+"\n";
			else
				head = head + "#Replicate\t"+c.getName()+"\t"+r.getName()+"\t"+r.getIndex()+"\t"+r.getSignal().getHitCount()+"\t"+r.getControl().getHitCount()+"\t"+String.format("%.3f",r.getControlScaling())+"\t"+String.format("%.3f",r.getSignalVsNoiseFraction())+"\n";
		}
		
		head = head + "#\n#Point";
		head = head +"\t"+c.getName()+"_Sig"+"\t"+c.getName()+"_Ctrl"+"\t"+c.getName()+"_log2Fold"+"\t"+c.getName()+"_log2P";
		
		for(ExperimentCondition c2 : experiments.getConditions()){
			if(c!=c2){
				head = head +"\t"+c.getName()+"_vs_"+c2.getName()+"_log2CPM"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2Fold"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2P";
			}
		}
		
		head = head+"\tSubtypePoint\tTau\tSubtypeName";

		if(config.isAddingAnnotations())
			head = head +"\tnearestGene\tdistToGene";
		return head;
	}
	
	/**
	 * Returns a string suitable for use as the header of a BED file whose rows are 
	 * the output of BindingEvent.getConditionBED().
	 * (i.e. single condition BED files)
	 */
	public static String conditionBEDHeadString(ExperimentCondition c){
		String head="track name=ChExMix-"+c.getName()+" description="+"ChExMix events "+c.getName();
		return head;
	}
	
	/**
	 * Returns a string suitable for use as the header of a BED file whose rows are 
	 * the output of BindingEvent.getConditionBED().
	 * (i.e. single condition BED files)
	 */
	public static String diffConditionBEDHeadString(ExperimentCondition c1, ExperimentCondition c2){
		String head="track name=ChExMix-"+c1.getName()+"-gt-"+c2.getName()+" description="+"ChExMix diff events "+c1.getName()+" > "+c2.getName();
		return head;
	}
	
	/**
	 * Returns a single-line string suitable for use as the header of a table of differential binding events
	 */
	public static String conditionShortHeadString(ExperimentCondition c){
		String head = "#Point";
		head = head +"\t"+c.getName()+"_Sig"+"\t"+c.getName()+"_Ctrl"+"\t"+c.getName()+"_log2Fold"+"\t"+c.getName()+"_log2P";
		
		for(ExperimentCondition c2 : experiments.getConditions()){
			if(c!=c2){
				head = head +"\t"+c.getName()+"_vs_"+c2.getName()+"_log2CPM"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2Fold"+"\t"+c.getName()+"_vs_"+c2.getName()+"_log2P";
			}
		}
		if(config.isAddingSequences());
			head = head +"\tSequence\tMotifScore";
		if(config.isAddingAnnotations())
			head = head +"\tnearestGene\tdistToGene";
		return head;
	}
	
	/**
	 * Print info on a single condition and associated inter-condition tests
	 */
	public String getConditionString(ExperimentCondition c) {
		//String gene = nearestGene == null ? "NONE" : nearestGene.getName();
		String out = point.getLocationString();	
		double logP = Math.log(getCondSigVCtrlP(c))/config.LOG2;
		double logF = Math.log(getCondSigVCtrlFold(c))/config.LOG2;
		//if(logP!=0.0)
		//	logP*=-1;
		out = out+"\t"+String.format("%.1f", getCondSigHits(c))+"\t"+String.format("%.1f", getCondCtrlHits(c))+"\t"+String.format("%.3f", logF)+"\t"+String.format("%.3f", logP); 
		
		for(ExperimentCondition c2 : experiments.getConditions()){
			if(c!=c2){
				double logPC = Math.log(getInterCondP(c, c2))/config.LOG2; 
				out = out+"\t"+String.format("%.2f", getInterCondScMean(c, c2))+"\t"+String.format("%.3f", getInterCondFold(c, c2))+"\t"+String.format("%.3f",logPC);
			}
		}	
		if (this.getTypePoints(c) != null){
			for (int bt=0; bt < this.getTypePoints(c).length;bt++){
				for (int s=0; s< 2;s++){
					out = out+"\t"+this.getTypePoints(c)[bt][s]+"\t"+String.format("%.2f", this.getTypeProbs(c)[bt][s]);
					if(config.isAddingSequences() && findMotif)
						if (this.getSequence(c)!=null)
							out = out+"\t"+this.getSequence(c)[bt][s]+"\t"+String.format("%.2f", this.getMotifScore(c)[bt][s]);
				}}}		
		
		if(config.isAddingAnnotations())
			if(nearestGene==null)
				out = out + "\tNONE\t"+config.getMaxAnnotDistance();
			else
				out = out + "\t"+nearestGene.getName()+"\t"+distToGene;
		//out = out+"\n";
		return out;
	}
	
	/**
	 * Print info on ML subtypes
	 */
	public String getSubtypeString(ExperimentCondition c){
		String out = point.getLocationString();	
		double logP = Math.log(getCondSigVCtrlP(c))/config.LOG2;
		double logF = Math.log(getCondSigVCtrlFold(c))/config.LOG2;
		//if(logP!=0.0)
		//	logP*=-1;
		out = out+"\t"+String.format("%.1f", getCondSigHits(c))+"\t"+String.format("%.1f", getCondCtrlHits(c))+"\t"+String.format("%.3f", logF)+"\t"+String.format("%.3f", logP); 
		double maxTau=0; int maxTypeIndex=0; int maxStrand=0;
		if (this.getTypeProbs(c)!= null){
			for (int bt=0; bt < this.getTypePoints(c).length;bt++){
				for (int s=0; s< 2;s++){
					if (this.getTypeProbs(c)[bt][s]> maxTau){
						maxTau=this.getTypeProbs(c)[bt][s];
						maxTypeIndex=bt; maxStrand=s;
					}
				}
			}
			out = out+"\t"+this.getTypePoints(c)[maxTypeIndex][maxStrand]+"\t"+String.format("%.2f", this.getTypeProbs(c)[maxTypeIndex][maxStrand])+"\tSubtype"+maxTypeIndex;
		}
		if(config.isAddingSequences() && findMotif)
			if (this.getSequence(c)!=null)
				out = out+"\t"+this.getSequence(c)[maxTypeIndex][maxStrand]+"\t"+String.format("%.2f", this.getMotifScore(c)[maxTypeIndex][maxStrand]);
		
		return out;
	}
	
	public Pair<Integer,StrandedPoint> getMaxSubtypePoint(ExperimentCondition c){
		double maxTau=0; int maxTypeIndex=0; int maxStrand=0;
		StrandedPoint p = null;
		if (this.getTypeProbs(c)!= null){
			for (int bt=0; bt < this.getTypePoints(c).length;bt++){
				for (int s=0; s< 2;s++){
					if (this.getTypeProbs(c)[bt][s]> maxTau){
						maxTau=this.getTypeProbs(c)[bt][s];
						maxTypeIndex=bt; maxStrand=s;
					}}}		
			p=this.getTypePoints(c)[maxTypeIndex][maxStrand];
		}
		return new Pair(maxTypeIndex,p);
	}
	
	
	/**
	 * Output in BED format
	 */
	public String getConditionBED(ExperimentCondition c){
		if (point != null) {
			double logP = Math.log(getCondSigVCtrlP(c))/config.LOG2;
			//BED score ranges from 0 to 1000
			//Arbitrary conversion:
			// score = 10 * -logP
			int score = (int)(Math.max(0, Math.min(-10*logP, 1000)));
			
			//Note that BED is 0-start and end-noninclusive, so we always shift the 
			// multiGPS binding event location back by one. 
			return new String("chr"+point.getChrom()+"\t"+(point.getLocation()-1)+"\t"+point.getLocation()+"\tChExMix\t"+score);
		}
		return null;
	}
	
	/**
	 * Returns the sequence window (length=win) centered on the point
	 */
	public String toSequence(int win){
		String seq;
		SequenceGenerator seqgen = new SequenceGenerator();
		Region peakWin=null;
		
		int start = point.getLocation()-((int)(win/2));
		if(start<1){start=1;}
		int end = point.getLocation()+((int)win/2)-1;
		if(end>point.getGenome().getChromLength(point.getChrom())){end =point.getGenome().getChromLength(point.getChrom());} 
		peakWin = new Region(point.getGenome(), point.getChrom(), start, end);
	
		String seqName = new String(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth()+"\t"+point.getLocation());
		String sq = seqgen.execute(peakWin);
		seq = seqName +"\n"+ sq+"\n";
		return seq;
	}
	
	/**
	 * Find the closest genes to the enriched features
	 * @param enriched
	 */
	public void addClosestGenes(){
		distToGene = config.getMaxAnnotDistance();
		for(AnnotationLoader loader : config.getGeneAnnotations()){
		    if (config.getAnnotOverlapOnly()) {
                for(Gene gene : loader.getGenes(point)){
                    if(gene.contains(point)){
                        nearestGene = gene;
                        distToGene = gene.distance(point);
                    }
                }
                if (nearestGene != null) {
                    distToGene = nearestGene.distance(point);
                }
            } else {
            	for(Gene gene : loader.getGenes(point)){
                    int distance = point.getLocation() - gene.getFivePrime();
                    if (gene.getStrand()=='-')
                    	distance = -distance;
                    if (Math.abs(distance) < Math.abs(distToGene)) {
                        nearestGene = gene;
                        distToGene = distance;                            
                    }                            
                }
            }
		}
	}
	
	public static void setExperimentManager(ExperimentManager e){experiments = e; sortingCondition = experiments.getConditions().get(0);}
	public static void setConfig(EventsConfig c){config = c;}
	public static void setGPSConfig(ChExMixConfig c){ gpsconfig=c;}
	public static void setNumBindingTypes(int[] bt){ numBindingTypes = bt;}
	public static void setSortingCond(ExperimentCondition c){sortingCondition = c;}
	public static void setSortingCondB(ExperimentCondition c){sortingConditionB = c;}
}
