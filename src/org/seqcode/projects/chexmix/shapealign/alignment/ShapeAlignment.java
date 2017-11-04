package org.seqcode.projects.chexmix.shapealign.alignment;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.stats.DistributionPercentile;
import org.seqcode.projects.galaxyexo.FeatureCountsLoader;


/**
 * ShapeAlignment: It holds two methods to perform ChIP-exo footprint alignment. 
 * - Pearson correlation
 * - Smith-Waterman overlap alignment with various similarity metrics
 * 
 * @author naomi yamada
 *
 */

public class ShapeAlignment {	
	protected GenomeConfig gconf;
	protected ExperimentManager manager;		
	protected ShapeAlignConfig shapeAlignConfig;
	protected List<StrandedRegion> strandedRegions;
	
	protected int window;	
	protected double error = 0;
	protected double totalNum = 0;	
	static final double MINIMUM_VALUE = -10000;
	
	protected Map<ControlledExperiment, Map<StrandedRegion, double[][]>> strandedRegionSampleCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>(); 	
	protected double [] offsetArray;
	protected double [][] pairwiseSimilarities; // This includes pairwise correlations
	protected double [][] distanceMatrix;
	protected Map<StrandedRegion, double[]> regionDistances;
	
	/*** Smith-Waterman specific parameters  ***/
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;	

	public ShapeAlignment(ShapeAlignConfig shapeAlignmentConf, GenomeConfig gcon,ExperimentManager man){	
		shapeAlignConfig = shapeAlignmentConf;
		gconf = gcon;
		manager = man;	
		window = shapeAlignmentConf.getWindowSize();
	}
	
	// getters
	public List<StrandedRegion> getStrandedRegions(){return strandedRegions;}
	public double[][] getSimilarityMatrix(){return pairwiseSimilarities;}
	
	public double[][] getDistanceMatrix(){
		distanceMatrix = new double[pairwiseSimilarities[0].length][pairwiseSimilarities.length];
		double maxVal = -1;
		for (int i = 0 ; i < pairwiseSimilarities.length; i++){
			for (int j = 0 ; j < pairwiseSimilarities[i].length; j++){
				if (pairwiseSimilarities[j][i] > maxVal){
					maxVal = pairwiseSimilarities[j][i];
				}
			}
		}
		for (int i = 0 ; i < pairwiseSimilarities.length; i++){
			for (int j = 0 ; j < pairwiseSimilarities[i].length; j++){
				distanceMatrix[j][i] = maxVal-pairwiseSimilarities[j][i];
			}
		}
		return distanceMatrix;
	}
	
	// return the matrix threshold by a percentile
	public double[][] thresholdDistance(){
		double[][] distanceM = getDistanceMatrix();		
		double[][] tDistances = new double[distanceM[0].length][distanceM.length];
		DistributionPercentile distr = new DistributionPercentile(distanceM);
		double half = distr.getPercentileValue(50);
		double TwentyFive = distr.getPercentileValue(25);
		double Fifteen = distr.getPercentileValue(15);
		double perVal = distr.getPercentileValue(shapeAlignConfig.getPercentile());
		double five = distr.getPercentileValue(5);
		
		System.out.println("within shapeAlignment class");
		System.out.println("max value is "+distr.getMaxVal());
		System.out.println("50% value is "+half);
		System.out.println("25% value is "+TwentyFive);
		System.out.println("15% value is "+Fifteen);
		System.out.println("10% value is "+perVal);
		System.out.println("5% value is "+five);
		
		for (int i = 0 ; i <distanceM.length; i ++){
			for (int j = 0 ; j < distanceM[i].length; j++){
				if (distanceM[j][i]< perVal){
					tDistances[j][i] = distanceM[j][i];
				}else{
					tDistances[j][i] = distr.getMaxVal();
				}
			}
		}
		return tDistances;
	}
	
	
	/**printing methods  **/
	public void printErrorRate(){System.out.println("error is "+error+ " totalNum is "+ totalNum+ " error rate is "+ (error/totalNum));}	
	public void printRegions() throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(shapeAlignConfig.getOutDir()+".regions","UTF-8");
		for (int i = 0 ; i < strandedRegions.size(); i++){
			writer.println(strandedRegions.get(i).getLocationString());
		}	
		writer.close();
	}
	public void printOffsetArray(){
		System.out.println("offset array");
		for (int i = 0; i <= window ; i++)
			System.out.print(offsetArray[i]+"\t");
		System.out.println();
	}
	// print similarity matrix
	public void printSimScores(){	
		for (int i = 0 ; i < strandedRegions.size(); i++){
			System.out.println(strandedRegions.get(i).getLocationString());
		}		
		for (int i = 0 ; i < strandedRegions.size();i++){
			for (int j = 0 ; j <strandedRegions.size(); j++){
				System.out.print(pairwiseSimilarities[i][j]+"\t");
			}
			System.out.println();	
		}
	}	
	// print distance matrix threshold by percentile distributions
	public void printThresDistances() throws FileNotFoundException, UnsupportedEncodingException{			
		PrintWriter writer = new PrintWriter(shapeAlignConfig.getOutDir()+".distances","UTF-8");	
		double[][] thresDistances = thresholdDistance();
		for (int i = 0 ; i < thresDistances.length;i++){
			for (int j = 0 ; j <thresDistances[0].length; j++){
				writer.print(thresDistances[i][j]+"\t");
			}
			writer.println();
		}
		writer.close();
	}
	
	public void excute() throws FileNotFoundException, UnsupportedEncodingException{		
		FeatureCountsLoader fcLoader = null;
		for (ExperimentCondition condition : manager.getConditions()){		
			for (ControlledExperiment rep: condition.getReplicates()){	
				if (!shapeAlignConfig.isSmithWaterman()){ // Pearson Correlation case
					fcLoader = new FeatureCountsLoader(gconf, shapeAlignConfig.getStrandedPoints(), window*2);
				}else{
					fcLoader = new FeatureCountsLoader(gconf, shapeAlignConfig.getStrandedPoints(), window);
				}
				strandedRegionSampleCounts.put(rep, fcLoader.strandedRegionSampleCounts(rep));
			}
		}
		strandedRegions = fcLoader.getStrandedRegions();
//		shapeAlignConfig.setStrandedRegions(strandedRegions);
		
		//initialize offsetArray
		offsetArray = new double [window+1];
		for (int i = 0 ; i <= window ; i++)
			offsetArray[i] = 0;	
		
		// initialize pairwiseDistance
		pairwiseSimilarities = new double [strandedRegions.size()][strandedRegions.size()];
		if (!shapeAlignConfig.isSmithWaterman()){ // Pearson Correlation case
			for (int i = 0 ; i< strandedRegions.size(); i++){
				for (int j = 0 ; j < strandedRegions.size(); j++)
					pairwiseSimilarities[i][j] = 1;		
			}
			// loop through all possible pairwise region to do Pearson correlation
			for (ControlledExperiment cexpt : strandedRegionSampleCounts.keySet()){
				for (int i = 0; i <strandedRegions.size();i++){	
					for (int j = i+1; j <strandedRegions.size();j++){
						double corr = pearsonCorrelation(cexpt, strandedRegions.get(i), strandedRegions.get(j));
						pairwiseSimilarities[i][j] = corr;
						pairwiseSimilarities[j][i] = corr;
					}
				}
			}
		}else{ // Smith-Waterman case
			for (int i = 0 ; i< strandedRegions.size(); i++){
				for (int j = 0 ; j < strandedRegions.size(); j++)
					pairwiseSimilarities[i][j] = 100;
			}
			for (ControlledExperiment cExpt : strandedRegionSampleCounts.keySet()){
				for (int i = 0; i <strandedRegions.size();i++){						
					for (int j = i+1; j <strandedRegions.size();j++){
						double maxScore = smithWatermanAlgorithm(cExpt, strandedRegions.get(i), strandedRegions.get(j));	
						pairwiseSimilarities[i][j] = maxScore;
						pairwiseSimilarities[j][i] = maxScore;
					}
				}
			}	
		}	
		if(shapeAlignConfig.printError() ){this.printErrorRate();}
		if(shapeAlignConfig.printOffsetArray()){this.printOffsetArray();}
		if(shapeAlignConfig.printThresDistances()){			
			this.printRegions();
			this.printThresDistances();
			}
		
	}	
	
	public double pearsonCorrelation(ControlledExperiment controlledExpt, Region regA, Region regB){
		//get counts
		double [][] regACounts = strandedRegionSampleCounts.get(controlledExpt).get(regA);
		double [][] regBCounts = strandedRegionSampleCounts.get(controlledExpt).get(regB);		
		//normalize the regAarray to set the max value 1
		// size of normRegACounts is [window_size][2]
		double maxA = MINIMUM_VALUE;
		for (int i = 0 ; i <=window ; i++){
			for (int s = 0 ; s <2 ; s++){
				if (regACounts[i+window/2][s] > maxA){maxA = regACounts[i+window/2][s];}
			}
		}
		double [][] normRegACounts = new double[window+1][2];
		for (int i = 0; i <=window; i++){
			for (int s = 0 ; s < 2 ; s++)
				normRegACounts[i][s] = regACounts[i+window/2][s]/maxA;
		}		
		// compute Pearson correlation with sliding window
		double maxCorr = -1;
		boolean reverse = false;
		int offset = 0; 		
		for (int slidingW = 0; slidingW <=window; slidingW++){			
			// copy and normalize B array
			double maxB = MINIMUM_VALUE;
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					if (regBCounts[slidingW+i][s] > maxB){maxB = regBCounts[slidingW+i][s];}
				}			
			}
			double [][] normRegBCounts = new double [window+1][2];
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					normRegBCounts[i][s] = regBCounts[slidingW+i][s]/maxB;
				}
			}			
			//calculate Pearson Correlation
			PearsonCorrelation pc = new PearsonCorrelation(normRegACounts,normRegBCounts);			
			double r = -1;
			boolean fromReverse = false;
			if (!shapeAlignConfig.weightedPC()){
				r = pc.computePearsonCorr();
				fromReverse = pc.isReverse();
			}else{
				r = pc.computeWeightedPearsonCorr();
				fromReverse = pc.isReverse();				
			}	
			if (r > maxCorr ){
				maxCorr = r;
				reverse = fromReverse;
				offset = slidingW - window/2;
			}					
		}
//		System.out.println("offsets is "+offset);
//		System.out.println("correlation is "+maxCorr);

		double[][] alignedRegB = new double[window+1][2];
		double max_b = MINIMUM_VALUE;
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				if (regBCounts[window/2+offset+i][s] > max_b){max_b = regBCounts[window/2+offset+i][s];}
			}
		}		
//		System.out.println("max_b is "+max_b);

		for (int i = 0; i <= window; i++){
			for (int s = 0; s <2 ; s++){
				alignedRegB[i][s] = regBCounts[window/2+offset+i][s]/max_b;
			}
		}
		
		double temp_max_b = MINIMUM_VALUE;
		double [][] norm_b = new double[window+1][2];
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				if (regBCounts[window/2+i][s] > temp_max_b){temp_max_b = regBCounts[window/2+i][s];}
			}
		}
		for (int i = 0; i <=window; i++){
			for (int s = 0; s <2 ; s++){
				norm_b[i][s] = regBCounts[window/2+i][s]/temp_max_b;
			}
		}		

		// increment offset array
		offsetArray[offset+window/2]++;				
		// incrementing error allowing offset of +-1
		totalNum += 1;
		if (Math.abs(offset) > 1){
			error += 1;
//			System.out.println(regA.getLocationString());
//			System.out.println(regB.getLocationString());
			}		
		return maxCorr;	
	}
	
	public double smithWatermanAlgorithm(ControlledExperiment cExpt, Region regA, Region regB){		
		//get counts
		double [][] regACounts = strandedRegionSampleCounts.get(cExpt).get(regA);
		double [][] regBCounts = strandedRegionSampleCounts.get(cExpt).get(regB);		
		//normalize the arrays to set the max value 1
		double maxA = MINIMUM_VALUE;
		double maxB = MINIMUM_VALUE;
		for (int i = 0; i <=window ; i++){
			for (int s = 0 ; s < 2 ; s++){
				if (regACounts[i][s] > maxA){maxA = regACounts[i][s];}
				if (regBCounts[i][s] > maxB){maxB = regBCounts[i][s];}
			}
		}		
//		System.out.println("max counts are "+maxA+" : "+maxB);
		
		double [][] normRegACounts = new double [window+1][2];
		double [][] normRegBCounts = new double [window+1][2];
		double [][] normRegBRevCounts = new double[window+1][2];
		for (int i = 0 ; i <= window; i++){
			for (int s = 0 ; s < 2 ; s++){
				normRegACounts[i][s] = regACounts[i][s]/maxA;
				normRegBCounts[i][s] = regBCounts[i][s]/maxB;
			}
		}		
		//reversing normRegBCounts
		for (int i = 0; i <= window; i++){
			for (int s = 0; s < 2 ;s++){
				normRegBRevCounts[window-i][1-s] = normRegBCounts[i][s];
			}
		}
		// align using two possible ways
		SmithWatermanAlignmentMatrix alignOne = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBCounts,shapeAlignConfig);
		SmithWatermanAlignmentMatrix alignTwo = new SmithWatermanAlignmentMatrix(normRegACounts,normRegBRevCounts,shapeAlignConfig);
		alignOne.buildMatrix();
		alignTwo.buildMatrix();
		
		Stack<Integer> traceBack = new Stack<Integer>();
		double[][] regBarray = new double[window+1][2];
		double maxScore = MINIMUM_VALUE;
		boolean reverseB;
		int s_x_coord = 0;
		int s_y_coord = 0;
		int e_x_coord = 0;
		int e_y_coord = 0;		
		if (alignOne.getMaxScore() > alignTwo.getMaxScore()){	
			regBarray = normRegBCounts;
			maxScore = alignOne.getMaxScore();
			traceBack = alignOne.getTraceBack();
			s_x_coord = alignOne.getStartX();
			s_y_coord = alignOne.getStartY();
			e_x_coord = alignOne.getEndX();
			e_y_coord = alignOne.getEndY();		
			reverseB = false;		
		}else{	
			maxScore = alignTwo.getMaxScore();
			regBarray = normRegBRevCounts;
			traceBack = alignTwo.getTraceBack();
			s_x_coord = alignTwo.getStartX();
			s_y_coord = alignTwo.getStartY();
			e_x_coord = alignTwo.getEndX();
			e_y_coord = alignTwo.getEndY();	
			reverseB = true;
		}		
		double x_mid = (s_x_coord + e_x_coord)/2;
		double y_mid = (s_y_coord + e_y_coord)/2;		
//		System.out.println("alignment start coordinates "+ s_x_coord + " : " + s_y_coord);
//		System.out.println("alignment end coordinates "+ e_x_coord + " : " + e_y_coord);		
//		System.out.println(Arrays.toString(traceBack.toArray()));
		
		double[][] alignedRegA = new double[e_x_coord-s_x_coord+1][2];
		double[][] alignedRegB = new double[e_y_coord-s_y_coord+1][2];
		int current_x = e_x_coord-1;
		int current_y = e_y_coord-1;
		
		// trace back for x 
		@SuppressWarnings("unchecked")
		Stack<Integer> xTraceBack = (Stack<Integer>) traceBack.clone();		
		for (int i = e_x_coord-s_x_coord; i >= 0 ; i--){				
			if (current_x >= 0){
				for (int s = 0 ; s <2; s++)
					alignedRegA[i][s] = normRegACounts[current_x][s];			
				if ( !xTraceBack.empty() ){			
					if (xTraceBack.peek() == DIAG || xTraceBack.peek() == LEFT){
						current_x --;
					}			
					xTraceBack.pop();
				}	
			}
		}				
		// trace back for y 
		@SuppressWarnings("unchecked")
		Stack<Integer> yTraceBack = (Stack<Integer>) traceBack.clone();;
		
		for (int i = e_y_coord-s_y_coord; i >= 0 ; i--){				
			if (current_y >= 0){			
				for (int s = 0 ; s <2; s++){
					if (reverseB == false){
						alignedRegB[i][s] = normRegBCounts[current_y][s];
					}else{
						alignedRegB[i][s] = normRegBRevCounts[current_y][s];
					}
				}
				if ( !yTraceBack.empty() ){				
					if (yTraceBack.peek() == DIAG || yTraceBack.peek() == UP){
						current_y --;	
					}
					yTraceBack.pop();			
				}
			}
		}	

		// increment offset array
		if (y_mid-x_mid >= - window/2 && y_mid-x_mid <=  window/2)
			offsetArray[(int) (y_mid-x_mid+window/2)]++;				
		// incrementing error allowing offset of +-1
		totalNum += 1;
		if ( traceBack.contains(LEFT) || traceBack.contains(UP) ){ // check that stack only contains DIAG
			error += 1;
		}else{
			if (Math.abs(y_mid-x_mid) >1){
				error +=1;
			}
		}
		return maxScore;
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException{
		
		ShapeAlignConfig shapeAlignConf = new ShapeAlignConfig(args);	
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);	
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);

		ShapeAlignment profile = new ShapeAlignment(shapeAlignConf, gconf, manager); 	
		profile.excute();
		profile.printErrorRate();
		
		manager.close();		
	}
}
