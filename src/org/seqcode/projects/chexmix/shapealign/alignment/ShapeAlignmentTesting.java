package org.seqcode.projects.chexmix.shapealign.alignment;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
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
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.math.probability.NormalDistribution;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.stats.DistributionPercentile;
import org.seqcode.projects.galaxyexo.FeatureCountsLoader;


/**
 * ShapeAlignmentTesting: 
 * 		Edited copy of ShapeAlignment so Shaun can test some clustering approaches. 
 * 
 * : It holds two methods to perform ChIP-exo footprint alignment. 
 * - Pearson correlation
 * - Smith-Waterman overlap alignment with various similarity metrics
 * 
 * @author naomi yamada
 *
 */

public class ShapeAlignmentTesting {	
	protected GenomeConfig gconf;
	protected ExperimentManager manager;		
	protected ShapeAlignConfig shapeAlignConfig;
	protected List<StrandedRegion> strandedRegions;
	
	protected int window;
	protected int slidingWinOffset=25;
	protected double error = 0;
	protected double totalNum = 0;	
	static final double MINIMUM_VALUE = -10000;
	protected double smoothingGaussianVariance =1;
	
	protected Map<ControlledExperiment, Map<StrandedRegion, double[][]>> strandedRegionSampleCounts = new HashMap<ControlledExperiment, Map<StrandedRegion,double[][]>>(); 	
	protected double [] offsetArray;
	protected double [][] pairwiseSimilarities; // This includes pairwise scores
	protected int [][] pairwiseOffsets; // This includes pairwise offsets
	protected boolean [][] pairwiseStrands; // This includes pairwise strands
	protected double [][] distanceMatrix;
	protected Map<StrandedRegion, double[]> regionDistances;
	protected double[] gaussian;
	
	protected boolean pearsonSim=false;
	protected boolean euclideanSim=true;
	protected boolean kldivSim=false;
	protected boolean normScores=true;
	//Temporary params
	private boolean currAlignIsReverse;
	private int currAlignOffset;
	
	/*** Smith-Waterman specific parameters  ***/
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;	

	public ShapeAlignmentTesting(ShapeAlignConfig shapeAlignmentConf, GenomeConfig gcon,ExperimentManager man){	
		shapeAlignConfig = shapeAlignmentConf;
		gconf = gcon;
		manager = man;	
		window = shapeAlignmentConf.getWindowSize();
		initGaussianKernel(smoothingGaussianVariance);
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
	
	//pre-calculate and store the Guassian kernel prob., for efficiency
	private void initGaussianKernel(double var){
		gaussian = new double[250]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, var*var);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
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
			System.out.print("\t"+strandedRegions.get(i).getLocationString());
		}System.out.print("\n");
		for (int i = 0 ; i < strandedRegions.size();i++){
			System.out.print(strandedRegions.get(i).getLocationString());
			for (int j = 0 ; j <strandedRegions.size(); j++){
				System.out.print("\t"+pairwiseSimilarities[i][j]);
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
	
	public void execute() throws FileNotFoundException, UnsupportedEncodingException{		
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
		
		//initialize offsetArray
		offsetArray = new double [window+1];
		for (int i = 0 ; i <= window ; i++)
			offsetArray[i] = 0;	
		
		// initialize pairwiseDistance
		pairwiseSimilarities = new double [strandedRegions.size()][strandedRegions.size()];
		pairwiseStrands = new boolean [strandedRegions.size()][strandedRegions.size()];
		pairwiseOffsets = new int [strandedRegions.size()][strandedRegions.size()];
		if (!shapeAlignConfig.isSmithWaterman()){ // Pearson Correlation case
			double maxScore=-Double.MAX_VALUE;
			double minScore=Double.MAX_VALUE;
			for (int i = 0 ; i< strandedRegions.size(); i++){
				for (int j = 0 ; j < strandedRegions.size(); j++){
					pairwiseSimilarities[i][j] = 1;
					pairwiseStrands[i][j] = true;		
					pairwiseOffsets[i][j] = 0;		
				}
			}
			// loop through all possible pairwise regions
			for (ControlledExperiment cexpt : strandedRegionSampleCounts.keySet()){
				for (int i = 0; i <strandedRegions.size();i++){	
					for (int j = 0; j <strandedRegions.size();j++){
						double score = slidingComparison(cexpt, strandedRegions.get(i), strandedRegions.get(j));
						pairwiseSimilarities[i][j] = score;
						pairwiseStrands[i][j] = !currAlignIsReverse;		
						pairwiseOffsets[i][j] = currAlignOffset;
						if(score>maxScore)
							maxScore=score;
						else if(score<minScore)
							minScore=score;
					}
				}
			}
			
			if(normScores){
				for (int i = 0; i <strandedRegions.size();i++){	
					for (int j = 0; j <strandedRegions.size();j++){
						pairwiseSimilarities[i][j] = (pairwiseSimilarities[i][j]-minScore)/(maxScore-minScore);
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
		if(shapeAlignConfig.printSimScore()){			
			this.printSimScores();
		}
		
		//Print similarity score
		System.out.println("similarity score");
		for (int i=0; i < pairwiseSimilarities.length; i++){
			for (int j=0; j < pairwiseSimilarities[i].length ;j++)
				System.out.print(pairwiseSimilarities[i][j]+",");
			System.out.println();
		}
		System.out.println("strand");
		for (int i=0; i < pairwiseStrands.length; i++){
			for (int j=0; j < pairwiseStrands[i].length ;j++)
				System.out.print(pairwiseStrands[i][j]+",");
			System.out.println();
		}
		System.out.println("offset");
		for (int i=0; i < pairwiseOffsets.length; i++){
			for (int j=0; j < pairwiseOffsets[i].length ;j++)
				System.out.print(pairwiseOffsets[i][j]+",");
			System.out.println();
		}
		
	}	
	
	public double slidingComparison(ControlledExperiment controlledExpt, Region regA, Region regB){
		//get counts
		double [][] regACounts = strandedRegionSampleCounts.get(controlledExpt).get(regA);
		double [][] regBCounts = strandedRegionSampleCounts.get(controlledExpt).get(regB);		

		
		//Transpose matrices... why are they defined as array[window_size][2]!?
		double [][]regACountsT = new double[2][1+window*2];
		double [][]regBCountsT = new double[2][1+window*2];
		for (int i = 0 ; i <=window*2 ; i++){
			for(int s=0; s<2; s++){
				regACountsT[s][i] = regACounts[i][s];
				regBCountsT[s][i] = regBCounts[i][s];
			}
		}

		//Run Gaussian smoothing kernels over both arrays here
		double[][] densitiesA = new double[2][];
		double[][] densitiesB = new double[2][];
		densitiesA[0] = StatUtil.symmetricKernelSmoother(regACountsT[0], gaussian);
		densitiesA[1] = StatUtil.symmetricKernelSmoother(regACountsT[1], gaussian);
		densitiesB[0] = StatUtil.symmetricKernelSmoother(regBCountsT[0], gaussian);
		densitiesB[1] = StatUtil.symmetricKernelSmoother(regBCountsT[1], gaussian);
		
		
		//normalize the regArray to contain fractions of total reaad counts
		// size of normRegACounts is [window_size][2]		
		double sumA = 0.0;
		for (int i = 0 ; i <=window ; i++){
			for (int s = 0 ; s <2 ; s++)
				sumA+=densitiesA[s][i+window/2];
		}
		double [][] normRegACounts = new double[window+1][2];
		for (int i = 0; i <=window; i++){
			for (int s = 0 ; s < 2 ; s++)
				normRegACounts[i][s] = densitiesA[s][i+window/2]/sumA;
		}		

		
		// compute scores with sliding window
		double maxSim = -1;
		boolean reverse = false;
		int offset = 0; 		
		for (int slidingW =(window/2)-slidingWinOffset ; slidingW <=(window/2)+slidingWinOffset; slidingW++){			
			// copy and normalize B array
			//double maxB = MINIMUM_VALUE;
			double sumB=0.0;
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					sumB+=densitiesB[s][slidingW+i];
				}			
			}
			double [][] normRegBCounts = new double [window+1][2];
			for (int i = 0; i <=window; i++){
				for (int s = 0 ; s < 2 ; s++){
					normRegBCounts[i][s] = densitiesB[s][slidingW+i]/sumB;
				}
			}

			
			//calculate Pearson Correlation
			double sim =0;
			boolean fromReverse = false;
		
			if(pearsonSim==true){
				PearsonCorrelation pc = new PearsonCorrelation(normRegACounts,normRegBCounts);			
				if (!shapeAlignConfig.weightedPC()){
					sim = pc.computePearsonCorr();
					fromReverse = pc.isReverse();
				}else{
					sim = pc.computeWeightedPearsonCorr();
					fromReverse = pc.isReverse();				
				}	
			}else if(euclideanSim==true){
				EuclideanSqDistance e = new EuclideanSqDistance(normRegACounts,normRegBCounts);			
				sim = -1*e.computeDistance();
				fromReverse = e.isReverse();
			}else if(kldivSim==true){
				KLDivergence kl = new KLDivergence(normRegACounts,normRegBCounts);			
				sim = -1*kl.computeDistance();
				fromReverse = kl.isReverse();
			}
			
			
			
			if (sim > maxSim ){
				maxSim = sim;
				reverse = fromReverse;
				offset = slidingW - window/2;
			}					
		}
		
		//System.out.println("maxSim = "+maxSim);
		//System.out.println("offset = "+offset);
		//System.out.println((reverse ? "reverse" : "forward"));

		double[][] alignedRegB = new double[window+1][2];
		double sumB=0.0;
		for (int i = 0; i <=window; i++){
			for (int s = 0 ; s < 2 ; s++){
				sumB+=densitiesB[s][(window/2)+offset+i];
			}			
		}
		for (int i = 0; i <=window; i++){
			for (int s = 0 ; s < 2 ; s++){
				alignedRegB[i][s] = densitiesB[s][(window/2)+offset+i]/sumB;
			}
		}

		
		// increment offset array
		offsetArray[offset+window/2]++;				
		// incrementing error allowing offset of +-1
		totalNum += 1;
		if (Math.abs(offset) > 1){
			error += 1;
		}	
		
		currAlignOffset = offset;
		currAlignIsReverse=reverse;
		return maxSim;	
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
	
	//Generate Composites based on alignments against an exemplar
	public double[][] getExemplarBasedAlignment(ControlledExperiment controlledExpt, int exemplarIndex, List<Integer> alignList){
		double[][] composite = new double[window+1][2];
		double[] coverage = new double[window+1];
		for(int i=0; i<=window; i++){
			coverage[i]=0.0;
			for (int s = 0 ; s < 2 ; s++)
				composite[i][s]=0.0;
		}
		
		StrandedRegion exReg = strandedRegions.get(exemplarIndex);
		double [][] regACounts = strandedRegionSampleCounts.get(controlledExpt).get(exReg);
		for (int i = 0; i<=window; i++){
			for (int s = 0 ; s < 2 ; s++)
				composite[i][s] += regACounts[i+window/2][s];
			coverage[i]++;
		}	

		for(Integer x : alignList){
			if(x!=exemplarIndex){
				int currOffset = pairwiseOffsets[exemplarIndex][x];
				boolean currStrand = pairwiseStrands[exemplarIndex][x];
				
				StrandedRegion bReg = strandedRegions.get(x);
				double [][] regBCounts = strandedRegionSampleCounts.get(controlledExpt).get(bReg);
				for (int i = 0; i<=window; i++){
					for (int s = 0 ; s < 2 ; s++){
						if(currStrand)
							composite[i][s] += regBCounts[currOffset+i+window/2][s];
						else
							composite[window-i][1-s] += regBCounts[currOffset+i+window/2][s];
					}
					coverage[i]++;
				}
			}
		}
		
		for (int i = 0; i<=window; i++){
			for (int s = 0 ; s < 2 ; s++)
				composite[i][s] /= coverage[i];
		}
		
		return composite;
	}
	
	//Generate Composites based on alignments against an exemplar
	public List<StrandedPoint> getExemplarBasedAlignedPoints(int exemplarIndex, List<Integer> alignList){
			
		List<StrandedPoint> aligned = new ArrayList<StrandedPoint>();
		
		StrandedPoint exPt = shapeAlignConfig.getStrandedPoints().get(exemplarIndex);
		StrandedRegion exReg = strandedRegions.get(exemplarIndex);
		if(!exReg.contains(exPt)){System.err.println("ERROR: misalignment of points & regions");}
		aligned.add(exPt);
		
		for(Integer x : alignList){
			if(x!=exemplarIndex){
				int currOffset = pairwiseOffsets[exemplarIndex][x];
				boolean currStrand = pairwiseStrands[exemplarIndex][x];
				
				StrandedPoint bPt = shapeAlignConfig.getStrandedPoints().get(x);
				StrandedRegion bReg = strandedRegions.get(x);
				if(!bReg.contains(bPt)){System.err.println("ERROR: misalignment of points & regions");}
	
				char strand = bPt.getStrand();
				if(!currStrand)
					if(strand=='+')
						strand='-';
					else
						strand='+';
				
				StrandedPoint shiftBPt = new StrandedPoint(bPt.getGenome(), bPt.getChrom(), 
						bPt.getLocation()+currOffset, strand);
					
				aligned.add(shiftBPt);
			}
		}
		
		return aligned;
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException{
		
		ShapeAlignConfig shapeAlignConf = new ShapeAlignConfig(args);	
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);	
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);

		ShapeAlignmentTesting profile = new ShapeAlignmentTesting(shapeAlignConf, gconf, manager); 	
		profile.execute();
		//profile.printErrorRate();
		
		manager.close();		
	}
}
