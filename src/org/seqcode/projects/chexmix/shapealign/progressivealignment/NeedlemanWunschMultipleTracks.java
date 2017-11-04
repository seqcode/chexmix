package org.seqcode.projects.chexmix.shapealign.progressivealignment;

import java.util.Stack;

/**
 * NeedlemanWunschMultipleTracks: A global alignment algorithm to align two sets of multi-track profiles.
 * 
 * @author naomi yamada
 *
 */
public class NeedlemanWunschMultipleTracks {
	protected double[][] watsonA;	// watson tags from region A is indexed by genomic positions and experiments
	protected double[][] crickA;	// crick tags from region A is indexed by genomic positions and experiments
	protected double[][] watsonB;	// watson tags from region B is indexed by genomic positions and experiments
	protected double[][] crickB;	// crick tags from region B is indexed by genomic positions and experiments
	protected int pwindow;
	protected int numExpt;			// number of experiments
	protected Stack<Integer> traceBackTable = new Stack<Integer>();
	protected double maxScore;
	protected int alignStartXCoord = 0;
	protected int alignStartYCoord = 0;
	protected int alignEndXCoord = 0;
	protected int alignEndYCoord = 0;	
	static final double GAP_OPEN = 10;	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;	
	static final double MINIMUM_VALUE = -10000;
	
	public NeedlemanWunschMultipleTracks(double[][] watsonA, double[][] crickA, double[][] watsonB, double[][] crickB, int profilewidth, int numExpt){
		this.watsonA= watsonA;
		this.crickA = crickA;
		this.watsonB = watsonB;
		this.crickB = crickB;
		pwindow=profilewidth;
		this.numExpt = numExpt;
	}
	
	// accessors
	public double getMaxScore(){return maxScore;}
	public Stack<Integer> getTraceBack(){return traceBackTable;}
	public int getStartX(){return alignStartXCoord;}
	public int getStartY(){return alignStartYCoord;}
	public int getEndX(){return alignEndXCoord;}
	public int getEndY(){return alignEndYCoord;}

	public void buildMatrix(){			
		//initialization of M matrix
		double [][] M = new double [pwindow+1][pwindow+1];
		for (int i = 0 ; i <= pwindow; i++)
			for (int j = 0 ; j <= pwindow; j++)
				M[i][j] = 0;
							
		for (int i = 1 ; i <= pwindow; i++){
			for (int j = 1 ; j <= pwindow ; j++){				
				double mScore = computeEclideanScore(watsonA[i-1],watsonB[j-1],crickA[i-1],crickB[j-1]);			
				double temp_M[] = new double[3];
				temp_M[0] = M[i-1][j-1] + mScore;
				temp_M[1] = M[i][j-1] - GAP_OPEN;
				temp_M[2] = M[i-1][j] - GAP_OPEN;		
				double max_I = temp_M[0];
				for (int k = 1 ; k < 3 ; k++){
					if (temp_M[k] > max_I){ max_I = temp_M[k];}
				}
				M[i][j] = max_I;
			}
		}		
//		System.out.println("printing M matrix");
//		for (int i = 0; i <=window ; i++){
//			for (int j = 0; j <=window; j++){
//				System.out.print(M[i][j]+" ");
//			}
//			System.out.println();
//		}	
		
		// initialize
		double maximumScore = MINIMUM_VALUE;
		// find the highest value
		for (int i = (int) Math.floor(pwindow/2); i <= pwindow; i++){
			if (M[i][pwindow] > maximumScore){
				maximumScore = M[i][pwindow];
				alignEndXCoord = i;
				alignEndYCoord = pwindow;	
			}
		}
		for (int j = (int) Math.floor(pwindow/2); j <= pwindow; j++){
			if (M[pwindow][j] > maximumScore){
				maximumScore = M[pwindow][j];
				alignEndXCoord = pwindow;
				alignEndYCoord = j;	
			}
		}		
//		System.out.println("max score is "+maxScore+" x coord is "+alignEndXCoord+" y coord is "+alignEndYCoord);
		
		maxScore = maximumScore;
		
		// back track to reconstruct the path
		double currentScore = maxScore;
		alignStartXCoord = alignEndXCoord;
		alignStartYCoord = alignEndYCoord;		
		int i = alignEndXCoord;
		int j = alignEndYCoord;
		
		while ( i != 0 && j != 0){			
			double mScore = computeEclideanScore(watsonA[i-1],watsonB[j-1],crickA[i-1],crickB[j-1]);			
			if ( M[i-1][j-1] + mScore == currentScore ){	// diagonal case
				traceBackTable.push(DIAG);		
				currentScore = M[i-1][j-1];
				i--;
				j--;			
			}else if( M[i][j-1]-GAP_OPEN == currentScore ){	// left case
				traceBackTable.push(LEFT);
				currentScore = M[i][j-1];
				j--;			
			}else{	// right case
				System.out.println("current score is "+currentScore);
				System.out.println("M value is"+M[i-1][j]);
				traceBackTable.push(UP);
				currentScore = M[i-1][j];
				i--;
			}			
			alignStartXCoord = i;
			alignStartYCoord = j;
		}
	}
	
	/**
	 * Compute score by eclidean distance
	 * @param x1
	 * @param x2
	 * @param y1
	 * @param y2
	 * @return eclidean score
	 */
	public double computeEclideanScore(double[] x1, double[] x2, double[] y1, double[] y2){
		double squareDist = 0;
		for (int e=0; e< numExpt; e++){
			squareDist += Math.pow(x1[e]-x2[e], 2);
			squareDist += Math.pow(y1[e]-y2[e], 2);
		}
		double penalty = 0;
		for (int e=0; e<numExpt; e++){
			penalty += Math.abs(x1[e]-x2[e]);
			penalty += Math.abs(y1[e]-y2[e]);
		}
		return (1- Math.sqrt(squareDist) - penalty);
	}
	
	/**
	 * Compute score by chi square distance
	 * @param x1
	 * @param x2
	 * @param y1
	 * @param y2
	 * @return chi square score
	 */
	public double computeChiSquareScore(double[] x1, double[] x2, double[] y1, double[] y2){
		double chisquare = 0;
		for (int e=0; e<numExpt; e++){
			if (x1[e]!=0 || x2[e] !=0)
				chisquare += Math.pow(x1[e]-x2[e], 2)/(x1[e]+x2[e]);
			if (y1[e]!=0 || y2[e] !=0)
				chisquare += Math.pow(y1[e]-y2[e], 2)/(y1[e]+y2[e]);
		}
		double penalty = 0;
		for (int e=0; e<numExpt; e++){
			penalty += Math.abs(x1[e]-x2[e]);
			penalty += Math.abs(y1[e]-y2[e]);
		}
		return (1 - chisquare - penalty);
	}	
}