package org.seqcode.projects.chexmix.shapealign.progressivealignment;

import java.util.Stack;

/**
 * NeedlemanWunsch: A global alignment algorithm to align two profiles
 * 
 * @author naomi yamada
 *
 */
public class NeedlemanWunsch {
	protected double[] watsonA;
	protected double[] crickA;
	protected double[] watsonB;
	protected double[] crickB;
	protected int pwindow;
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
	
	public NeedlemanWunsch(double[] watsonA, double[] crickA, double[] watsonB, double[] crickB, int profilewidth){
		this.watsonA= watsonA;
		this.crickA = crickA;
		this.watsonB = watsonB;
		this.crickB = crickB;
		pwindow=profilewidth;		
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
		for (int i = 0 ; i <= pwindow; i++){
			for (int j = 0 ; j <= pwindow; j++)
				M[i][j] = 0;
		}					
		
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
	 * @return score
	 */
	public double computeEclideanScore(double x1, double x2, double y1, double y2){
		double score = 1 - Math.sqrt(Math.pow(x1-x2, 2) + Math.pow(y1-y2, 2)) - Math.abs(x1-x2)-Math.abs(y1-y2);		
		return score;
	}
	
	/**
	 * Compute score by chi square distance
	 * @param x1
	 * @param x2
	 * @param y1
	 * @param y2
	 * @return score
	 */
	public double computeChiSquareScore(double x1, double x2, double y1, double y2){
		double score = 0;		
		if (x1 == x2 && x1 == 0){
			score = 1 - Math.pow(y1-y2, 2)/(y1+y2) - Math.abs(y1-y2);
		}else if (y1 == y2 && y2 == 0){
			score = 1 - Math.pow(x1-x2, 2)/(x1+x2) - Math.abs(x1-x2);
		}else{
			score =  1 - (Math.pow(x1-x2, 2)/(x1+x2) + Math.pow(y1-y2, 2)/(y1+y2)) - Math.abs(x1-x2)-Math.abs(y1-y2);
		}	
		return score;
	}
	
}