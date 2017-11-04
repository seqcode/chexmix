package org.seqcode.projects.chexmix.shapealign.alignment;

import java.util.Stack;

import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;

/**
 * SmithWatermanAlignmentMatrix: builds M and I matrix using SmithWaterman algorithms
 * 
 * @author naomi yamada
 *
 */

public class SmithWatermanAlignmentMatrix {
	protected ShapeAlignConfig shapeAlignConfig;
	protected SimilarityScore similarity_s;	
	protected double [][] regACounts;
	protected double [][] regBCounts;
	protected int window;	
	protected Stack<Integer> traceBackTable = new Stack<Integer>();
	protected double maxScore = MINIMUM_VALUE;
	protected int alignStartXCoord = 0;
	protected int alignStartYCoord = 0;
	protected int alignEndXCoord = 0;
	protected int alignEndYCoord = 0;	
	static final double GAP_OPEN = 10000;	
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;	
	static final double MINIMUM_VALUE = -10000;
	
	public SmithWatermanAlignmentMatrix (double[][] regAarray, double[][] regBarray, ShapeAlignConfig shapeConfig){		
		regACounts = regAarray;
		regBCounts = regBarray;
		window = regAarray.length;	
		shapeAlignConfig = shapeConfig;
		similarity_s = new SimilarityScore(shapeConfig);
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
		double [][] M = new double [window+1][window+1];
		for (int i = 0 ; i <= window; i++){
			for (int j = 0 ; j <= window; j++)
				M[i][j] = 0;
		}					
		
		for (int i = 1 ; i <= window; i++){
			for (int j = 1 ; j <= window ; j++){				
				double mScore = similarity_s.computeScore(regACounts[i-1][0],regBCounts[j-1][0],regACounts[i-1][1],regBCounts[j-1][1]);			
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

		// find the highest value
		for (int i = (int) Math.floor(window/2); i <= window; i++){
			if (M[i][window] > maxScore){
				maxScore = M[i][window];
				alignEndXCoord = i;
				alignEndYCoord = window;	
			}
		}
		for (int j = (int) Math.floor(window/2); j <= window; j++){
			if (M[window][j] > maxScore){
				maxScore = M[window][j];
				alignEndXCoord = window;
				alignEndYCoord = j;	
			}
		}		
//		System.out.println("max score is "+maxScore+" x coord is "+alignEndXCoord+" y coord is "+alignEndYCoord);
		
		// back track to reconstruct the path
		double currentScore = maxScore;
		alignStartXCoord = alignEndXCoord;
		alignStartYCoord = alignEndYCoord;		
		int i = alignEndXCoord;
		int j = alignEndYCoord;
		
		while ( i != 0 && j != 0){			
			double mScore = similarity_s.computeScore(regACounts[i-1][0],regBCounts[j-1][0],regACounts[i-1][1],regBCounts[j-1][1]);			
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
}