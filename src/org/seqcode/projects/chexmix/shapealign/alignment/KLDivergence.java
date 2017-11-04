package org.seqcode.projects.chexmix.shapealign.alignment;

import org.seqcode.math.stats.StatUtil;

/**
 * KLDivergence
 * 
 *
 */

public class KLDivergence {	
	protected double [][] normAarray;
	protected double [][] normBarray;
	protected double [][] revNormBarray;
	protected int window;	
	protected boolean revAlign = false; 
	
	public KLDivergence(double arrayA[][], double arrayB[][]){	
		normAarray = arrayA;
		normBarray = arrayB;
		window = arrayB.length-1;		
		revNormBarray = new double [window+1][2];		
		for (int i = 0; i <= window; i++){
			for (int s = 0; s < 2 ;s++){
				revNormBarray[window-i][1-s] = arrayB[i][s];
			}
		}	
	}
	
	// accessor 
	public boolean isReverse(){return revAlign;}
	
	public double computeDistance(){		
		double [] catAarray = new double [2*(window+1)]; 
		double [] catBarray = new double [2*(window+1)]; 
		double [] catRevBarray = new double [2*(window+1)]; 
		
		// concatenate forward and reverse strands to make 1 dimensional arrays
		for (int i = 0 ; i <= window; i++){			
			catAarray[i] = normAarray[i][0];
			catAarray[i+window+1] = normAarray[i][1];
			
			catBarray[i] = normBarray[i][0];
			catBarray[i+window+1] = normBarray[i][1];
			
			catRevBarray[i] = revNormBarray[i][0];
			catRevBarray[i+window+1] = revNormBarray[i][1];	
		}

				
		// Make sure that p and q are all proper probability values
		StatUtil.mutate_normalize(catAarray);
		StatUtil.mutate_normalize(catBarray);
		StatUtil.mutate_normalize(catRevBarray);
		
		// calculate squared distances from forward and reverse pair
		double distF=0.0, distR=0.0;
		for (int i = 0; i < catAarray.length; i++){
			distF += catAarray[i]*(Math.log(catAarray[i])-Math.log(catBarray[i]));
			distR += catAarray[i]*(Math.log(catAarray[i])-Math.log(catRevBarray[i]));
		}
		
		if (distF < distR){
			revAlign = false;
			return distF;
		}else{
			revAlign = true;
			return distR;
		}
	}
		
}
