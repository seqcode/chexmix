package org.seqcode.projects.chexmix.shapealign.alignment;

/**
 * PearsonCorrelation : pearson correlation and weighted pearson correlation
 * weights are determined by the average counts of the position
 * 
 * @author naomi yamada
 *
 */

public class PearsonCorrelation {	
	protected double [][] normAarray;
	protected double [][] normBarray;
	protected double [][] revNormBarray;
	protected int window;	
	protected boolean revAlign = false; 
	
	public PearsonCorrelation(double arrayA[][], double arrayB[][]){	
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
	
	public double computePearsonCorr(){		
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

		double sumA = 0, sumB = 0;
		double aveA = 0, aveB = 0;		
		for (int i = 0 ; i < catAarray.length; i++){
			sumA += catAarray[i];
			sumB += catBarray[i];
		}
		aveA = sumA/(double)catAarray.length;
		aveB = sumB/(double)catBarray.length;
				
		// calculate Pearson correlation from forward and reverse pair
		double covf = 0, covr = 0;
		double varA = 0, varBf = 0, varBr = 0;
		double corrF = 0, corrR = 0;
		for (int i = 0; i < catAarray.length; i++){
			double ai = catAarray[i] - aveA;
			double bi_f = catBarray[i] - aveB;
			double bi_r = catRevBarray[i] - aveB;
			
			covf += ai*bi_f;
			covr += ai*bi_r;
			varA += ai*ai;
			varBf += bi_f*bi_f;
			varBr += bi_r*bi_r;
		}
		corrF = covf/(Math.sqrt(varA)*Math.sqrt(varBf));
		corrR = covr/(Math.sqrt(varA)*Math.sqrt(varBr));
		
		if (corrF > corrR){
			revAlign = false;
			return corrF;
		}else{
			revAlign = true;
			return corrR;
		}
	}
	
	public double computeWeightedPearsonCorr(){
		
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
		
		double [] fweight = new double [2*(window+1)];
		double [] rweight = new double [2*(window+1)];
		for (int i = 0; i <catAarray.length ; i++){
			fweight[i] = (catAarray[i]+catBarray[i])/2;
			rweight[i] = (catAarray[i]+catRevBarray[i])/2;
		}		

		double weightedSumA = 0, weightedSumB = 0, weightedSumB_r = 0;
		double weightSum = 0;
		double weightedAveA = 0, weightedAveB = 0, weightedAveB_r = 0;		
		for (int i = 0 ; i < catAarray.length; i++){
			weightedSumA += catAarray[i]*fweight[i];
			weightedSumB += catBarray[i]*fweight[i];
			weightedSumB_r += catRevBarray[i]*rweight[i];
			weightSum += fweight[i];
		}
		weightedAveA = weightedSumA/weightSum;
		weightedAveB = weightedSumB/weightSum;
		weightedAveB_r = weightedSumB_r/weightSum;
				
		// calculate Pearson correlation from forward and reverse pair
		double covf_w = 0, covr_w = 0;
		double varA_w = 0, varBf_w = 0, varBr_w = 0;
		double corrF_w = 0, corrR_w = 0;
		for (int i = 0; i < catAarray.length; i++){
			double ai = catAarray[i] - weightedAveA;
			double bi_f = catBarray[i] - weightedAveB;
			double bi_r = catRevBarray[i] - weightedAveB_r;
			
			covf_w += fweight[i]*ai*bi_f;
			covr_w += rweight[i]*ai*bi_r;
			varA_w += fweight[i]*ai*ai;
			varBf_w += fweight[i]*bi_f*bi_f;
			varBr_w += rweight[i]*bi_r*bi_r;
		}
		corrF_w = covf_w/(Math.sqrt(varA_w)*Math.sqrt(varBf_w));
		corrR_w = covr_w/(Math.sqrt(varA_w)*Math.sqrt(varBr_w));
		
		if (corrF_w > corrR_w){
			revAlign = false;
			return corrF_w;
		}else{
			revAlign = true;
			return corrR_w;
		}
	}
	
}
