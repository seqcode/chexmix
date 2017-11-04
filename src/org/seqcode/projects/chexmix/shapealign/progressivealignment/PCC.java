package org.seqcode.projects.chexmix.shapealign.progressivealignment;

/**
 * PCC : computes Pearson Correlation with the forward and reverse profile
 * @author naomi
 */
public class PCC {	
	protected int pwindow;
	protected int calcwindow;
	protected double score;
	protected boolean revAlign = false; 
	protected int offset =0;
	protected double[] w_profile;	// Hold a merged watson profile using an offset at the maximum correlation
	protected double[] c_profile;	// Hold a merged crick profile using an offset at the maximum correlation
	
	//Accessors
	public double getMaxScore(){return score;}
	public boolean isReverse(){return revAlign;}
	public int getOffset(){return offset;}
	public double[] getWatsonProfile(){return w_profile;}
	public double[] getCrickProfile(){return c_profile;}
	
	public PCC (double[] watsonA, double[] crickA, double[] watsonB, double[] crickB, int profilewidth, int window){
		
		pwindow=profilewidth;	
		calcwindow = window;
		double maxScore = -1;
		int maxOffset = 0;
		boolean isReverse = false; 		
		
		System.out.println("within PCC");
		
		// Array A is fixed
		double[] normarrayA = new double[calcwindow*2];
		for (int i=0; i< calcwindow; i++){
			normarrayA[i] = watsonA[(pwindow-calcwindow)/2+i];
			normarrayA[calcwindow+i] = crickA[(pwindow-calcwindow)/2+i];
		}
		//Normalize and get average
		double sumA=0, sumAn=0, aveA=0;
		for (int j=0; j< calcwindow*2; j++)
			sumA += normarrayA[j];
		for (int j=0; j< calcwindow*2; j++)
			normarrayA[j]/=sumA;
		for (int j=0; j< calcwindow*2; j++)
			sumAn += normarrayA[j];
		aveA = sumAn/((double) calcwindow*2);
		
		for (int start = 0; start < (pwindow-calcwindow); start++){	
			int end = start+calcwindow;
			int off = start-(pwindow-calcwindow)/2;		
			// Copy
			double[] normarrayB = new double[calcwindow*2];
			double[] revNormarrayB = new double[calcwindow*2];
			int i=0;
			for (int w=start; w<end; w++){
				normarrayB[i] = watsonB[w];
				normarrayB[calcwindow+i] = crickB[w];
				revNormarrayB[calcwindow-i-1] = crickB[w];
				revNormarrayB[calcwindow*2-i-1] = watsonB[w];
				i++;
			}						
			// Normalize and get average
			double sumB=0, sumBn=0, aveB=0;
			for (int j=0; j< calcwindow*2; j++)
				sumB += normarrayB[j];
			for (int j=0; j< calcwindow*2; j++){
				normarrayB[j]/= sumB;
				revNormarrayB[j]/=sumB;
			}
			for (int j=0; j< calcwindow*2; j++)
				sumBn += normarrayB[j];
			aveB = sumBn/((double) calcwindow*2);
			
			// calculate Pearson correlation from forward and reverse pair
			double covf = 0, covr = 0;
			double varA = 0, varBf = 0, varBr = 0;
			for (int j = 0; j < calcwindow; j++){
				double ai = normarrayA[j] - aveA;
				double bi_f = normarrayB[j] - aveB;
				double bi_r = revNormarrayB[j] - aveB;				
				covf += ai*bi_f;
				covr += ai*bi_r;
				varA += ai*ai;
				varBf += bi_f*bi_f;
				varBr += bi_r*bi_r;
			}	
			
			double corrF = covf/(Math.sqrt(varA)*Math.sqrt(varBf));
			double corrR = covr/(Math.sqrt(varA)*Math.sqrt(varBr));			
			if (corrF > maxScore){
				maxScore = corrF;
				maxOffset = off; 
				isReverse = false;
			}
			if (corrR > maxScore){
				maxScore = corrR;
				maxOffset = off-1; 
				isReverse = true;
			}		
		}
		
		w_profile = new double[pwindow];
		c_profile = new double[pwindow];
		for (int i=0; i< pwindow;i++){
			w_profile[i]=0;
			c_profile[i]=0;		
		}
		
		for (int i=0; i < pwindow; i++){
			w_profile[i] += watsonA[i];
			c_profile[i] += crickA[i];
		}
		
		// merge tags
		if (!isReverse){
			for (int i=0; i< pwindow;i++){
				if ((i+maxOffset) >= 0 && (i+maxOffset) < pwindow){
					w_profile[i] += watsonB[i+maxOffset];
					c_profile[i] += crickB[i+maxOffset];
				}				
			}			
		}else{
			for (int i=0; i< pwindow;i++){
				if ((pwindow-i+maxOffset) >= 0 && (pwindow-i+maxOffset) < pwindow){
					
					// need to check this !!					
					w_profile[i] += crickB[pwindow-i+maxOffset];
					c_profile[i] += watsonB[pwindow-i+maxOffset];
				}				
			}	
		}		
		
		score = maxScore;
		revAlign = isReverse;
		offset = maxOffset;		
		
		System.out.println("revAlign? "+revAlign);
	}
}
