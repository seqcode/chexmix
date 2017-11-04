package org.seqcode.projects.chexmix.stats;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

/**
 * DistributionPercentile: Determine the percentile of the doubles  
 * 
 * @author naomi yamada
 *
 */

public class DistributionPercentile {
	
	protected double[] data;
	protected double percentile;
	protected double percentileVal;
	
	public DistributionPercentile(double[] array){
		data = array;
	}
	public DistributionPercentile(double[][] matrix){
		data = new double[matrix.length*matrix[0].length];
		int k = 0; 
		for (int i = 0 ; i < matrix.length; i++){
			for (int j = 0 ; j < matrix[0].length ; j++){
				data[k] = matrix[i][j];
				k++;
			}			
		}
	}
	
	public double getPercentileValue(double p){
		percentile = p;
		Percentile dist = new Percentile();
		return dist.evaluate(data, percentile);
	}
	
	public double getMaxVal(){
		double maxVal= -Double.MAX_VALUE;;
		for (int i=0; i< data.length; i++){
			if (data[i]>maxVal)
				maxVal=data[i];
		}
		return maxVal;
	}
}
