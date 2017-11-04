package org.seqcode.projects.chexmix.stats;

public class DistanceMetrics {

	public static double Euclidean_Distance(double[] x1, double[] x2, double[] y1, double[] y2){
		if (x1.length!=y1.length || x2.length!=y2.length)
			System.err.println("unequal vector length!!");		
		double max_x=0; double max_y=0;
		for (int i=0; i< x1.length; i++){
			if (x1[i]> max_x){ max_x=x1[i];}
			if (x2[i]> max_x){ max_x=x2[i];}
		}
		for (int i=0; i< x1.length; i++){
			x1[i]=x1[i]/max_x;
			x2[i]=x2[i]/max_x;
		}
		for (int i=0; i< y1.length; i++){
			if (y1[i] > max_y){ max_y=y1[i];}
			if (y2[i] > max_y){ max_y=y2[i];}
		}
		for (int i=0; i< y1.length; i++){
			y1[i]=y1[i]/max_y;
			y2[i]=y2[i]/max_y;
		}
		double squaredEuclideanDistance =0;
		for (int i=0; i< x1.length; i++)
			squaredEuclideanDistance += Math.pow(x1[i]-y1[i], 2);
		for (int i=0; i< y1.length; i++)
			squaredEuclideanDistance += Math.pow(x2[i]-y2[i], 2);
		
		return Math.sqrt(squaredEuclideanDistance);
	}
}
