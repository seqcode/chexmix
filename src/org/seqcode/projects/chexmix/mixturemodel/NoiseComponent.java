package org.seqcode.projects.chexmix.mixturemodel;

import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.genome.location.Region;

/**
 * NoiseComponents are used in mixture models to soak up noise reads.
 * NoiseComponents have an emission probability but are position-less and generate reads along the region 
 * either uniformly or according to an empirical distribution. 
 *  
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class NoiseComponent {
	protected double pi; //emission probability
	protected double sum_resp=0; //sum of responsibilities
	protected double[][] distrib=null; //Expected distribution of noise reads, defined per replicate
	protected boolean uniformDistrib=true;
	protected Region reg;
	protected int regionWidth=0;
	
	/**
	 * Constructor
	 * @param emissionProb: probability of noise in the examined region (fixed)
	 * @param noiseDistrib: relative distribution of noise in the region (must be same length as width of region, indexed by replicate)
	 * @param reg: region that this noise component is being used in (only required for the width and the starting coordinate) 
	 */
	public NoiseComponent(double emissionProb, double[][] noiseDistrib, Region r, int numReps){
		pi = emissionProb; 
		reg = r;
		regionWidth=reg.getWidth();
		distrib = new double[numReps][];
		for(int h=0; h<numReps; h++){
			if(noiseDistrib[h]==null){
				uniformDistrib=true;
				distrib[h] = new double[r.getWidth()];
				double w = 1/(double)regionWidth;
				for(int i=0; i<regionWidth; i++)
					distrib[h][i]=w;
			}else{
				uniformDistrib=false;
				distrib[h] = noiseDistrib[h];
				if(distrib[h].length!=regionWidth){
					System.err.println("Noise distribution in "+reg+" is not appropriate width."); System.exit(1);
				}
			}
		}
	}
	
	public double getPi(){return pi;}
	public double getSumResponsibility(){return sum_resp;}
	
	public void setPi(double p){pi=p;}
	public void setSumResponsibility(double r){sum_resp=r;}
	
	//Scoring functions
	public double scoreHit(ReadHit h, int rep){
		int i = h.getFivePrime()-reg.getStart();
		if(i<0){i=0;}
		else if(i>=regionWidth){i=regionWidth;}
		return distrib[rep][i];
	}
	public double scoreBase(StrandedBaseCount b, int rep){
		int i= b.getCoordinate()-reg.getStart();
		if(i<0){i=0;}
		else if(i>=regionWidth){i=regionWidth;}
		return distrib[rep][i];
	}
	public double scorePosition(int pos, int rep){
		int i = pos-reg.getStart();
		if(i<0){i=0;}
		else if(i>=regionWidth){i=regionWidth;}
		return distrib[rep][i];
	}
	public double score(int i, int rep){
		if(i<0){i=0;}
		else if(i>=regionWidth){i=regionWidth;}
		return distrib[rep][i];
	}
	
	public String toString(){
		return "N\t"+String.format("%.3f", pi)+"\t"+regionWidth;
	}
	public String distribString(){
		String ds = new String();
		for(int d=0; d<regionWidth; d++)
			ds = ds+distrib[d]+"\n";
		return ds;
	}
}
