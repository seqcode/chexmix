package org.seqcode.projects.chexmix.mixturemodel;

import org.seqcode.genome.location.Point;

/**
 * BindingComponents are used in mixture models to represent potential binding events.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class BindingSubComponents implements Comparable<BindingSubComponents>{

	protected Point coord;  //Event coordinate
	protected Point[][] coords;  //Event coordinates
	protected int position; //Position without the chromosome name (for computational convenience)
	protected int[][] positions; //Position without the chromosome name (for computational convenience)
	protected double[][] tau; //Subtype probabilities
	protected double pi; //Emission probability
	protected double sumRespW=0, sumRespC=0; //Sum of read responsibilities, per strand
	protected double[][] readProfile_plus;  //Read profiles are left uninitialized, as they are only used to store read distributions at the end of training. 
	protected double[][] readProfile_minus; //The read distributions stored in readProfiles are centered on the binding position. 
	protected double[][][][] sub_readProfile_plus; 
	protected double[][][][] sub_readProfile_minus; 
	protected int index=0;
	
	public BindingSubComponents(Point pos, int numReps){
		coord=pos;
		position = coord.getLocation();
		positions = null;
		pi = 1;
		readProfile_plus= new double[numReps][];
		readProfile_minus= new double[numReps][];
		for(int r=0; r<numReps; r++){
			readProfile_plus[r]=null;
			readProfile_minus[r]=null;
		}	
		sub_readProfile_plus= new double[numReps][][][];
		sub_readProfile_minus= new double[numReps][][][];
		for(int r=0; r<numReps; r++){
			sub_readProfile_plus[r]=null;
			sub_readProfile_minus[r]=null;
		}
	}
	
	//Accessors
	public double getPi(){return pi;}
	public Point getCoord(){return coord;}
	public int getPosition(){return position;}
	public int[][] getPositions(){return positions;}
	public int getIndex(){return index;}
	public double[][] getTau(){return tau;}
	
	public boolean isNonZero(){return pi>0;}

	public double getSumResponsibility(){return sumRespW+sumRespC;}
	
	public double[] getReadProfile_plus(int repIndex){
		return readProfile_plus[repIndex];
	}
	
	public double[] getReadProfile_minus(int repIndex){
		return readProfile_minus[repIndex];
	}	
	
	public double getReadProfile(int repIndex,int i, char strand){
		double result=0;
		if (strand=='+')
			result=readProfile_plus[repIndex][i];
		else if (strand=='-')
			result=readProfile_minus[repIndex][i];
		return result;
	}
	
	public int getMaxType(){
		double maxProb=0; int maxType=0;
		for (int bt=0; bt< tau.length; bt++){
            double typeTau=0;
            for (int s=0; s< tau[bt].length; s++)
                typeTau+=tau[bt][s];
            if (typeTau> maxProb){ maxProb=typeTau ; maxType=bt;}
        }
		return maxType;
	}
	
	public char getMaxTypeStrand(){
		if (tau[getMaxType()][0]> tau[getMaxType()][1]){ return '+';}
		else { return '-';}
	}
	
	public double getMaxTypeProbability(){
		return (tau[getMaxType()][0]+tau[getMaxType()][1]);
	}
	
	//Mutators
	public void setPi(double p){pi=p;}
	public void setPosition(int p){position = p; updateCoordFromLocation();}
	public void setCoord(Point p){coord=p; position=p.getLocation();}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setTau(double[][] t){tau=t;}
	public void setSumResponsibilities(double sumRespW, double sumRespC) { this.sumRespW=sumRespW; this.sumRespC=sumRespC; }
	public void setUnstrandedSumResponsibilities(double sumResp) { this.sumRespW=sumResp/2; this.sumRespC=sumResp/2; }
	
	public void setPositions(int[][] p){
		positions = p;
		coords = new Point[p.length][2];
		for (int bt=0; bt< p.length;bt++){
			for (int s=0; s< p[bt].length; s++){
				Point newCoord = new Point(coord.getGenome(), coord.getChrom(), positions[bt][s]); 
				coords[bt][s]= newCoord;
			}
		}
	}
	
	public void uniformInit(double initValue){
		pi=initValue;
	}
	
	public void setReadProfile(int repIndex, double[] profile, char strand){
		if (readProfile_plus[repIndex]==null){
			readProfile_plus[repIndex] = new double[profile.length];
			readProfile_minus[repIndex] = new double[profile.length];
		}
		if (strand=='+')
			readProfile_plus[repIndex] = profile;
		else if (strand=='-')
			readProfile_minus[repIndex] = profile;
	}
	
	public void setSubReadProfile(int repIndex, double[][][] profile, char strand){
		if (sub_readProfile_plus[repIndex]==null){
			sub_readProfile_plus[repIndex] = new double[profile.length][profile[0].length][profile[0][0].length];
			sub_readProfile_minus[repIndex] = new double[profile.length][profile[0].length][profile[0][0].length];
		}
		if (strand=='+')
			sub_readProfile_plus[repIndex] = profile;
		else if (strand=='-')
			sub_readProfile_minus[repIndex] = profile;
	}

	//Comparable default method
	public int compareTo(BindingSubComponents m) {
		return getCoord().compareTo(m.getCoord());
	}
	
	//Compare by responsibility
	public int compareByResp(BindingSubComponents m) {
		return Double.compare(sumRespW+sumRespC, m.sumRespW+m.sumRespC);
	}
	
	//Compare by tau probabilities
	public int compareByTauProb(BindingSubComponents m){
		return Double.compare(getMaxTypeProbability(), m.getMaxTypeProbability());
	}
	
	public String toString(){
		String out ="B\t"+coord.getLocationString()+"\t"+String.format("%.3f",pi)+"\t"+String.format("%.3f", sumRespW+sumRespC)+"\t"+index;
		if (positions!=null)
			for (int bt=0; bt< positions.length; bt++)
				for (int s=0; s< 2; s++)
					out += "\t"+coords[bt][s]+"\t"+String.format("%.3f",tau[bt][s]);
		return out;
	}
}//end of BindingComponent class
