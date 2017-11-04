package org.seqcode.projects.chexmix.shapealign.progressivealignment;

import org.seqcode.genome.location.Point;

/**
 * Profile represents a merged ChIP-exo reads from regions or other profiles
 * @author nuy11
 *
 */
public class Profile implements Comparable<Profile> {
	
	protected Point coord; //Event coordinate
	protected int position; //Position without the chromosome name (for computational convenience)
	protected char strand; //Event Strand 
	protected int index;
	protected int offset;
	protected double[] watsonTags; // watson tag counts
	protected double[] crickTags; // crick tag counts
	protected double[] normWatson; // normalized watson tag counts
	protected double[] normCrick; // normalized crick tag counts
	protected int indexA;	// smaller index is stored at indexA
	protected int indexB;	// indexB is bigger than indexA
	protected double similarity;
	protected boolean reverse;
	protected boolean isProfile=false;	//boolean to indicate if the profile is a profile or original region
	
	public Profile(double[] watson, double[] crick, int index){
		watsonTags = watson;
		crickTags = crick;
		this.index = index;	
	}
	
	//Accessors
	public Point getCoord(){return coord;}
	public int getPosition(){return position;}
	public char getStrand(){return strand;}
	public int getIndex(){return index;}
	public int getOffset(){return offset;}
	public double getSimilarity(){return similarity;}
	public int getIndexA(){return indexA;}
	public int getIndexB(){return indexB;}
	public double[] getWatsonProfile(){return watsonTags;}
	public double[] getCrickProfile(){return crickTags;}
	public boolean isReverse(){return reverse;}
	public boolean isProfile(){return isProfile;}
	
	//Mutators
	public void setPairProfile(int indexOne, int indexTwo, double sim){
		this.indexA = indexOne < indexTwo ? indexOne : indexTwo;
		this.indexB = indexOne < indexTwo ? indexTwo : indexOne;
		similarity = sim;
	}
	
	public void setPosition(int p){position = p; updateCoordFromLocation();}
	public void setCoord(Point p){coord=p; position=p.getLocation();}
	public void setStrand(char s){strand =s;}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setOffset(int off){offset=off;}
	public void setReverse(boolean isReverse){reverse=isReverse;}
	public void setProfile(boolean profile){isProfile=profile;}
	
	// normalize tag profile 
	public void normalizeProfile(){
		normWatson = new double[watsonTags.length];
		normCrick = new double[crickTags.length];		
		double sum = 0;
		for (int i=0; i < watsonTags.length; i++)
			sum += watsonTags[i];
		for (int i=0; i < crickTags.length; i++)
			sum += crickTags[i];
		// Normalize
		for (int i=0; i < watsonTags.length; i++)
			normWatson[i] = watsonTags[i]/sum;
		for (int i=0; i < crickTags.length; i++)
			normCrick[i] = crickTags[i]/sum;
	}
	
	//Comparable default method
	public int compareTo(Profile pair){
		double sim = pair.getSimilarity();
		if (sim > similarity)
			return 1;
		else
			return -1;
	}
	
}
