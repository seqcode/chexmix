package org.seqcode.projects.chexmix.composite;

import java.util.Map;

/**
 * CompositeModelComponent: a component in ChExMix mixtures
 * 
 * TODO: add per-replicate support
 * @author mahony
 *
 */
public class CompositeModelComponent  implements Comparable<CompositeModelComponent>{
	
	protected int index; //Index of the component
	protected String label; //Label is only used in reporting - has no other consequence in algorithm
	protected int position; //Positions wrt the left edge of the composite
	protected double pi; //Emission probability
	protected double sumRespW=0, sumRespC=0; //Sum of read responsibilities, per strand
	protected TagProbabilityDensity distrib; //Tag distribution associated with this model component
	protected boolean updatablePosition = true;
	protected boolean updatablePi = true;
	protected double[] tagProfileWatson=null;  //Read profiles are left uninitialized, as they are only used to store read distributions at the end of training. 
	protected double[] tagProfileCrick=null; //The read distributions stored in readProfiles are centered on the binding position. 
	
	
	public CompositeModelComponent(TagProbabilityDensity d, int pos, int index, String label, boolean updatablePosition, boolean updatablePi){
		distrib = d;
		pi=1;
		position = pos;
		this.index = index;
		this.label = label;
		this.updatablePosition=updatablePosition;
		this.updatablePi=updatablePi;
	}
	
	
	//Accessors
	public int getIndex(){return index;}
	public String getLabel(){return label;}
	public double getPi(){return pi;}
	public int getPosition(){return position;}
	public double getSumResponsibility(){return sumRespW+sumRespC;}
	public double getWatsonResponsibility(){return sumRespW;}
	public double getCrickResponsibility(){return sumRespC;}
	public TagProbabilityDensity getTagDistribution(){return distrib;}
	public boolean isNonZero(){return pi>0;}
	public boolean hasUpdatablePosition(){return updatablePosition;}
	public boolean hasUpdatablePi(){return updatablePi;}
	public void setPi(double p){pi=p;}
	public void setPosition(int p){position = p;}
	public void setSumResponsibilities(double sumRespW, double sumRespC) { this.sumRespW=sumRespW; this.sumRespC=sumRespC; }
	
	public double[] getTagProfile(boolean watsonStrand){
		if(watsonStrand)
			return tagProfileWatson;
		else
			return tagProfileCrick;
	}
	public void setTagProfile(double[] profile, boolean watsonStrand){
		if (tagProfileWatson==null){
			tagProfileWatson = new double[profile.length];
			tagProfileCrick = new double[profile.length];
		}
		if (watsonStrand)
			tagProfileWatson = profile;
		else 
			tagProfileCrick = profile;
	}
	
	public int compareTo(CompositeModelComponent x){
		return position - x.position;
	}
	
	public String toString(int offset){
		return String.format("%d\t%.5f", position-offset, pi);
	}
	
	/**
	 * Encode variables in a string
	 * @return
	 */
	public String saveString(){
		String out="#CompositeModelComponent,"+index+","+label+","+position+","+pi+","+distrib.getIndex()+","+
				updatablePosition+","+updatablePi+","+sumRespW+","+sumRespC+",\n";		
		return out;
	}
	
	/**
	 * Load a CompositeModelComponent using a String in the same format as used in saveString
	 * @param line: String to be converted
	 * @param tagDensities: TagProbabilityDensities mapped by index
	 * @return
	 */
	public static CompositeModelComponent load(String line, Map<Integer, TagProbabilityDensity> tagDensities){
		CompositeModelComponent cmc = null;
		String[] bits = line.split(",");
		if(bits.length!=10 || !bits[0].equals("#CompositeModelComponent")){
			System.err.println("CompositeModelComponent.load(): Unexpected format");
			System.exit(1);
		}else{
			Integer index = new Integer(bits[1]);
			String label = bits[2];
			Integer pos = new Integer(bits[3]);
			Double pi = new Double(bits[4]);
			Integer tagDistIndex = new Integer(bits[5]);
			Boolean upPos = new Boolean(bits[6]);
			Boolean upPi = new Boolean(bits[7]);
			Double sumRW = new Double(bits[8]);
			Double sumRC = new Double(bits[9]);
			if(!tagDensities.containsKey(tagDistIndex)){
				System.err.println("CompositeModelComponent.load(): TagProbabilityDensity not defined");
				System.exit(1);
			}else{
				TagProbabilityDensity tagDist = tagDensities.get(tagDistIndex);
				cmc = new CompositeModelComponent(tagDist, pos, index, label, upPos, upPi);
				cmc.setPi(pi);
				cmc.setSumResponsibilities(sumRW, sumRC);
			}
		}
		return cmc;
	}
}
