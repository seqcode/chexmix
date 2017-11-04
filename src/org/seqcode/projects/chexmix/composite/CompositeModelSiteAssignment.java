package org.seqcode.projects.chexmix.composite;

import org.seqcode.genome.location.StrandedPoint;

public class CompositeModelSiteAssignment{

	protected StrandedPoint point;
	protected double totalTags; //Per-condition total tag counts for this point
	protected int[] modelComponentIndices;  //Indices of components that are active in model
	protected double[] modelComponentResponsibilities; //Per-condition, per-component responsibility tag counts
	
	public CompositeModelSiteAssignment(StrandedPoint pt, double total, int[] modelIndices, double[] modelResp){
		point = pt;
		totalTags = total;
		modelComponentIndices = modelIndices;
		modelComponentResponsibilities = modelResp;
	}
	
	public StrandedPoint getPoint(){return point;}
	
	public double getCompResponsibility(int componentIndex){
		int index=-1;
		for(int i=0; i<modelComponentIndices.length; i++)
			if(modelComponentIndices[i]==componentIndex)
				index = i;
		
		if(index==-1){
			System.err.println("CompositeModelSiteAssignment: error: component not indexed");
			System.exit(1);
		}
		return (modelComponentResponsibilities[index]);
	}
	
	public String toString(){
		String out = point.getLocationString()+String.format("\t%.2f", totalTags);
		for(int i=0; i<modelComponentIndices.length; i++)
			out = out+String.format("\t%.2f", modelComponentResponsibilities[i]);
		return out;
	}
	
}
