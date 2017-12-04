package org.seqcode.projects.chexmix.events;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.chexmix.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;

/**
 * BindingSubtype class represent binding subtypes that are associated with read distribution and optionally motif.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class BindingSubtype {
	ExperimentCondition condition;
	TagProbabilityDensity[] repBindingModel;
	List<StrandedPoint> modelReferences;
	WeightMatrix motif ;
	WeightMatrix freqMatrix;
	int motifOffset;
	int numEvents;
	boolean hasMotif;
	boolean reverseMotif;
	boolean isClusterProfile;
	
	public BindingSubtype(ExperimentCondition cond, List<StrandedPoint> modelRefs, int bindingModelWidth) throws Exception{
		condition = cond;
		modelReferences=modelRefs;
		motif = null;
		freqMatrix=null;
		motifOffset=0;
		hasMotif=false;
		reverseMotif=false;		
		isClusterProfile=false;		
		numEvents=modelRefs.size();
		// Make bindingModel
		makeTagDistribution(bindingModelWidth);	
	}
	
	public BindingSubtype(ExperimentCondition cond, TagProbabilityDensity model, int events){
		condition = cond;
		modelReferences=null;
		motif = null;
		freqMatrix=null;
		motifOffset=0;
		hasMotif=false;
		reverseMotif=false;		
		isClusterProfile=false;	
		numEvents=events;
		setBindingModel(model);
	}
	
	public TagProbabilityDensity getBindingModel(int repIndex){return repBindingModel[repIndex];}
	public WeightMatrix getMotif(){return motif;}
	public WeightMatrix getFreqMatrix(){return freqMatrix;}
	public int getMotifOffset(){return motifOffset;}
	public boolean reverseMotif(){return reverseMotif;}
	public int getNumEvents(){return numEvents;}
	public boolean hasMotif(){return hasMotif;}
	public boolean isClustered(){return isClusterProfile;}
	
	public void setClusteredProfile(boolean clust){ isClusterProfile=clust;}
	public void setMotif(WeightMatrix wm, WeightMatrix fm){ motif=wm; freqMatrix = fm; hasMotif=true;}
	public void setMotifOffset(int offset){motifOffset=offset;}
	public void setReverseMotif(boolean reverse){reverseMotif=reverse;}

	public void makeTagDistribution(int modelWidth) throws Exception{
		//Build the composite distribution(s)
		CompositeTagDistribution signalComposite = new CompositeTagDistribution(modelReferences, condition, modelWidth,true);							
		TagProbabilityDensity model = new TagProbabilityDensity(signalComposite.getWinSize()-1);
		model.loadData(signalComposite.getCompositeWatson(), signalComposite.getCompositeCrick());
		setBindingModel(model);
	}
	
	public void setBindingModel(TagProbabilityDensity model){
		int repSize = condition.getReplicates().size();
		repBindingModel = new TagProbabilityDensity[repSize];
		for (int r=0; r< repSize; r++)
			repBindingModel[r]=model;
		
	}

}
