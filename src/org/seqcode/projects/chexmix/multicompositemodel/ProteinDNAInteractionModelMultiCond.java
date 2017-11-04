package org.seqcode.projects.chexmix.multicompositemodel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.CompositeModelComponent;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;

/**
 * ProteinDNAInteractionModel: defines a collection of composite model components 
 * 
 *   Coordinates are 0-based, but the midpoint of the model is defined by the centerOffset variable. 
 * @author mahony
 *
 */
public class ProteinDNAInteractionModelMultiCond {

	protected int modelWidth;
	protected int centerOffset;
	protected int numConditions; // number of conditions compositeProfiles are based on
	protected double[][] compositeWatson; //per-condition compositeProfile
	protected double[][] compositeCrick; //per-condition compositeProfile
	protected ChExMixConfig cmConfig;
	protected List<List<CompositeModelComponent>> allComponents;
	protected List<List<CompositeModelComponent>> XLComponents;
	protected List<Map<CompositeModelComponent, TagProbabilityDensity>> XLComponentDensities;
	protected List<CompositeModelComponent> CSComponent;
	protected List<CompositeModelComponent> backgroundComponent;
	protected int numXLComponents; 
	
	/**
	 * Initialize a new ProteinDNAInteractionModel
	 * @param config
	 * @param width
	 * @param initXLdistrib
	 * @param initCSdistrib
	 * @param initBackDistrib
	 * @param noisePi
	 */
	public ProteinDNAInteractionModelMultiCond(ChExMixConfig config, CompositeTagDistribution composite, TagProbabilityDensity initXLdistrib, TagProbabilityDensity initCSdistrib, 
			TagProbabilityDensity initBackDistrib, double noisePi){
		modelWidth = composite.getWinSize();
		centerOffset = composite.getCenterOffset();
		numConditions = composite.getNumConditions();
		compositeWatson = composite.getCompositeWatson();
		compositeCrick = composite.getCompositeCrick();
		cmConfig = config;
		
		allComponents = new ArrayList<List<CompositeModelComponent>>();
		XLComponents = new ArrayList<List<CompositeModelComponent>>();
		XLComponentDensities = new ArrayList<Map<CompositeModelComponent,TagProbabilityDensity>>();
		CSComponent = new ArrayList<CompositeModelComponent>();
		backgroundComponent = new ArrayList<CompositeModelComponent>();
		
		// Initialize
		for (int c=0; c < numConditions; c++){
			allComponents.add(new ArrayList<CompositeModelComponent>());
			XLComponents.add(new ArrayList<CompositeModelComponent>());
			XLComponentDensities.add(new HashMap<CompositeModelComponent,TagProbabilityDensity>());
		}		
		
		numXLComponents = (initCSdistrib.getInfluenceRange()/2)/cmConfig.getXLComponentSpacing();
		
		for (int c=0; c < numConditions; c++){
			//Background component
			backgroundComponent.add(new CompositeModelComponent(initBackDistrib, centerOffset, 0, "Back",  false, true));

			//ChIP-seq component
			CSComponent.add(new CompositeModelComponent(initCSdistrib, centerOffset, 1, "CS", false, true));
				
			//XL components
			int xlPos = centerOffset-(numXLComponents*cmConfig.getXLComponentSpacing())/2;
			for(int i=0; i<numXLComponents; i++){
				TagProbabilityDensity XLTagDist = initXLdistrib.clone();
				CompositeModelComponent newComp = new CompositeModelComponent(XLTagDist, xlPos, i+2, "XL", true, true);
				XLComponents.get(c).add(newComp);
				XLComponentDensities.get(c).put(newComp, XLTagDist);
				xlPos += cmConfig.getXLComponentSpacing(); 
			}
			
			//All components
			allComponents.get(c).add(backgroundComponent.get(c));
			allComponents.get(c).add(CSComponent.get(c));
			allComponents.get(c).addAll(XLComponents.get(c));
		}
		
		//Set initial pi values
		setInitialPi(noisePi);
	}
	
	/**
	 * Initialize a ProteinDNAInteractionModel using existing components (e.g. when loading saved model)
	 * @param config
	 * @param width
	 * @param CSComp
	 * @param backComp
	 * @param XLComps
	 */
	public ProteinDNAInteractionModelMultiCond(ChExMixConfig config, int width, int centerOff, int numConds, double[][] compWatson, double[][] compCrick, List<CompositeModelComponent> CSComp, List<CompositeModelComponent> backComp, List<List<CompositeModelComponent>> XLComps){
		modelWidth = width;
		centerOffset = centerOff;
		numConditions = numConds;
		compositeWatson = compWatson;
		compositeCrick = compCrick;
		cmConfig = config;
		backgroundComponent = backComp;
		CSComponent = CSComp;
		XLComponents = XLComps;
		// Hack here
		numXLComponents = XLComponents.get(0)==null ? 0 : XLComponents.get(0).size();

		allComponents = new ArrayList<List<CompositeModelComponent>>();
		for (int c=0; c < numConditions; c++){
			allComponents.add(new ArrayList<CompositeModelComponent>());
			allComponents.get(c).add(backgroundComponent.get(c));
			allComponents.get(c).add(CSComponent.get(c));
		}
		XLComponentDensities = new ArrayList<Map<CompositeModelComponent,TagProbabilityDensity>>();
		for (int c=0; c < numConditions; c++){
			XLComponentDensities.add(new HashMap<CompositeModelComponent,TagProbabilityDensity>());
			if(numXLComponents>0){
				allComponents.get(c).addAll(XLComponents.get(c));
				for(CompositeModelComponent xl : XLComponents.get(c))
					XLComponentDensities.get(c).put(xl, xl.getTagDistribution());
			}			
		}
	}
	
	//Accessors
	public int getWidth(){return modelWidth;}
	public int getCenterOffset(){return centerOffset;}
	public TagProbabilityDensity getXLTagDistribution(int c, CompositeModelComponent comp){return XLComponentDensities.get(c).get(comp);}
	public TagProbabilityDensity getCSTagDistribution(int c){return CSComponent.get(c).getTagDistribution();}
	public TagProbabilityDensity getBackTagDistribution(int c){return backgroundComponent.get(c).getTagDistribution();}
	public int getNumComponents(){return numXLComponents+2;}
	public List<CompositeModelComponent> getAllComponents(int c){return allComponents.get(c);}
	public List<CompositeModelComponent> getXLComponents(int c){return XLComponents.get(c);}
	public CompositeModelComponent getCSComponent(int c){return CSComponent.get(c);}
	public CompositeModelComponent getBackgroundComponent(int c){return backgroundComponent.get(c);}
	
	/**
	 * Return the non-zero components
	 * @return
	 */
	public List<CompositeModelComponent> getNonZeroComponents(int cond){
		List<CompositeModelComponent> comps = new ArrayList<CompositeModelComponent>();
		for(CompositeModelComponent c : allComponents.get(cond))
			if(c.isNonZero())
				comps.add(c);
		return comps;
	}
	
	/**
	 * Set the initial pi values for all components. 
	 * Assumes model is initialized. 
	 */
	protected void setInitialPi(double noisePi){
		for (int c=0; c < numConditions; c++){
			backgroundComponent.get(c).setPi(noisePi);
			double bindingPi = 1-noisePi;
			double initCS = cmConfig.noXL() ? bindingPi : Math.max(bindingPi*cmConfig.INIT_CS_TO_XL_RATIO, cmConfig.MIN_CS_PI);
			//Temporarily removed CS component
//			double initCS = 0;
			CSComponent.get(c).setPi(initCS);
			for(CompositeModelComponent xl : XLComponents.get(c)){
				double xlPi = cmConfig.noXL() ? 0.0 : bindingPi*(1-initCS)/(double)numXLComponents;
				xl.setPi(xlPi);
			}
		}
	}
	
	
	/**
	 * toString
	 */
	public String toString(int c){
		String out = "";
		Collections.sort(XLComponents.get(c));
		out = out+"ProteinDNAInteractionModel:\n";
		for(CompositeModelComponent xl : XLComponents.get(c)){
			if(xl.isNonZero())
				out = out+"\t"+xl.getIndex()+"\tXL:\t"+xl.toString(centerOffset)+"\n";
		}
		out = out+"\t"+CSComponent.get(c).getIndex()+"\tCS:\t"+CSComponent.get(c).toString(centerOffset)+"\n";
		out = out+"\t"+backgroundComponent.get(c).getIndex()+"\tBack:\t"+backgroundComponent.get(c).toString(centerOffset)+"\n";
		return out;
	}
	
	/**
	 * Save the entire model to a String
	 * @return
	 */
	public String saveString(int c){
		String out = "#ProteinDNAInteractionModel,"+modelWidth+","+centerOffset+","+numConditions+",\n";
		//Composite profiles
		out = out+"#CompositeWatson,"+c+",";
		for(int x=0; x<modelWidth; x++)
			out=out+compositeWatson[c][x]+",";
		out = out+"\n";
		out = out+"#CompositeCrick,"+c+",";
		for(int x=0; x<modelWidth; x++)
			out=out+compositeCrick[c][x]+",";
		out = out+"\n";
			
		//Tag distributions & Components
		// Here we need to either put all components out here (incl non-zero) or re-index non-zero components
		for(CompositeModelComponent comp : getAllComponents(c))
			out = out+comp.getTagDistribution().saveString();
		for(CompositeModelComponent comp : getAllComponents(c))
			out = out+comp.saveString();					
		return out;
	}
	
	/**
	 * Save the model to a file in the format specified by saveString()
	 * @param filename
	 */
	public void saveToFile(String filename, int c){
		try {
			FileWriter fout = new FileWriter(filename);
			fout.write(this.saveString(c));
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
}
