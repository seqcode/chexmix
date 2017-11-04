package org.seqcode.projects.chexmix.multicompositemodel;

import java.io.File;
import java.util.List;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;


public class CompositeXLFinderMultiCond {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig cconfig;
	protected CompositeModelMixtureMultiCond mixtureModel;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite;
	protected List<StrandedPoint> compositePoints;
	
	
	public CompositeXLFinderMultiCond(GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeChExMixOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public void execute(){
		//Load appropriate options
		compositePoints = cconfig.getCompositePoints();
		int winSize = cconfig.getCompositeWinSize();
		
		//Build the composite distribution(s)
		signalComposite = new CompositeTagDistribution(compositePoints, manager, winSize, true);
		//controlComposite = new CompositeTagDistribution(points, manager, winSize, false);
		controlComposite =null;
		
		//Initialize the mixture model 
		mixtureModel = new CompositeModelMixtureMultiCond(signalComposite, controlComposite, gconfig, econfig, cconfig, manager);
		
		//Train EM
		System.err.println("EM training");
		mixtureModel.trainEM();
		
		
		//ML assignment
		System.err.println("ML assignment");
		mixtureModel.assignML(false);
		
		//Report
		for(ExperimentCondition cond : manager.getConditions()){
			String compositeFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
			
			String perSiteRespFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_site-component-ML."+cond.getName()+".txt";
			mixtureModel.printPerSiteComponentResponsibilitiesToFile(cond, perSiteRespFileName);
			
			//Print component responsibility profiles
			String compProfileFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
			+"_component-profile."+cond.getName()+".txt";
			mixtureModel.printComponentProfilesToFile(cond, compProfileFileName);
		}
		mixtureModel.saveCompositePlots();
		
		//Save the model
		for(ExperimentCondition cond : manager.getConditions()){
			String modelFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()+"."+cond.getName()+".chexmix";
			mixtureModel.getModel().saveToFile(modelFileName, cond.getIndex());
		}
	}
	
	//Main
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("ChExMix version: "+ChExMixConfig.version);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		//econ.setLoadRead2(false);//Enforce for chip-exo
		ChExMixConfig ccon = new ChExMixConfig(gcon, args);
		if(ccon.helpWanted()){
			System.err.println(ccon.getArgsList());
		}else{
			CompositeXLFinderMultiCond xlFinder = new CompositeXLFinderMultiCond(gcon, econ, ccon);
			xlFinder.execute();
		}
		
	}
}
