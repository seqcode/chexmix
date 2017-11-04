package org.seqcode.projects.chexmix.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;
import org.seqcode.math.stats.StatUtil;
import org.seqcode.projects.chexmix.composite.ProteinDNAInteractionModel;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;

/**
 * Classify binding subtypes using KL divergence using tag landscapes of a given window
 * @author nuy11
 *
 */
public class BindingTypeClassifier {
	protected ExperimentManager exptMan;	
	protected List<StrandedPoint> points;
	protected int win;
	protected int smoothWin;
	protected int centerOffset;
	protected int numConditions;
	protected int numPoints;
	protected List<TagProbabilityDensity> modelDensities;
	
	public BindingTypeClassifier(List<StrandedPoint> points, ExperimentManager eMan, int win, int smoothingWindow, List<TagProbabilityDensity> tagProbDensities ){
		exptMan = eMan;
		this.win = win;
		smoothWin = smoothingWindow;
		centerOffset = win/2;
		this.numConditions=exptMan.getNumConditions();
		this.points = points;
		numPoints = points.size();
		modelDensities = tagProbDensities;
	}
	
	public void execute(){
		
		CompositeTagDistribution maker = new CompositeTagDistribution(points, exptMan, win, true);		
		for(ExperimentCondition cond : exptMan.getConditions()){
			String compositeFileName = "out_composite."+cond.getName()+".txt";
			maker.printProbsToFile(cond, compositeFileName);
		}
		
		for(ExperimentCondition cond : exptMan.getConditions()){			
			String KLFileName = "out_kl."+cond.getName()+".txt";
			try {
				FileWriter fout = new FileWriter(KLFileName);
				for (StrandedPoint p : points){
					double[] watsonTags = maker.getPointWatson(p, cond);
					double[] crickTags = maker.getPointCrick(p,  cond);
					if (smoothWin >0){	// Smooth tags
						watsonTags = StatUtil.gaussianSmoother(watsonTags, smoothWin);
						crickTags = StatUtil.gaussianSmoother(crickTags, smoothWin);
					}
					// Measure KL divergence with sliding window
					double[] minLogKLs = new double[modelDensities.size()];
					for (int i=0; i<modelDensities.size() ; i++)
						minLogKLs[i]=Double.MAX_VALUE;
					for (int modelType=0; modelType< modelDensities.size(); modelType++){
						int modelWidth = modelDensities.get(modelType).getWinSize();
						double[] modelW = modelDensities.get(modelType).getWatsonProbabilities();
						double[] modelC = modelDensities.get(modelType).getCrickProbabilities();
						double currLogKL = 0;     			
						for (int start=0; start<=(win-modelWidth); start++){
							currLogKL = 0;		
							//copy current window
							double[] currW = new double[modelWidth];
							double[] currC = new double[modelWidth];
							for (int offset=0; offset< modelWidth; offset++){
								currW[offset]=watsonTags[offset+start];
								currC[offset]=crickTags[offset+start];
							}							
							
							//Calc KL
							currLogKL += StatUtil.log_KL_Divergence(modelW, currW) + StatUtil.log_KL_Divergence(currW, modelW);
							currLogKL += StatUtil.log_KL_Divergence(modelC, currC) + StatUtil.log_KL_Divergence(currC, modelC);

							if (currLogKL < minLogKLs[modelType])
								minLogKLs[modelType]=currLogKL;
							
							// calculate divergence in reversed tags
        					currLogKL = 0;
							//copy current window
        					double[] rcurrW = new double[modelWidth];
        					double[] rcurrC = new double[modelWidth];
        					for (int i=0; i < modelWidth; i++){
        						rcurrW[i]=currW[modelWidth-i-1];
        						rcurrC[i]=currC[modelWidth-i-1];
        					}
					
        					//Calc KL
        					currLogKL += StatUtil.log_KL_Divergence(modelW, rcurrW) + StatUtil.log_KL_Divergence(rcurrW, modelW);
        					currLogKL += StatUtil.log_KL_Divergence(modelC, rcurrC) + StatUtil.log_KL_Divergence(rcurrC, modelC);

        					if (currLogKL < minLogKLs[modelType])
        						minLogKLs[modelType]=currLogKL; 						
						} 
					}	
					// print KL values
					double minLogKLVal=Double.MAX_VALUE;
					int minLogKLType=0;
					for (int modelType=0; modelType< modelDensities.size(); modelType++)
						fout.write(minLogKLs[modelType]+"\t");						
					for (int modelType=0; modelType< modelDensities.size(); modelType++){	
						if (minLogKLs[modelType] < minLogKLVal){
							minLogKLVal=minLogKLs[modelType];
							minLogKLType=modelType;
						}	
					}
					fout.write(p.toString()+"\t"+minLogKLType);						
					fout.write("\n");
				}
				fout.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}		
	}

	//Main method to make new composite distributions
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if(args.length==0){
			System.err.println("BindingTypeClassifier:"+
					"\t--cpoints <stranded point file>"+
					"\t--cwin <window around points>"+
					"Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"\t--smooth <gaussian smoothing step>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			int w = Args.parseInteger(args, "cwin", 1100);
			int swindow = Args.parseInteger(args, "smooth", 0);
			String pFile = Args.parseString(args, "cpoints", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
			
			ChExMixConfig config = new ChExMixConfig(gcon, args);			
			String mFileA = Args.parseString(args, "modelA", null);
			String mFileB = Args.parseString(args, "modelB", null);
			
			ProteinDNAInteractionModel modelA = ProteinDNAInteractionModel.loadFromFile(config, new File(mFileA));
			ProteinDNAInteractionModel modelB = ProteinDNAInteractionModel.loadFromFile(config, new File(mFileB));
			
			List<TagProbabilityDensity> tagProbDensities = new ArrayList<TagProbabilityDensity>();		
			tagProbDensities.add(modelA.makeTagProbabilityDensityFromAllComponents());
			tagProbDensities.add(modelB.makeTagProbabilityDensityFromAllComponents());			
			
			BindingTypeClassifier classifier = new BindingTypeClassifier(pts, manager, w, swindow, tagProbDensities);
			classifier.execute();
						
			manager.close();
		}
	}	
}
