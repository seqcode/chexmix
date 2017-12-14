package org.seqcode.projects.chexmix.shapealign.alignment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.ml.clustering.Clusterable;
import org.seqcode.ml.clustering.affinitypropagation.APCluster;
import org.seqcode.ml.clustering.affinitypropagation.MatrixSimilarityMeasure;
import org.seqcode.ml.clustering.affinitypropagation.SimilarityMeasure;
import org.seqcode.ml.clustering.affinitypropagation.SimilarityMeasure.APExemplar;

public class ShapeClusterTesting {

	protected double preferenceValue = -0.5;
	
	ShapeAlignConfig shapeAlignConf;
	GenomeConfig gconf;
	ExperimentManager manager;
	
	public ShapeClusterTesting(ShapeAlignConfig shapeAlignConf, GenomeConfig gconf, ExperimentManager manager){
		this.shapeAlignConf = shapeAlignConf;
		this.gconf = gconf;
		this.manager = manager;
		
	}
	
	public void setPreferenceValue(double p){
		preferenceValue=p;
	}
	
	public void execute(){
		//Run the aligner
		ShapeAlignmentTesting profile = new ShapeAlignmentTesting(shapeAlignConf, gconf, manager); 	
		try {
			profile.execute();
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}

		//Get the alignment matrix
		double[][] alignmentScores = profile.getSimilarityMatrix();
		List<StrandedRegion> regs = profile.getStrandedRegions();
		List<String> regNames= new ArrayList<String>();
		for(StrandedRegion sr : regs)
			regNames.add(sr.getLocationString());
		
		//Run AffinityPropagation
		MatrixSimilarityMeasure<Clusterable> msm = new MatrixSimilarityMeasure<Clusterable>(regNames, 
				alignmentScores, preferenceValue);
		double netsim = APCluster.cluster(msm.objects(), msm, 0.5, 50, 500);
		List<SimilarityMeasure<Clusterable>.APExemplar> exemplars = msm.getExemplars();
		List<SimilarityMeasure<Clusterable>.APAssignment> assignments = msm.getAssignments();
		
		int numExemplars = exemplars.size();
		Map<Integer, List<Integer>> clusters = 
				new HashMap<Integer, List<Integer>>();
		for(SimilarityMeasure<Clusterable>.APExemplar e : exemplars){
			clusters.put(e.index, new ArrayList<Integer>());
			clusters.get(e.index).add(e.index);
		}
		
		for(SimilarityMeasure<Clusterable>.APAssignment a : assignments){
			//System.out.println(a.index+"\t"+a.name+"\t"+a.exemplar.index);
			clusters.get(a.exemplar.index).add(a.index);
		}
		
		System.out.println("Affinity Propagation exemplars:");
		for(SimilarityMeasure<Clusterable>.APExemplar e : exemplars){
			System.out.println(e.index+"\t"+e.name+"\t"+clusters.get(e.index).size()+" members");
		}
		
		
		//OUTPUT
		
		//make outdir
		shapeAlignConf.makeOutputDir();
		File outDir = shapeAlignConf.getOutDir();
		
		//Get composites for each cluster
		int window = shapeAlignConf.getWindowSize();
		for(Integer c : clusters.keySet()){
			for (ExperimentCondition condition : manager.getConditions()){		
				for (ControlledExperiment rep: condition.getReplicates()){	
					double[][] composite = profile.getExemplarBasedAlignment(rep, c, clusters.get(c));
					
					try {
						File compFlie = new File(outDir.getAbsolutePath()+File.separator+"cluster"+c+".composite.txt");
						FileWriter fout = new FileWriter(compFlie);
						fout.write("#Cluster"+c+"\n");
						for(int i=0; i<=window; i++){
							fout.write(composite[i][0]+"\t"+composite[i][1]+"\n");
						}
						fout.close();
					} catch (IOException e1) {
						e1.printStackTrace();
					}
					
				}
			}
		}
		
		//Get aligned points for each cluster
		for(Integer c : clusters.keySet()){
			List<StrandedPoint> spts = profile.getExemplarBasedAlignedPoints(c, clusters.get(c));
			
			try {
				File compFlie = new File(outDir.getAbsolutePath()+File.separator+"cluster"+c+".points");
				FileWriter fout = new FileWriter(compFlie);
				fout.write("#Cluster"+c+"\n");
				for(StrandedPoint s :spts)
					fout.write(s.getLocationString()+"\n");
				fout.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException{
		
		ShapeAlignConfig shapeAlignConf = new ShapeAlignConfig(args);	
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);	
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);

		
		ShapeClusterTesting cluster = new ShapeClusterTesting(shapeAlignConf, gconf, manager);
		double preferenceValue = -0.5;
		ArgParser ap = new ArgParser(args);
		if(ap.hasKey("preferenceValue")){
			preferenceValue = new Double(ap.getKeyValue("preferenceValue"));
			cluster.setPreferenceValue(preferenceValue);
		}
		cluster.execute();
		
		manager.close();		
	}
	
}
