package org.seqcode.projects.chexmix.shapealign.misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.Pair;
import org.seqcode.motifs.MemeER;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignment;

import net.sf.javaml.clustering.FarthestFirst;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.*;
import net.sf.javaml.distance.EuclideanDistance;


/**
 * ClusterMotifFinder uses java-lm to do Kmeans clustering on shape features
 * This implementation was slow.  I'm leaving this as an example of how to use java-ml.
 * input : reference point
 * 
 * @author naomi yamada
 */

public class ClusterMotifFinderML {
	protected ExperimentManager manager;
	protected ShapeAlignConfig shapeAlignConfig;
	protected GenomeConfig gconfig;
	protected SequenceGenerator<Region> seqgen;
	protected ShapeAlignment shapeAlign;
	
//	protected double[][] distanceMatrix;
	protected List<StrandedRegion> strandedRegions = new ArrayList<StrandedRegion>();
	protected int numClusters;
//	protected Map<Instance, StrandedRegion> shapeDistanceMap = new HashMap<Instance, StrandedRegion>(); 
	protected Set<HashSet<StrandedRegion>> regionClusterSet = new HashSet<HashSet<StrandedRegion>>();
	
	protected File basedir_meme;	
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	
	/** The ROC value of the identified motifs that indicates their significance */
	protected HashMap<String, Double> motifsRocs = new HashMap<String,Double>();
	protected String[] randomSequences = new String[MemeER.MOTIF_FINDING_NEGSEQ];
	
	public ClusterMotifFinderML(ExperimentManager man, ShapeAlignConfig shapeConfig, GenomeConfig gcon) throws FileNotFoundException, UnsupportedEncodingException{
		manager = man;
		shapeAlignConfig = shapeConfig;
		numClusters = shapeConfig.getK();		
		gconfig = gcon;
		seqgen = new SequenceGenerator<Region>(gcon.getGenome());
		shapeAlign = new ShapeAlignment(shapeAlignConfig, gconfig, manager); 
		shapeAlign.excute();
		strandedRegions = shapeAlign.getStrandedRegions();
	}
	
	public void execute(){
		Dataset data = loadData(shapeAlign.getSimilarityMatrix());
		System.out.println("loaded the data.");
		
		// test
		System.out.println(data.get(0));
		
		Dataset[] clusters= Kmeans(data);
		System.out.println("clustering done.");
		printClusterSize(clusters);
		runMEME(clusters);
	}
	
	public void setRandomRegs(){
		List<Region> randomRegions = MemeER.randomRegionPick(gconfig.getGenome(), null, MemeER.MOTIF_FINDING_NEGSEQ,shapeAlignConfig.getWindowSize());
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = seqgen.execute(randomRegions.get(i));
		}
	}
	
	public Dataset loadData(double[][] distanceMatrix){		
		Dataset data = new DefaultDataset();
		for (int i = 0 ; i < distanceMatrix.length; i++){
			data.add(new DenseInstance(distanceMatrix[i], strandedRegions.get(i)));
		}	
		return data;
	}
	
	public Dataset[] Kmeans(Dataset data){		
		KMeans kmeans = new KMeans(numClusters);
		System.out.println("initialized Kmeans with K "+numClusters);
		System.out.println("data size "+data.size());
		Dataset[] clusters = kmeans.cluster(data);	
		return clusters;
	}
	
	public Dataset[] KCenter(Dataset data){	//aka Farthest-First
		FarthestFirst ff = new FarthestFirst(numClusters, new EuclideanDistance());
		Dataset[] clusters = ff.cluster(data);			
		return clusters;
	}
	
	public void printClusterSize(Dataset[] clusters){
		for (int i = 0 ; i < clusters.length;i++){
			System.out.print(clusters[i].size()+"\t");
		}
		System.out.println();
	}
	
	public void runMEME(Dataset[] clusters){
		
		List<ArrayList<String>> seqClusterList = new ArrayList<ArrayList<String>>();
		// get regions and sequences from clustered instances
		for (int i = 0 ; i < clusters.length;i++){
			ArrayList<StrandedRegion> clusterRegions = new ArrayList<StrandedRegion>();
			ArrayList<String> clusterSeq = new ArrayList<String>();
			clusters[i].classes();
			for (Instance inst : clusters[i]){
				System.out.println("instance val is : "+inst.classValue());
				clusterRegions.add((StrandedRegion) inst.classValue());
				System.out.println("seq is : "+seqgen.execute((StrandedRegion) inst.classValue()));
				clusterSeq.add(seqgen.execute((StrandedRegion) inst.classValue()));
			}
			regionClusterSet.add(new HashSet<StrandedRegion>(clusterRegions));
			seqClusterList.add(new ArrayList<String>(new HashSet<String>(clusterSeq)));
		}		
		// run meme using MEMEer
		String memeargs = shapeAlignConfig.getMEMEargs();
		MemeER meme = new MemeER(shapeAlignConfig.getMemePath(), memeargs);		
		basedir_meme = new File(shapeAlignConfig.getOutDir().getAbsoluteFile()+File.separator);
		basedir_meme.mkdirs();
		
		for (int i = 0 ; i < seqClusterList.size(); i ++){
			List<String> seqcluster = seqClusterList.get(i);
			File meme_outFile = new File(basedir_meme+File.separator+"meme_clust"+(i+1));
			System.out.println("output directory is "+meme_outFile);
			System.out.println("Running meme now for cluster "+(i+1));
			Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqcluster, meme_outFile, false);
			System.out.println("Finished running meme for cluster "+(i+1)+", Now evaluating motif significance");
			List<WeightMatrix> wm = matrices.car();
			List<WeightMatrix> fm = matrices.cdr();
			// Index for all the selected motifs
			int motInd = 0;
			
			if(wm.size()>0){
				//Evaluate the significance of the  motifs
				double rocScores[] = meme.motifROCScores(wm,seqcluster,randomSequences);
				System.err.println("MEME results for cluster "+(i+1));
				for(int w=0; w<fm.size(); w++){
					if(fm.get(w)!=null){
						System.err.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
					}
					if(rocScores[w] > meme.getMotifMinROC()){
						//selectedMotifs.add(fm.get(w));
						fm.get(w).setName("cluster"+(i+1)+"_"+Integer.toString(motInd));
						motifs.add(fm.get(w));
						motifsRocs.put("cluster"+(i+1)+"_"+Integer.toString(motInd), rocScores[w]);
					}
					motInd++;
				}
			}
		}		
	}
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		ShapeAlignConfig shapeAlignConf = new ShapeAlignConfig(args);	
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);	
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);
//		ShapeAlignment profile = new ShapeAlignment(shapeAlignConf, gconf, manager); 	
		ClusterMotifFinderML motifFinder = new ClusterMotifFinderML(manager,shapeAlignConf,gconf);
		motifFinder.execute();
				
		manager.close();		
	}

}
