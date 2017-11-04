package org.seqcode.projects.chexmix.shapealign.alignment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.Pair;
import org.seqcode.ml.clustering.ClusterRepresentative;
import org.seqcode.ml.clustering.PairwiseElementMetric;
import org.seqcode.ml.clustering.kmeans.KMeansClustering;
import org.seqcode.ml.clustering.vectorcluster.DefaultVectorClusterElement;
import org.seqcode.ml.clustering.vectorcluster.EuclideanDistance;
import org.seqcode.ml.clustering.vectorcluster.Mean;
import org.seqcode.ml.clustering.vectorcluster.VectorClusterElement;
import org.seqcode.motifs.MemeER;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignConfig;
import org.seqcode.projects.chexmix.shapealign.alignment.ShapeAlignment;

/**
 * ClusterMotifFinder uses KMeansClustering from code base to cluster profile similarities
 * input : reference point
 * 
 * @author naomi yamada
 */

public class ClusterMotifFinder {
	protected ExperimentManager manager;
	protected ShapeAlignConfig shapeAlignConfig;
	protected GenomeConfig gconfig;
	protected SequenceGenerator<Region> seqgen;
	protected ShapeAlignment shapeAlign;
	
	private PairwiseElementMetric<VectorClusterElement> metric = new EuclideanDistance<VectorClusterElement>();
	private ArrayList<VectorClusterElement> profiles = new ArrayList<VectorClusterElement>();
	
	protected List<StrandedRegion> strandedRegions;
	protected int K;
	private int TRAIN_REPEATS=10;
	private int MAX_ITER = 100; // default 100 iterations
	private int numRandRegions = 1000;
	protected ArrayList<HashSet<StrandedRegion>> regionClusterList = new ArrayList<HashSet<StrandedRegion>>();
	
	protected File basedir_meme;	
	protected List<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
	
	/** The ROC value of the identified motifs that indicates their significance */
	protected HashMap<String, Double> motifsRocs = new HashMap<String,Double>();
	protected String[] randomSequences = new String[numRandRegions];
	
	// getters
	public List<WeightMatrix> getMotifs(){return motifs;}
	public HashMap<String, Double> getMotifsROCs(){return motifsRocs;}	
	
	public ClusterMotifFinder(ExperimentManager man, ShapeAlignConfig shapeConfig, GenomeConfig gcon) throws FileNotFoundException, UnsupportedEncodingException{
		manager = man;
		shapeAlignConfig = shapeConfig;
		K = shapeConfig.getK();		
		gconfig = gcon;
		seqgen = new SequenceGenerator<Region>(gcon.getGenome());
		shapeAlign = new ShapeAlignment(shapeAlignConfig, gconfig, manager); 
		shapeAlign.excute();
		strandedRegions = shapeAlign.getStrandedRegions();
	}
	
	// get cluster assignment and MEME
	public void execute() throws IOException{
		List<Integer> clusAssignment = cluster();
		System.out.println("Clustering finished");
		setRandomRegs();
		System.out.println("Loading randome regions finished");
		runMEME(clusAssignment);
		System.out.println("number of motifs found "+motifs.size());
	}
	
	public void setRandomRegs(){
		List<Region> randomRegions = MemeER.randomRegionPick(gconfig.getGenome(), null, numRandRegions,shapeAlignConfig.getWindowSize());
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = seqgen.execute(randomRegions.get(i));
		}
	}
	
	/**
	 * The method should be executed after initiating the class object
	 * @throws IOException 
	 */
	public List<Integer> cluster() throws IOException {
		
		// Initialize
		Vector<VectorClusterElement> bestClusterMeans=null;
		Double SSD = Double.MAX_VALUE;
		ClusterRepresentative<VectorClusterElement> crep = new Mean();
		double thresholdDistances[][] = shapeAlign.thresholdDistance();
		
		for (int i= 0 ; i<10; i++){
			for (int j=0; j<10 ; j++){
				System.out.print(thresholdDistances[i][j]);
			}
			System.out.println();
		}
		
		// converting double[] to vector
		for (int i = 0; i <thresholdDistances.length; i++){
			DefaultVectorClusterElement v = new DefaultVectorClusterElement(thresholdDistances[i]);
			profiles.add(v);
		}
		
		for(int i=0; i<=TRAIN_REPEATS; i++){			
			//Random starts
			Random generator = new Random();
			ArrayList<Integer> pNum = new ArrayList<Integer>();
			Vector<VectorClusterElement> starts = new Vector<VectorClusterElement>();
			int s=0;
			do{
				int r = generator.nextInt(profiles.size());
				if (!pNum.contains(r)){	// make sure that random number picked is not previously included
					starts.add(profiles.get(r));
					pNum.add(r);
					s++;
				}
			}while (s<K);		
			System.out.println("iteration round "+i);
		
			//Initialize clustering
			KMeansClustering<VectorClusterElement> kmc = new KMeansClustering<VectorClusterElement>(metric, crep, starts);
			kmc.setIterations(MAX_ITER);			
			// Cluster !
			kmc.clusterElements(profiles,0.01);
			
			Double currSSD = kmc.sumOfSquaredDistance();
			if(currSSD<SSD){
				SSD = currSSD;
				bestClusterMeans = kmc.getClusterMeans();
			}System.out.println("SSD: "+currSSD);	
		}
		
		return clusterList(bestClusterMeans);
	}
	
	private int getClusterAssignment(VectorClusterElement pfl, Vector<VectorClusterElement> clusMeans){
		int minCluster = -1;
        double minDist = 0.0;
        for(int i = 0; i < clusMeans.size(); i++) { 
            double clustDist = metric.evaluate(pfl, clusMeans.get(i));
            if(minCluster == -1 || clustDist < minDist) { 
                minDist = clustDist;
                minCluster = i;
            }
        }
        return minCluster;
	}
	
	private List<Integer> clusterList(Vector<VectorClusterElement> clusMeans) throws IOException{
		List<Integer> clusterAssignment = new ArrayList<Integer>();

		for(int p=0; p<profiles.size(); p++){
			int memebership = getClusterAssignment(profiles.get(p),clusMeans);
			clusterAssignment.add(memebership);
		}
		return clusterAssignment;
	}
	
	
	public void runMEME(List<Integer> clusterAssignment) throws FileNotFoundException{
		
		List<ArrayList<String>> seqClusterList = new ArrayList<ArrayList<String>>();
		// get regions and sequences from clusters
		for (int k = 0; k <K;k++){
			ArrayList<StrandedRegion> clusterRegions = new ArrayList<StrandedRegion>();
			ArrayList<String> clusterSeq = new ArrayList<String>();			
			for (int i = 0 ; i < clusterAssignment.size(); i++){
				if (clusterAssignment.get(i) == k){
					clusterRegions.add(strandedRegions.get(i));
					clusterSeq.add(seqgen.execute(strandedRegions.get(i)));
				}
			}
			regionClusterList.add(new HashSet<StrandedRegion>(clusterRegions));
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
			// writing clustered regions
			PrintWriter region_out = new PrintWriter(basedir_meme+File.separator+"clust"+(i+1)+".regions");			
			for (Region reg : regionClusterList.get(i))
				region_out.println(reg.getLocationString());		
			region_out.close();	
			
			System.out.println("output directory is "+meme_outFile);
			System.out.println("Running meme now for cluster "+(i+1));
			System.out.println("Number of sequences in cluster is "+seqcluster.size());
			Pair<List<WeightMatrix>,List<WeightMatrix>> matrices = meme.execute(seqcluster, meme_outFile, false);
			System.out.println("Finished running meme for cluster "+(i+1)+", Now evaluating motif significance");
			List<WeightMatrix> wm = matrices.car();
			List<WeightMatrix> fm = matrices.cdr();
			// Index for all the selected motifs
			int motInd = 0;
			
			if(wm.size()>0){
				//Evaluate the significance of the  motifs
				double rocScores[] = meme.motifROCScores(wm,seqcluster,randomSequences);
				System.out.println("MEME results for cluster "+(i+1));
				for(int w=0; w<fm.size(); w++){
					if(fm.get(w)!=null){
						System.out.println("\t"+fm.get(w).getName()+"\t"+ WeightMatrix.getConsensus(fm.get(w))+"\tROC:"+String.format("%.2f",rocScores[w]));
					}
					if(rocScores[w] > MemeER.MOTIF_MIN_ROC){
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
	
	public void getMotifLocations(){
		
	}
	
	public static void main(String[] args) throws IOException {
		ShapeAlignConfig shapeAlignConf = new ShapeAlignConfig(args);	
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);	
		econf.setPerBaseReadFiltering(false);
		ExperimentManager manager = new ExperimentManager(econf);
//		ShapeAlignment profile = new ShapeAlignment(shapeAlignConf, gconf, manager); 	
		ClusterMotifFinder motifFinder = new ClusterMotifFinder(manager,shapeAlignConf,gconf);
		motifFinder.execute();
				
		manager.close();		
	}
}
