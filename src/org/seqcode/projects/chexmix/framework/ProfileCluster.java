package org.seqcode.projects.chexmix.framework;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Vector;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.gseutils.Pair;
import org.seqcode.ml.clustering.ClusterRepresentative;
import org.seqcode.ml.clustering.PairwiseElementMetric;
import org.seqcode.ml.clustering.kmeans.KMeansClustering;
import org.seqcode.ml.clustering.vectorcluster.DefaultVectorClusterElement;
import org.seqcode.ml.clustering.vectorcluster.EuclideanDistance;
import org.seqcode.ml.clustering.vectorcluster.Mean;
import org.seqcode.ml.clustering.vectorcluster.VectorClusterElement;
import org.seqcode.projects.chexmix.mixturemodel.BindingSubComponents;

/**
 * ProfileCluster: produces PCC pairwise similarity scores and performs K-means clustering on the similarity matrix.
 * 
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class ProfileCluster {
	
	protected XOGPSConfig config;
	protected ExperimentManager manager;	
	private PairwiseElementMetric<VectorClusterElement> metric = new EuclideanDistance<VectorClusterElement>();
//	private PairwiseElementMetric<VectorClusterElement> metric = new ManhattanDistance<VectorClusterElement>();
	protected int profileWidth;
	
	public ProfileCluster(XOGPSConfig c, ExperimentManager man){
		config=c;
		manager=man;
		profileWidth = config.MAX_BINDINGMODEL_WIDTH;
	}
	
	public ProfileCluster(XOGPSConfig c){
		config=c;
		profileWidth = config.MAX_BINDINGMODEL_WIDTH;
	}
	
	public List<Pair<double[][][], Integer>> execute(ExperimentCondition cond, List<BindingSubComponents> typeComps, boolean doClustering) throws IOException{
		int numReps = cond.getReplicates().size();
		int numComps=typeComps.size();
		PCC pcc = new PCC(config, numReps);
		// store pairwise similarity scores
		double[][] simScores = new double[numComps][numComps];
		double[][] distanceMatrix = new double[numComps][numComps];
		for (int i=0; i< numComps; i++){
			for (int j=0; j< numComps; j++){
				simScores[i][j]=1;
				distanceMatrix[i][j]=0;
			}
		}
				
		// store all profiles to array
		double[][][] all_profile_plus = new double[numComps][numReps][profileWidth];
		double[][][] all_profile_minus = new double[numComps][numReps][profileWidth];
		for (int i=0; i< numComps;i++){
			for (int r=0; r < numReps; r++){
				for (int w=0; w< profileWidth; w++){
					all_profile_plus[i][r][w]=0;
					all_profile_minus[i][r][w]=0;
				}}}		
		
		// store component indexes
		Map<BindingSubComponents, Integer> comps2index = new HashMap<BindingSubComponents,Integer>();
		for (int i=0; i< numComps;i++){
			BindingSubComponents currComp = typeComps.get(i);
			comps2index.put(currComp, i);
			for (ControlledExperiment rep : cond.getReplicates()){
				int x = rep.getIndex();
				for (int w=0; w< typeComps.get(i).getReadProfile_plus(x).length; w++){
					all_profile_plus[i][x][w]=currComp.getReadProfile_plus(x)[w];
					all_profile_minus[i][x][w]=currComp.getReadProfile_minus(x)[w];
				}}}
		
		// Grouped components based on clustering or list of comps if not clustering
		List<List<BindingSubComponents>> clustComps = new ArrayList<List<BindingSubComponents>>();		
		if (doClustering){
			// Fill in similarity matrix
			for (int i=0; i< numComps;i++){
				for (int j=i+1; j< numComps; j++){
					pcc.execute(all_profile_plus[i], all_profile_minus[i], all_profile_plus[j], all_profile_minus[j]);
					simScores[i][j] = pcc.getCorrelation(); 
					simScores[j][i] = pcc.getCorrelation();
				}}			
			// Convert similarities to distance
			distanceMatrix = similarity2DistanceMatrix(simScores, config.CORR_THRES);		
			// Testing only
			print2dmatrix(distanceMatrix);
			
			// Perform K-means
			List<Integer> clusAssignment=null;
		
		
		try {
			clusAssignment = cluster(distanceMatrix);  //clusAssignment stores a list of integer indicating a membership
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
			//K medoid clustering	
//			clusAssignment = doClustering(distanceMatrix);
		
			// group components based on cluster assignment 
			for (int k=0; k < config.KMEANS_K; k++)
				clustComps.add(new ArrayList<BindingSubComponents>());
			for (int j=0; j< clusAssignment.size(); j++)
				clustComps.get(clusAssignment.get(j)).add(typeComps.get(j));
		}else{
			clustComps.add(typeComps);
		}
		
		List<Pair<double[][][], Integer>> BMupdateProfiles= new ArrayList<Pair<double[][][], Integer>>();		
		int cnum=0;
		for (List<BindingSubComponents> currClustComps : clustComps){
			if (currClustComps.size() > config.getMinComponentsForBMUpdate()){			
				//Sort by responsibilities
				Collections.sort(currClustComps, new Comparator<BindingSubComponents>(){
					public int compare(BindingSubComponents o1, BindingSubComponents o2) {return o1.compareByResp(o2);}
				});
				Collections.reverse(currClustComps);
				
				int indexA = comps2index.get(currClustComps.get(0));
				int indexB = comps2index.get(currClustComps.get(1));
				pcc.execute(all_profile_plus[indexA],all_profile_minus[indexA], all_profile_plus[indexB], all_profile_minus[indexB]);
				
				for (int j=2; j< currClustComps.size(); j++){
					/**
					double[] curr_profile_plus = new double[config.MAX_BINDINGMODEL_WIDTH];
					double[] curr_profile_minus = new double[config.MAX_BINDINGMODEL_WIDTH];
					for (int i=0; i< config.MAX_BINDINGMODEL_WIDTH; i++){
						curr_profile_plus[i]=mergedProfile[0][i];
						curr_profile_minus[i]=mergedProfile[1][i];
					} does this work below?
					**/
					pcc.execute(pcc.getMergedWatsonProfile(),pcc.getMergedCrickProfile(), all_profile_plus[comps2index.get(currClustComps.get(j))], all_profile_minus[comps2index.get(currClustComps.get(j))]);					
				}
				//copy
				double[][][] mergedProfile = new double [2][numReps][profileWidth];
				for (int r=0; r < numReps; r++){
					for (int w=0; w< profileWidth; w++){
						mergedProfile[0][r][w] = pcc.getMergedWatsonProfile()[r][w];
						mergedProfile[1][r][w] = pcc.getMergedCrickProfile()[r][w];
					}
				}
				
				BMupdateProfiles.add(new Pair<double[][][], Integer>(mergedProfile, currClustComps.size()));
				
				//Print coords associated with cluster: Testing only
				printClusterPoints(currClustComps, cnum);				
				cnum++; // increment cluster group number
				
			}else{
				System.err.println("The "+cond.getName()+" read distributions cannot be updated due to too few binding components ("+currClustComps.size()+"<"+config.getMinComponentsForBMUpdate()+").");			}
		}		
		return BMupdateProfiles;
	}
	
	// Convert similarity matrix to distance matrix with a given threshold
	public double[][] similarity2DistanceMatrix(double[][] inputScores, double perc){
		double[][] distanceMatrix = new double[inputScores.length][inputScores[0].length];
		double maxVal = -1;
		for (int i = 0 ; i < inputScores.length; i++)
			for (int j = 0 ; j < inputScores[i].length; j++)
				if (inputScores[i][j] > maxVal)
					maxVal = inputScores[i][j];
		for (int i = 0 ; i < inputScores.length; i++)
			for (int j = 0 ; j < inputScores[i].length; j++)
				distanceMatrix[i][j] = maxVal-inputScores[i][j];
			
		double[] data = new double[distanceMatrix.length*distanceMatrix[0].length];
		int index=0;
		for (int i=0; i< distanceMatrix.length; i++){
			for (int j=0; j< distanceMatrix[0].length; j++){
				data[index] = distanceMatrix[i][j];
				index++;
			}
		}		
		Percentile distr = new Percentile();
		double percVal = distr.evaluate(data, perc);		
		double half = distr.evaluate(data,50);
		double twentyFive = distr.evaluate(data,25);
		double fifteen = distr.evaluate(data,15);
		double ten = distr.evaluate(data,10);
		double five = distr.evaluate(data,5);
		
		/**
		System.out.println(perc+ "% value is "+percVal);
		System.out.println("50% value is "+half);
		System.out.println("25% value is "+twentyFive);
		System.out.println("15% value is "+fifteen);
		System.out.println("10% value is "+ten);
		System.out.println("5% value is "+five);
		**/
		
		double[][] tdist = new double[inputScores.length][inputScores[0].length];
		for (int i=0; i< distanceMatrix.length; i++){
			for (int j=0; j< distanceMatrix[0].length; j++){
				if (distanceMatrix[i][j] <= percVal)
					tdist[i][j] = distanceMatrix[i][j];
				else
					tdist[i][j] = maxVal;			
			}
		}
		return tdist;
	}
	/**
	 * Prints matrix to a file for testing
	 * @throws IOException 
	 */
	public void print2dmatrix(double[][] matrix) throws IOException{
		String filename=config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_cluster.txt";
		FileWriter fout = new FileWriter(filename);
		for (int i=0; i < matrix.length ; i++){
			for (int j=0; j < matrix[i].length; j++)
				fout.write(matrix[i][j]+"\t");
			fout.write("\n");
		}
		fout.close();
	}
	
	/**
	 * Prints cluster coordinates
	 * @throws IOException 
	 */	
	public void printClusterPoints(List<BindingSubComponents> clusterComps, int num) throws IOException{
		String filename=config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_cluster_"+num+".points";
		FileWriter fout = new FileWriter(filename);
		for (BindingSubComponents comp : clusterComps)
			fout.write(comp.getCoord().toString()+"\n");
		fout.close();
	}
	
	// K-medoid clustering
	public List<Integer> doClustering(double[][] distanceMatrix){
		int numData = distanceMatrix.length;
		int numClusters = config.KMEANS_K;
		List<Integer> clusterAssignment = new ArrayList<Integer>();
		List<Integer> bestClusterMedoids=null;
		Double SD = Double.MAX_VALUE;
		for(int j=0; j<=config.KMEANS_TRAIN_REPEATS; j++){	
			List<Integer> clusterMedoids = kMedoidClustering(distanceMatrix);
			double totalDist =0.0;
			for (int e=0; e < numData; e++){
				int minCluster =-1;
				double minDist =0.0;
				for (int i=0; i < numClusters; i++){
					double clustDist = distanceMatrix[e][clusterMedoids.get(i)]; 
					if (minCluster ==-1 || clustDist < minDist){
						minDist = clustDist;
						minCluster = i;
					}
				}totalDist += minDist;	
			}
			if (totalDist < SD){
				SD = totalDist;
				bestClusterMedoids = clusterMedoids;
			}
		}
		System.out.println("total distance "+SD+" medoids "+bestClusterMedoids.toString());
		for (int e=0; e < numData; e++){
			int membership = getClusterAssignment(distanceMatrix, e, bestClusterMedoids);
			clusterAssignment.add(membership);
		}
		return clusterAssignment;
	}
		
		
	public List<Integer> kMedoidClustering(double[][] distanceMatrix){			
		int numData = distanceMatrix.length;
		int numClusters = config.KMEANS_K;		
		// Select random points as starting medoids
		Random rand = new Random();
		ArrayList<Integer> clusterMedoids = new ArrayList<Integer>();
		int s=0;
		do{ 
			int r = rand.nextInt(numData);
			if (!clusterMedoids.contains(r)){
				clusterMedoids.add(r);
				s++;
			}
		} while ( s < numClusters);
		
		// Assign data points to medoids	
		List<List<Integer>> clusters = new ArrayList<List<Integer>>();
		for (int i=0; i < numClusters; i++)
			clusters.add(new ArrayList<Integer>());	
		boolean converged=false;
		for(int j = 0; j < config.KMEANS_MAX_ITER && !converged; j++) {
			for (int i=0; i < numClusters; i++)
				clusters.get(i).clear();				
			List<Integer> oldClusterMedoids = new ArrayList<>(clusterMedoids);	
			// K-medoids 
			//assign datapoints to clusters
			for (int e=0; e < numData; e++){
				int clustAssignment = getClusterAssignment(distanceMatrix, e, clusterMedoids);
				clusters.get(clustAssignment).add(e);
			}			
			// get cluster medoids
			clusterMedoids.clear();
			for (int i=0; i < numClusters; i++){
				int newMedoid =-1;
				double clustMinDist = 0.0;
				for (Integer m : clusters.get(i)){ // Candidate medoids
					double clustTotalDist =0.0;
					for (Integer e : clusters.get(i))
						clustTotalDist += distanceMatrix[e][m];
					if (newMedoid ==-1 || clustTotalDist < clustMinDist){
						clustMinDist = clustTotalDist;
						newMedoid = m;
					}
				}
				clusterMedoids.add(newMedoid);
			}
			
			// Check convergence
			Collections.sort(oldClusterMedoids);
			Collections.sort(clusterMedoids);
			if (oldClusterMedoids.equals(clusterMedoids))
				converged=true;
		}
		return clusterMedoids;
	}
	
	private int getClusterAssignment(double[][] distanceMatrix, int elm, List<Integer> clusterMedoids){
		int minCluster = -1;
		double minDist = 0.0;
		for (int i=0; i < clusterMedoids.size(); i++){
			double clustDist = distanceMatrix[elm][clusterMedoids.get(i)];
			if (minCluster == -1 || clustDist <minDist){
				minDist = clustDist;
				minCluster = i;
			}
		}
		return minCluster;
	}	
	
	/**
	 * The method should be executed after initiating the class object
	 * @throws IOException 
	 */
	public List<Integer> cluster(double[][] thresholdDistances) throws IOException {
		
		ArrayList<VectorClusterElement> profiles = new ArrayList<VectorClusterElement>();
		
		// Initialize
		Vector<VectorClusterElement> bestClusterMeans=null;
		Double SSD = Double.MAX_VALUE;
		ClusterRepresentative<VectorClusterElement> crep = new Mean();
		
		// converting double[] to vector
		for (int i = 0; i <thresholdDistances.length; i++){
			DefaultVectorClusterElement v = new DefaultVectorClusterElement(thresholdDistances[i]);
			profiles.add(v);
		}
		
		for(int i=0; i<=config.KMEANS_TRAIN_REPEATS; i++){			
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
			}while (s<config.KMEANS_K);		
			//System.out.println("iteration round "+i);
		
			//Initialize clustering
			KMeansClustering<VectorClusterElement> kmc = new KMeansClustering<VectorClusterElement>(metric, crep, starts);
//			KMedoidsClustering<VectorClusterElement> kmc = new KMedoidsClustering<VectorClusterElement>(metric, starts);
			kmc.setIterations(config.KMEANS_MAX_ITER);			
			// Cluster !
			kmc.clusterElements(profiles,config.KMEANS_CONVERGENCE_THRES);
			
			Double currSSD = kmc.sumOfSquaredDistance();
//			Double currSSD = kmc.sumOfDistance();
			if(currSSD<SSD){
				SSD = currSSD;
				bestClusterMeans = kmc.getClusterMeans();
//				bestClusterMeans = kmc.getClusterMedoids();
			}//System.out.println("SSD: "+currSSD);	
		}
		
		return clusterList(bestClusterMeans, profiles);
	}
	
	private List<Integer> clusterList(Vector<VectorClusterElement> clusMeans, ArrayList<VectorClusterElement> profiles) throws IOException{
		List<Integer> clusterAssignment = new ArrayList<Integer>();

		for(int p=0; p<profiles.size(); p++){
			int memebership = getClusterAssignment(profiles.get(p),clusMeans);
			clusterAssignment.add(memebership);
		}
		return clusterAssignment;
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


	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		XOGPSConfig config = new XOGPSConfig(gcon, args);
		
		double[][] array = new double[6][6];
		for (int i=0; i < 6; i++)
			array[i][i]=0.0;
		array[0][1]=1;array[1][0]=1;
		array[0][2]=2;array[2][0]=2;
		array[0][3]=4;array[3][0]=4;
		array[0][4]=3;array[4][0]=3;
		array[0][5]=5;array[5][0]=5;
		array[1][2]=1;array[2][1]=1;
		array[1][3]=2;array[3][1]=2;
		array[1][4]=2;array[4][1]=2;
		array[1][5]=3;array[5][1]=3;
		array[2][3]=1;array[3][2]=1;
		array[2][4]=1;array[4][2]=1;
		array[2][5]=2;array[5][2]=2;
		array[3][4]=1;array[4][3]=1;
		array[3][5]=2;array[5][3]=2;
		array[4][5]=1;array[5][4]=1;
		
		ProfileCluster cluster = new ProfileCluster(config);
		List<Integer> clusAssignment=cluster.doClustering(array);
		for (int i : clusAssignment)
			System.out.println(i+"\t");
	
	}
}
