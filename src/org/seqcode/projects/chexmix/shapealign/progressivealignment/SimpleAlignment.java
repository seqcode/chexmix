package org.seqcode.projects.chexmix.shapealign.progressivealignment;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeNode;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;
import org.seqcode.projects.chexmix.utilities.ArrayListSorter;

/**
 * SimpleAlignment progressively merges profiles from the highest to the lowest occupancy.
 * input : Stranded points
 * output : Stranded points with aligned profiles
 * @author naomi
 */
public class SimpleAlignment {
	protected ExperimentManager exptMan;	
	protected List<StrandedPoint> points;
	protected List<StrandedPoint> alignedPoints;
	protected List<Map<StrandedPoint, Integer>> condPoints2index;
	protected int win;
	protected int pccWindow =200;
	protected int numPoints; 
	protected String filename;
	protected CompositeTagDistribution maker;
	protected DefaultMutableTreeNode root;	
	static final double MINIMUM_VALUE = -10000;
	static final int DIAG = 1;
	static final int LEFT = 2;
	static final int UP = 4;
	
	public SimpleAlignment(List<StrandedPoint> points, ExperimentManager eMan, int win, String filename){
		this.points = points;
		exptMan = eMan;
		this.win = win;
		this.filename = filename;
		numPoints = points.size();
		maker = new CompositeTagDistribution(points, exptMan, win, true);
		for(ExperimentCondition cond : exptMan.getConditions()){
			String compositeFileName = this.filename+"_composite."+cond.getName()+".txt";
			maker.printProbsToFile(cond, compositeFileName);
		}
		alignedPoints = new ArrayList<StrandedPoint>();
		condPoints2index = new ArrayList<Map<StrandedPoint, Integer>>();
	}
	
	public List<List<ArrayListSorter>> sortPointsByTags(){				
		// Sort per location tags based on occupancies		
		List<List<ArrayListSorter>> sortedlist = new ArrayList<List<ArrayListSorter>>();
		for(ExperimentCondition cond : exptMan.getConditions()){
			List<ArrayListSorter> occupancies = new ArrayList<ArrayListSorter>();
			Map<StrandedPoint, Integer> points2index = new HashMap<StrandedPoint, Integer>();
			//Iterate through points
			int p=0;
			for (StrandedPoint pt : points){
				ArrayListSorter elem = new ArrayListSorter();
				elem.setKey(pt);
				double[] watsontags = maker.getPointWatson(pt,cond);
				double[] cricktags = maker.getPointCrick(pt, cond);
				float sum=0;
				for (int i=0; i<watsontags.length; i++)
					sum =+ (float) watsontags[i];
				for (int i=0; i<cricktags.length;i++)
					sum += (float) cricktags[i];
				elem.setVal(sum);
				occupancies.add(elem);
				points2index.put(pt, p);
				p++;
			}
			Collections.sort(occupancies);
			sortedlist.add(occupancies);
			condPoints2index.add(points2index);
		}
		return sortedlist;
	}

	public void execute(boolean useNW){
		
		if (useNW)
			System.out.println("executing NW");
		else
			System.out.println("executing PCC");
		
		List<List<ArrayListSorter>> sortedlist = sortPointsByTags();
		for(ExperimentCondition cond : exptMan.getConditions()){			
			// Index of current profile
			int currIndex=numPoints;			
			for (int pt=1; pt<numPoints;pt++){
				StrandedPoint firstp = null;
				if (pt==1)
					firstp = (StrandedPoint) sortedlist.get(cond.getIndex()).get(0).getKey();
				StrandedPoint currpt = (StrandedPoint) sortedlist.get(cond.getIndex()).get(pt).getKey();
				if (useNW)
					AddNodeNW(cond, firstp , currpt, currIndex);
				else
					AddNodePCC(cond, firstp , currpt, currIndex);
				currIndex++;
			}						
		}
		
		System.out.println("child count: "+root.getChildCount()+"\tdepth: "+root.getDepth()+"\tleaf count: "+root.getLeafCount());
		
		for(ExperimentCondition cond : exptMan.getConditions()){	
			
			// Walk through a tree in pre-order
			Enumeration<DefaultMutableTreeNode> en = root.preorderEnumeration();
			while (en.hasMoreElements()){
				DefaultMutableTreeNode node = en.nextElement();
			
				Profile currprofile = (Profile) node.getUserObject();
			
				TreeNode[] path = node.getPath();
				System.out.println((node.isLeaf() ? " -" : "+ ") + path[path.length-1]+"\t"+currprofile.getIndex());
			
				if (!currprofile.isProfile()){	
					
					StrandedPoint origp=null;
					for (StrandedPoint pt : condPoints2index.get(cond.getIndex()).keySet())
						if(condPoints2index.get(cond.getIndex()).get(pt)==currprofile.getIndex())
							origp=pt;
					
					System.out.println(origp.toString());
					
					char strand;
					int location = 0;					
					if (currprofile.isReverse()){
						System.out.println("is reverse");
						strand = origp.getStrand()=='+' ? '-' : '+';
						location = origp.getLocation()-currprofile.getOffset();
					}else{
						System.out.println("is not reverse");
						strand = origp.getStrand();
						if (strand=='-')
							location = origp.getLocation()-currprofile.getOffset();
						else
							location = origp.getLocation()+currprofile.getOffset();
					}				
					alignedPoints.add(new StrandedPoint(origp.getGenome(), origp.getChrom(), location, strand));				
				}
			}
		}
		printAlingedPointsToFile();
	}
	
	public void printAlingedPointsToFile(){
		String compositeFileName = filename+"_alinged.points";
		try{
			FileWriter fout = new FileWriter(compositeFileName);		
			for (StrandedPoint pt : alignedPoints){
				fout.write(pt.toString()+"\n");
			}	
			fout.close();		
		} catch (IOException e) {
			e.printStackTrace();
		}			
	}
	
	public void AddNodeNW(ExperimentCondition cond, StrandedPoint firstp, StrandedPoint currpt, int currIndex){
		
		double[] watsonA = null;
		double[] crickA = null;	
		int profileIndex = 0;
		if (firstp != null){
			// Get number of tags from initial regions to compare
			watsonA = maker.getPointWatson(firstp, cond);
			crickA = maker.getPointCrick(firstp, cond);
		}else{
			Profile prevprofile = (Profile) root.getUserObject();
			watsonA = prevprofile.getWatsonProfile(); 
			crickA = prevprofile.getCrickProfile();
			profileIndex = prevprofile.getIndex();
		}
		double[] watsonB = maker.getPointWatson(currpt, cond);
		double[] crickB = maker.getPointCrick(currpt, cond);
		
		//Normalize A and B tags by max points
		double maxA=0, maxB=0;
		for (int i=0; i<win; i++){
			if (watsonA[i] > maxA)
				maxA = watsonA[i];
			if (crickA[i] > maxA)
				maxA = crickA[i];
			if (watsonB[i] > maxB)
				maxB = watsonB[i];
			if (crickB[i] > maxB)
				maxB = crickB[i];
		}
		for (int i=0; i< win; i++){
			watsonA[i] /= maxA;
			crickA[i] /= maxA;
			watsonB[i] /= maxB;
			crickB[i] /= maxB;
		}			
		
		/**
		//Normalize A and B tags
		double sumA=0, sumB=0;
		for (int i=0; i<win; i++){
			sumA += watsonA[i];
			sumA += crickA[i];
			sumB += watsonB[i];
			sumB += crickB[i];
		}
		for (int i=0; i< win; i++){
			watsonA[i] /= sumA;
			crickA[i] /= sumA;
			watsonB[i] /= sumB;
			crickB[i] /= sumB;
		}	
		**/
		
		//Reverse B tags
		double[] rWatsonB = new double[win];
		double[] rCrickB = new double[win];
		for (int i=0; i< win; i++){
			rWatsonB[win-i-1] = crickB[i];
			rCrickB[win-i-1] = watsonB[i];
		}
			
		// align using two possible ways
		NeedlemanWunsch alignOne = new NeedlemanWunsch(watsonA, crickA, watsonB, crickB, win);
		NeedlemanWunsch alignTwo = new NeedlemanWunsch(watsonA, crickA, rWatsonB, rCrickB, win);
		alignOne.buildMatrix();
		alignTwo.buildMatrix();
			
		Stack<Integer> traceBack = new Stack<Integer>();
		double[] watson_b;
		double[] crick_b;
		double maxScore = MINIMUM_VALUE;
		boolean reverseB;
		int s_x_coord = 0;
		int s_y_coord = 0;
		int e_x_coord = 0;
		int e_y_coord = 0;		
		if (alignOne.getMaxScore() > alignTwo.getMaxScore()){	
			watson_b = watsonB;
			crick_b = crickB;
			maxScore = alignOne.getMaxScore();
			traceBack = alignOne.getTraceBack();
			s_x_coord = alignOne.getStartX();
			s_y_coord = alignOne.getStartY();
			e_x_coord = alignOne.getEndX();
			e_y_coord = alignOne.getEndY();		
			reverseB = false;		
		}else{	
			watson_b = rWatsonB;
			crick_b = rCrickB;
			maxScore = alignTwo.getMaxScore();
			traceBack = alignTwo.getTraceBack();
			s_x_coord = alignTwo.getStartX();
			s_y_coord = alignTwo.getStartY();
			e_x_coord = alignTwo.getEndX();
			e_y_coord = alignTwo.getEndY();	
			reverseB = true;
		}		
		double x_mid = (s_x_coord + e_x_coord)/2;
		double y_mid = (s_y_coord + e_y_coord)/2;		
		int offset = (int) (y_mid - x_mid);
		System.out.println("alignment start coordinates "+ s_x_coord + " : " + s_y_coord);
		System.out.println("alignment end coordinates "+ e_x_coord + " : " + e_y_coord);		
		System.out.println(Arrays.toString(traceBack.toArray()));
		System.out.println("offset is "+offset);
		System.out.println("reverse ? "+ reverseB);
			
		// mereged profiles
		double[] w_profile = new double[win+1];
		double[] c_profile = new double[win+1];
		for (int i=0; i < win; i++){
			w_profile[i] = 0;
			c_profile[i] = 0;
		}
			
		int current_x = e_x_coord-1;
		int current_y = e_y_coord-1;
			
		// trace back for x : profile A
		@SuppressWarnings("unchecked")
		Stack<Integer> xTraceBack = (Stack<Integer>) traceBack.clone();		
		for (int i = e_x_coord; i >= s_x_coord ; i--){				
			if (current_x >= 0){
				w_profile[i] += watsonA[current_x];
				c_profile[i] += crickA[current_x];			
				if ( !xTraceBack.empty() ){			
					if (xTraceBack.peek() == DIAG || xTraceBack.peek() == LEFT)
						current_x --;			
					xTraceBack.pop();
				}	
			}
		}				
		// trace back for y : profile B 
		@SuppressWarnings("unchecked")
		Stack<Integer> yTraceBack = (Stack<Integer>) traceBack.clone();;			
		for (int i = e_y_coord; i >= s_y_coord ; i--){				
			if (current_y >= 0){
				w_profile[i] += watson_b[current_y];
				c_profile[i] += crick_b[current_y];
				if ( !yTraceBack.empty() ){				
					if (yTraceBack.peek() == DIAG || yTraceBack.peek() == UP)
						current_y --;	
					yTraceBack.pop();			
				}
			}
		}
			
		// Look up index for the points to track profiles/tags
		int indexB = condPoints2index.get(cond.getIndex()).get(currpt);			
		// Store tags as profiles
		Profile profileB = new Profile(watsonB, crickB, indexB);
			
		// set profile B offset relative to profile A
		profileB.setOffset(offset);
		profileB.setReverse(reverseB);
		
		System.out.println("max score from NW "+maxScore);
//		System.out.println("offset "+offset);
		System.out.println("merged profiles, watson");
		for (int i=0; i< w_profile.length;i++)
			System.out.print(w_profile[i]+",");
		System.out.println("\nmerged profiles, crick");
		for (int i=0; i< c_profile.length;i++)
			System.out.print(c_profile[i]+",");
		System.out.println();
						
		// Create merged profile
		Profile profile = new Profile(w_profile, c_profile, currIndex);
		profile.setProfile(true);
		
		// if this is the first comparison
		if (firstp != null){
			// Look up index for the points to track profiles/tags
			int indexA = condPoints2index.get(cond.getIndex()).get(firstp);
			// Store tags as profiles
			Profile profileA = new Profile(watsonA, crickA, indexA);
			// Store index
			profile.setPairProfile(indexA, indexB, maxScore);
			
			// Construct an initial tree
			DefaultMutableTreeNode initNode = new DefaultMutableTreeNode(profile);
			initNode.add(new DefaultMutableTreeNode(profileA));
			initNode.add(new DefaultMutableTreeNode(profileB));		
			root = initNode;
		}else{
			profile.setPairProfile(profileIndex, indexB, maxScore);
			DefaultMutableTreeNode currnode = new DefaultMutableTreeNode(profile);
			currnode.add(root);		
			currnode.add(new DefaultMutableTreeNode(profileB));
			// move a reference to root
			root=currnode;			
		}		
	}
		
	public void AddNodePCC(ExperimentCondition cond, StrandedPoint firstp, StrandedPoint currpt, int currIndex){
		//if there is previous profile added
		if (firstp!=null){
			// Get number of tags from initial regions to compare
			double[] watsonA = maker.getPointWatson(firstp, cond);
			double[] crickA = maker.getPointCrick(firstp, cond);
			double[] watsonB = maker.getPointWatson(currpt, cond);
			double[] crickB = maker.getPointCrick(currpt, cond);
			
			// Perform Pearson correlations on the two regions
			PCC corr = new PCC(watsonA, crickA, watsonB, crickB, win, pccWindow);
			// Look up index for the points to track profiles/tags
			int indexA = condPoints2index.get(cond.getIndex()).get(firstp);
			int indexB = condPoints2index.get(cond.getIndex()).get(currpt);			
			// Store tags as profiles
			Profile profileA = new Profile(watsonA, crickA, indexA);
			Profile profileB = new Profile(watsonB, crickB, indexB);
			// set profile B offset relative to profile A
			profileB.setOffset(corr.getOffset());
			profileB.setReverse(corr.isReverse());
			
			System.out.println("correlation "+corr.getMaxScore());
			System.out.println("offset "+corr.getOffset());
			System.out.println("merged profiles, watson");
			for (int i=0; i< corr.getWatsonProfile().length;i++)
				System.out.print(corr.getWatsonProfile()[i]+",");
			System.out.println("\nmerged profiles, crick");
			for (int i=0; i< corr.getCrickProfile().length;i++)
				System.out.print(corr.getCrickProfile()[i]+",");
			System.out.println();
			
			// Create merged profile
			Profile profile = new Profile(corr.getWatsonProfile(), corr.getCrickProfile(), currIndex);
			// Store index
			profile.setPairProfile(indexA, indexB, corr.getMaxScore());
			profile.setProfile(true);

			// Construct an initial tree
			DefaultMutableTreeNode initNode = new DefaultMutableTreeNode(profile);
			initNode.add(new DefaultMutableTreeNode(profileA));
			initNode.add(new DefaultMutableTreeNode(profileB));		
			root = initNode;
			
		}else{			
			double[] currWatson = maker.getPointWatson(currpt, cond);
			double[] currCrick = maker.getPointCrick(currpt, cond);
			
			Profile prevprofile = (Profile) root.getUserObject();			
			PCC currpcc = new PCC(prevprofile.getWatsonProfile(), prevprofile.getCrickProfile(), currWatson, currCrick, win, pccWindow);
			int currptIndex = condPoints2index.get(cond.getIndex()).get(currpt);
			
			Profile currTagProfile = new Profile(currWatson, currCrick, currptIndex);
			currTagProfile.setOffset(currpcc.getOffset());
			currTagProfile.setReverse(currpcc.isReverse());
			
			System.out.println("correlation "+currpcc.getMaxScore());
			System.out.println("offset "+currpcc.getOffset());
			System.out.println("merged profiles, watson");
			for (int i=0; i< currpcc.getWatsonProfile().length;i++)
				System.out.print(currpcc.getWatsonProfile()[i]+",");
			System.out.println("\nmerged profiles, crick");
			for (int i=0; i< currpcc.getCrickProfile().length;i++)
				System.out.print(currpcc.getCrickProfile()[i]+",");
			System.out.println();
			
			// create merged profiles
			Profile mergedprofile = new Profile(currpcc.getWatsonProfile(), currpcc.getCrickProfile(),currIndex);
			// Store index
			mergedprofile.setPairProfile(prevprofile.getIndex(),currptIndex, currpcc.getMaxScore());
			mergedprofile.setProfile(true);
			
			// Add node to a tree and point to the root of tree
			DefaultMutableTreeNode currnode = new DefaultMutableTreeNode(mergedprofile);
			currnode.add(root);		
			currnode.add(new DefaultMutableTreeNode(currTagProfile));
			// move a reference to root
			root=currnode;
		}		
	}	
	
	//Main method to make new composite distributions
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
		if(args.length==0){
			System.err.println("CompositeTagDistribution:"+
					"\t--cpoints <stranded point file>"+
					"\t--cwin <window around points>"+
					"\t--out <output file name>"+
					"Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			// parse command line arguments
			int w = Args.parseInteger(args, "cwin", 600);
			String pFile = Args.parseString(args, "cpoints", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
			String filename = Args.parseString(args, "out", "out");
			boolean useNW = Args.parseFlags(args).contains("nw") ? true : false;
			
			SimpleAlignment alignment = new SimpleAlignment(pts, manager, w, filename);			
			alignment.execute(useNW);
			manager.close();
		}
	}	
}