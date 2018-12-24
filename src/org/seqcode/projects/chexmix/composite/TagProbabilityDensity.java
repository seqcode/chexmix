package org.seqcode.projects.chexmix.composite;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;


/**
 * TagProbabilityDensity defines a (probabilistic) model of sequenced tag occurrences around a binding event.
 * The probability is directional (i.e. stranded). The probability distributions can be asymmetric on Watson & Crick strands. 
 * Given a signed distance from the event, the TagProbabilityDensity should return 
 * a relative probability of seeing a read at that distance (given that the tag is on Watson or Crick).
 * 
 * TODO: support multi-rep TagDistributions
 * 
 * @author shaunmahony
 *
 */
public class TagProbabilityDensity {
	
	private static int count=0;
	
	protected final double LOG2 = Math.log(2);
	protected final double TAGDISTRIB_MIN_PROB = 1e-100; //Minimum probability allowed in any tag distribution
	protected int index;
	protected int winSize;
	protected int left, right; //relative left & right positions
	protected double[] watsonData, crickData; //Data landscape should be based on (typically tag counts)
	protected double[] watsonProbs, crickProbs; //Probability landscapes
	protected double[] watsonLogProbs, crickLogProbs; //Log probabilities
	protected int watsonSummit, crickSummit;		// relative positions of highest probs
	protected int influenceRange; //95% probability range (over both strands)
	protected double bgProb, logBgProb;
	protected boolean isGaussian=false;
	protected double gaussOffsetW=0, gaussOffsetC=0, gaussSigmaW=-1,gaussSigmaC=-1; 
	
	public TagProbabilityDensity(int size){
		index=count; count++;
		winSize=size;
		init(-(winSize/2), winSize-(winSize/2));
	}
	public TagProbabilityDensity(int l, int r){
		index=count; count++;
		init(l, r);
	}
	
	public TagProbabilityDensity(File distFile){
		index=count; count++;
		int min=Integer.MAX_VALUE, max=Integer.MIN_VALUE;
		try {
			List<Pair<Integer,Double>> empiricalDistribution = new LinkedList<Pair<Integer,Double>>(); 
			BufferedReader reader = new BufferedReader(new FileReader(distFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
		              Integer dist = new Integer(words[0]);
		              
		              if(dist.intValue() < min) 
		            	  min = dist.intValue();
		              if(dist.intValue() >max)
		            	  max = dist.intValue();
		                
		              Pair<Integer,Double> p = new Pair<Integer,Double>(dist, new Double(words[1]));
		              if (p.cdr().doubleValue()>=0)	// should be non-negative value
		            	  empiricalDistribution.add(p);
		              else {
		            	  System.err.println("\nTag distribution file contains negative probability(count) values!"); 
		            	  System.exit(1);
		              }
	            }
	        }
	        winSize = max-min+1;
	        init(min, max);
	        
	        loadData(empiricalDistribution, null);
			makeProbabilities();
			
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public TagProbabilityDensity(List<Pair<Integer,Double>> empiricalWatson, List<Pair<Integer,Double>> empiricalCrick){
		index=count; count++;
		int min=Integer.MAX_VALUE, max=Integer.MIN_VALUE;
		try{
			for(Pair<Integer,Double> p : empiricalWatson){
				if(p.car() < min)
					min = p.car();
				if(p.car() > max)
					max = p.car();
			}
			if(empiricalCrick!=null){
				for(Pair<Integer,Double> p : empiricalCrick){
					if(p.car() < min)
						min = p.car();
					if(p.car() > max)
						max = p.car();
				}
			}
			winSize = max-min+1;
	        init(min, max);
	        loadData(empiricalWatson, empiricalCrick);
			makeProbabilities();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//Accessors
	public int getIndex(){return index;}
	public void setIndex(int i){index=i;}
	public int getLeft(){return left;}
	public int getRight(){return right;}
	public int getWinSize(){return winSize;}
	public int getWatsonSummit(){return watsonSummit;}
	public int getCrickSummit(){return crickSummit;}
	public int getInfluenceRange(){return influenceRange;}
	public double[] getWatsonProbabilities(){return  watsonProbs.clone();}
	public double[] getCrickProbabilities(){return  crickProbs.clone();}
	public boolean isGaussian(){return isGaussian;}
	public double getGaussOffsetW(){return gaussOffsetW;}
	public double getGaussOffsetC(){return gaussOffsetC;}
	public double getGaussSigmaW(){return gaussSigmaW;}
	public double getGaussSigmaC(){return gaussSigmaC;}
	
	/**
	 *  Look up the probability corresponding to a distance and strand relationship to the central position.
	 *  Distance should be defined as (tag position - central position)
	 */
	public double probability(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
			return(bgProb);
		}else{
			return(watsonStrand ? 
					watsonProbs[distance-left] : crickProbs[distance-left]);
		}
	}
	/**
	 *  Look up the log probability corresponding to a distance and strand relationship to the central position.
	 *  Distance should be defined as (tag position - central position)
	 */
	public double logProbability(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
		  return(logBgProb);
		}else{
			return(watsonStrand ? 
					watsonLogProbs[distance-left] : crickLogProbs[distance-left]);
		}	  
	}
	
	/**
	 * Load data & make probabilities
	 * Input lists are pairs of relative positions and tag counts.
	 * Assumes the input lists are sorted according to increasing position.
	 * If the crick list is null, the mirror image of the watson list is used.   
	 * @param watsonTagDist: list of paired ints and doubles
	 * @param crickTagDist: list of paired ints and doubles. If null, uses mirror of watsonTagDist
	 * @throws Exception 
	 */
	public void loadData(List<Pair<Integer, Double>> watsonTagDist, List<Pair<Integer, Double>> crickTagDist) throws Exception{		
		//Find left, right values first
		for(Pair<Integer, Double> p : watsonTagDist){
			if(p.car()<left){ left=p.car();}
			if(p.car()>right){right=p.car();}
		}
		if(crickTagDist!=null){
			for(Pair<Integer, Double> p : crickTagDist){
				if(p.car()<left){ left=p.car();}
				if(p.car()>right){right=p.car();}
			}
		}
		
		//Initialize arrays
		init(left, right);
		
		//Populate the data array (assumes sorted)
		//WATSON
		int last=left-1;
		for(Pair<Integer, Double> p : watsonTagDist){
			int index = p.car();
			double val = p.cdr();
			//if list is not properly sorted, throw exception
			if(index-last<0)
				throw new Exception("Incorrectly sorted binding read distribution data!"); 
			//if unevenly spaced, linearly interpolate between values
			if(index-last>1){
				double lastVal=dataVal(last, true);
				double step = (val-lastVal)/(double)(index-last);
				for(int i=1; i<(index-last); i++){
					watsonData[(last+i)-left]=lastVal+(step*(double)i);
				}
			}
			watsonData[index-left]=val;
			last = p.car();
		}
		//CRICK
		if(crickTagDist==null){
			for(int x=0; x<watsonData.length; x++)
				crickData[crickData.length-x-1] = watsonData[x];
		}else{
			last=left-1;
			for(Pair<Integer, Double> p : crickTagDist){
				int index = p.car();
				double val = p.cdr();
				//if list is not properly sorted, throw exception
				if(index-last<0)
					throw new Exception("Incorrectly sorted binding read distribution data!"); 
				//if unevenly spaced, linearly interpolate between values
				if(index-last>1){
					double lastVal=dataVal(last, false);
					double step = (val-lastVal)/(double)(index-last);
					for(int i=1; i<(index-last); i++){
						crickData[(last+i)-left]=lastVal+(step*(double)i);
					}
				}
				crickData[index-left]=val;
				last = p.car();
			}
		}
		makeProbabilities();
	}
	
	/**
	 * Directly load arrays
	 * @param w
	 * @param c
	 * @throws Exception 
	 */
	public void loadData(double[] w, double[] c) throws Exception{
		if(w.length!=winSize || c.length!=winSize){
			throw new Exception("TagProbabilityDensity: trying to load data of unmatched size");
		}
		for(int i=0; i<w.length; i++){
			watsonData[i] = w[i];
			crickData[i] = c[i];
		}
		makeProbabilities();
	}
	
	/**
	 * Define a symmetric pair of gaussian distributions around the origin.
	 * @param offset: mean of gaussian this number of bp upstream for Watson strand and this number of bp downstream for Crick strand
	 * @param gaussianSigma 
	 */
	public void loadGaussianDistrib(double offsetWatson, double gaussianSigmaWatson, double offsetCrick, double gaussianSigmaCrick){
		NormalDistribution wNorm = new NormalDistribution(offsetWatson, gaussianSigmaWatson);
		NormalDistribution cNorm = new NormalDistribution(offsetCrick, gaussianSigmaCrick);
		for(int i=left; i<=right; i++){
			watsonData[i-left] = wNorm.density(i);
			crickData[i-left] = cNorm.density(i);
		}
		makeProbabilities();
		isGaussian = true;
		gaussOffsetW = offsetWatson;
		gaussOffsetC = offsetCrick;
		gaussSigmaW = gaussianSigmaWatson;
		gaussSigmaC = gaussianSigmaCrick;
	}
	
	/**
	 * Define a flat distribution across the full window.
	 * @param offset: mean of gaussian this number of bp upstream for Watson strand and this number of bp downstream for Crick strand
	 * @param gaussianSigma 
	 */
	public void loadFlatDistrib(){
		for(int i=left; i<=right; i++){
			watsonData[i-left] = 1;
			crickData[i-left] = 1;
		}
		makeProbabilities();
	}
	
	
	/**
	 * Set a probability landscape according to the data. 
	 */
	protected void makeProbabilities(){
		double totalW=0, totalC=0, minProb=Double.MAX_VALUE;
		for(int i=left; i<=right; i++){
			totalW+=Math.max(TAGDISTRIB_MIN_PROB,dataVal(i, true));
			totalC+=Math.max(TAGDISTRIB_MIN_PROB,dataVal(i, false));
		}
		for(int i=left; i<=right; i++){
			watsonProbs[i-left] = Math.max(TAGDISTRIB_MIN_PROB,dataVal(i, true))/totalW;
			crickProbs[i-left] = Math.max(TAGDISTRIB_MIN_PROB,dataVal(i, false))/totalC;
			watsonLogProbs[i-left] = Math.log(watsonProbs[i-left])/LOG2;
			crickLogProbs[i-left] = Math.log(crickProbs[i-left])/LOG2;
			if(watsonProbs[i-left]<minProb){minProb = watsonProbs[i-left];}
			if(crickProbs[i-left]<minProb){minProb = crickProbs[i-left];}
		}
		Pair<Double, TreeSet<Integer>> wSorted = StatUtil.findMax(watsonProbs);
		Pair<Double, TreeSet<Integer>> cSorted = StatUtil.findMax(crickProbs);
		watsonSummit = wSorted.cdr().first()+left;
		crickSummit = cSorted.cdr().first()+left;
		
		bgProb = Math.max(TAGDISTRIB_MIN_PROB, minProb/100);
		logBgProb = Math.log(bgProb)/LOG2;

		updateInfluenceRange();
	}
	
	protected void smooth(int splineStepSize, int avgStepSize){
		watsonProbs=StatUtil.cubicSpline(watsonProbs, splineStepSize, avgStepSize);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(watsonProbs);
		watsonSummit = sorted.cdr().first()+left;
		crickProbs=StatUtil.cubicSpline(crickProbs, splineStepSize, avgStepSize);
		sorted = StatUtil.findMax(crickProbs);
		crickSummit = sorted.cdr().first()+left;
	}
	protected void smoothGaussian (int kernelWidth){
		watsonProbs=StatUtil.gaussianSmoother(watsonProbs, kernelWidth);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(watsonProbs);
		watsonSummit = sorted.cdr().first()+left;
		crickProbs=StatUtil.gaussianSmoother(crickProbs, kernelWidth);
		sorted = StatUtil.findMax(crickProbs);
		crickSummit = sorted.cdr().first()+left;
	}	
	
	public void printDensityToFile(String filename){
		try {
			FileWriter writer = new FileWriter(filename);
			writer.write("#watson\n");
			for (int i=0; i <watsonProbs.length ; i++)
				writer.write(watsonProbs[i]+",");
			writer.write("\n");
			writer.write("#crick\n");
			for (int i=0; i< crickProbs.length; i++)
				writer.write(crickProbs[i]+",");
			writer.write("\n");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	//Initialize the data structures
	protected void init(int left, int right){
		this.left = left;
		this.right = right;
		winSize = right-left+1;
		watsonData = new double [winSize];
		crickData = new double [winSize];
		watsonProbs = new double [winSize];
		crickProbs = new double [winSize];
		watsonLogProbs = new double [winSize];
		crickLogProbs = new double [winSize];
		for(int w=0; w<winSize; w++){
			watsonData[w]=0; crickData[w]=0;
			watsonProbs[w]=0; crickProbs[w]=0;
			watsonLogProbs[w]=Double.NEGATIVE_INFINITY; crickLogProbs[w]=Double.NEGATIVE_INFINITY;
		}
	}
	
	//Return a pair of distances corresponding to the central probability interval provided
	protected Pair<Integer,Integer> probIntervalDistances(double prob){
		double ends=(1-prob)/2;
		double probSum=0;
		boolean firstFound=false, secondFound=false;
		int first=left, second=right;
		for(int i=left; i<=right; i++){
			probSum+=probability(i, true)+probability(i, false);
			if(!firstFound && probSum>ends){
				firstFound=true;
				first=i;
			}else if(!secondFound && probSum>(1-ends)){
				secondFound=true;
				second=i;
			}
		}
		Pair<Integer,Integer> intervalDists = new Pair<Integer, Integer>(first, second);
		return intervalDists;
	}
	
	//Update the influence range
	protected void updateInfluenceRange(){
		Pair<Integer,Integer> intervals = probIntervalDistances(0.95);
		int longest = Math.max(Math.abs(intervals.car()), Math.abs(intervals.cdr()));
		influenceRange = longest;
	}
	
	//Look up the data corresponding to a distance
	protected double dataVal(int distance, boolean watsonStrand){
		if(distance<left || distance>right){
			return(0.0);
		}else{
			return(watsonStrand ? 
					watsonData[distance-left] : crickData[distance-left]);
		}
	}
	
	/**
	 * Clone this TagProbabilityDensity
	 */
	public TagProbabilityDensity clone(){
		TagProbabilityDensity newTDS = null;
		try {
			newTDS = new TagProbabilityDensity(this.left, this.right);
			newTDS.isGaussian=this.isGaussian;
			newTDS.gaussOffsetW = this.gaussOffsetW;
			newTDS.gaussOffsetC = this.gaussOffsetC;
			newTDS.gaussSigmaW = this.gaussSigmaW;
			newTDS.gaussSigmaC = this.gaussSigmaC;
			newTDS.loadData(watsonData, crickData);
			newTDS.makeProbabilities();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return newTDS;
	}
	
	/**
	 * Reverse complement
	 * @return
	 */
	public TagProbabilityDensity reverseComplement(){
		TagProbabilityDensity newTDS = null;
		try {
			newTDS = new TagProbabilityDensity(this.left, this.right);
			newTDS.isGaussian=this.isGaussian;
			newTDS.gaussOffsetC = this.gaussOffsetW;
			newTDS.gaussOffsetW = this.gaussOffsetC;
			newTDS.gaussSigmaW = this.gaussSigmaC;
			newTDS.gaussSigmaC = this.gaussSigmaW;
			double []revW = new double[watsonData.length];
			double []revC = new double[crickData.length];
			for(int i=0; i<watsonData.length; i++){
				revW[i] = watsonData[watsonData.length-1-i];
				revC[i] = crickData[crickData.length-1-i];
			}
			newTDS.loadData(revC, revW);
			newTDS.makeProbabilities();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return newTDS;
	}
	
	/**
	 * Save this TagProbabilityDensity
	 */
	public String toString(){
		String out="";
		for(int w=0; w<winSize; w++){
			out = out + w +"\t"+watsonProbs[w]+"\t"+crickProbs[w]+"\n";
		}
		return out;
	}
	
	/**
	 * Save this TagProbabilityDensity
	 */
	public String saveString(){
		String out="#TagProbabilityDensity,"+index+","+left+","+right+","+winSize+",\n";
		if(this.isGaussian){
			out = out+"GaussianW,"+gaussOffsetW+","+gaussSigmaW+",\n";
			out = out+"GaussianC,"+gaussOffsetC+","+gaussSigmaC+",\n";
		}else{
			out = out+"EmpiricalW,";
			for(int w=0; w<winSize; w++){
				out = out + watsonProbs[w]+",";
			}out = out +"\n";
			out = out+"EmpiricalC,";
			for(int w=0; w<winSize; w++){
				out = out + crickProbs[w]+",";
			}out = out +"\n";
		}
		return out;
	}
	
	/**
	 * Load a TagProbabilityDensity using the same format as output in TagProbabilityDensity.saveString()
	 * @param lines: List of 3 lines describing the TPD
	 * @return
	 */
	public static TagProbabilityDensity load(List<String> lines){
		TagProbabilityDensity tpd=null;
		if(lines.size()!=3){
			System.err.println("TagProbabilityDensity.load(): Unexpected format");
			System.exit(1);
		}else{
			Integer index, l, r, w;
			String[] bits = lines.get(0).split(",");
			String[] bits2 = lines.get(1).split(",");
			String[] bits3 = lines.get(2).split(",");
			if(!bits[0].equals("#TagProbabilityDensity") || bits.length!=5){
				System.err.println("TagProbabilityDensity.load(): Unexpected format");
				System.exit(1);
			}
			index = new Integer(bits[1]);
			l = new Integer(bits[2]);
			r = new Integer(bits[3]);
			w = new Integer(bits[4]);
			tpd = new TagProbabilityDensity(l, r);
			tpd.setIndex(index);
			if(bits2[0].equals("GaussianW") && bits3[0].equals("GaussianC")){
				Double gow = new Double(bits2[1]);
				Double gsw = new Double(bits2[2]);
				Double goc = new Double(bits3[1]);
				Double gsc = new Double(bits3[2]);
				tpd.loadGaussianDistrib(gow, gsw, goc, gsc);
			}else if(bits2[0].equals("EmpiricalW") && bits3[0].equals("EmpiricalC") && bits2.length>w && bits3.length>w){
				double[] watson = new double[w];
				double[] crick = new double[w];
				for(int i=0; i<w; i++){
					watson[i] = new Double(bits2[i+1]);
					crick[i] = new Double(bits3[i+1]);
				}
				try {
					tpd.loadData(watson, crick);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}else{
				System.err.println("TagProbabilityDensity.load(): Unexpected format");
				System.exit(1);
			}
		}
		return tpd;
	}
	
	//Main
	public static void main(String[] args){
		int win = 200;
		TagProbabilityDensity td = new TagProbabilityDensity(win);
		td.loadGaussianDistrib(-6, 1, 6, 1);
		
		double[] watson = td.getWatsonProbabilities();
		double[] crick = td.getCrickProbabilities();
		for(int w=0; w<win; w++){
			int pos = (w+td.getLeft());
			System.out.println(pos+"\t"+watson[w]+"\t"+crick[w]);
		}
	}
	
	
	@SuppressWarnings("unchecked")
	public static final List<Pair<Integer, Double>> defaultChipSeqEmpiricalDistribution =
	        Arrays.asList(new Pair<Integer,Double>(-250, 2.1453194893774935E-4),
	        		new Pair<Integer,Double>(-249, 2.152548879235857E-4),
	        		new Pair<Integer,Double>(-248, 2.1598911849088055E-4),
	        		new Pair<Integer,Double>(-247, 2.167459322210923E-4),
	        		new Pair<Integer,Double>(-246, 2.1753662069567956E-4),
	        		new Pair<Integer,Double>(-245, 2.1837247549610073E-4),
	        		new Pair<Integer,Double>(-244, 2.1926478820381436E-4),
	        		new Pair<Integer,Double>(-243, 2.2022485040027885E-4),
	        		new Pair<Integer,Double>(-242, 2.2126395366695278E-4),
	        		new Pair<Integer,Double>(-241, 2.223933895852946E-4),
	        		new Pair<Integer,Double>(-240, 2.2362444973676274E-4),
	        		new Pair<Integer,Double>(-239, 2.2496842570281572E-4),
	        		new Pair<Integer,Double>(-238, 2.2643660906491208E-4),
	        		new Pair<Integer,Double>(-237, 2.280402914045103E-4),
	        		new Pair<Integer,Double>(-236, 2.2979076430306873E-4),
	        		new Pair<Integer,Double>(-235, 2.31699319342046E-4),
	        		new Pair<Integer,Double>(-234, 2.337772481029006E-4),
	        		new Pair<Integer,Double>(-233, 2.360358421670909E-4),
	        		new Pair<Integer,Double>(-232, 2.3848639311607544E-4),
	        		new Pair<Integer,Double>(-231, 2.4114019253131268E-4),
	        		new Pair<Integer,Double>(-230, 2.440085319942612E-4),
	        		new Pair<Integer,Double>(-229, 2.4710270308637945E-4),
	        		new Pair<Integer,Double>(-228, 2.504339973891258E-4),
	        		new Pair<Integer,Double>(-227, 2.5401370648395883E-4),
	        		new Pair<Integer,Double>(-226, 2.5785312195233717E-4),
	        		new Pair<Integer,Double>(-225, 2.6196353537571906E-4),
	        		new Pair<Integer,Double>(-224, 2.6635623833556296E-4),
	        		new Pair<Integer,Double>(-223, 2.710425224133276E-4),
	        		new Pair<Integer,Double>(-222, 2.760336791904713E-4),
	        		new Pair<Integer,Double>(-221, 2.8134100024845263E-4),
	        		new Pair<Integer,Double>(-220, 2.8697577716872995E-4),
	        		new Pair<Integer,Double>(-219, 2.929493015327619E-4),
	        		new Pair<Integer,Double>(-218, 2.9927286492200685E-4),
	        		new Pair<Integer,Double>(-217, 3.059577589179232E-4),
	        		new Pair<Integer,Double>(-216, 3.130152751019698E-4),
	        		new Pair<Integer,Double>(-215, 3.2045670505560477E-4),
	        		new Pair<Integer,Double>(-214, 3.282933403602868E-4),
	        		new Pair<Integer,Double>(-213, 3.3653647259747414E-4),
	        		new Pair<Integer,Double>(-212, 3.4519739334862553E-4),
	        		new Pair<Integer,Double>(-211, 3.542873941951994E-4),
	        		new Pair<Integer,Double>(-210, 3.6381776671865406E-4),
	        		new Pair<Integer,Double>(-209, 3.737998025004482E-4),
	        		new Pair<Integer,Double>(-208, 3.842447931220402E-4),
	        		new Pair<Integer,Double>(-207, 3.9516403016488875E-4),
	        		new Pair<Integer,Double>(-206, 4.0656880521045186E-4),
	        		new Pair<Integer,Double>(-205, 4.184704098401885E-4),
	        		new Pair<Integer,Double>(-204, 4.3087821722982814E-4),
	        		new Pair<Integer,Double>(-203, 4.4379392693218514E-4),
	        		new Pair<Integer,Double>(-202, 4.5721732009434523E-4),
	        		new Pair<Integer,Double>(-201, 4.7114817786339373E-4),
	        		new Pair<Integer,Double>(-200, 4.8558628138641643E-4),
	        		new Pair<Integer,Double>(-199, 5.005314118104988E-4),
	        		new Pair<Integer,Double>(-198, 5.159833502827266E-4),
	        		new Pair<Integer,Double>(-197, 5.319418779501851E-4),
	        		new Pair<Integer,Double>(-196, 5.484067759599599E-4),
	        		new Pair<Integer,Double>(-195, 5.653778254591369E-4),
	        		new Pair<Integer,Double>(-194, 5.828548075948015E-4),
	        		new Pair<Integer,Double>(-193, 6.008375035140394E-4),
	        		new Pair<Integer,Double>(-192, 6.193256943639358E-4),
	        		new Pair<Integer,Double>(-191, 6.383191612915764E-4),
	        		new Pair<Integer,Double>(-190, 6.578176854440472E-4),
	        		new Pair<Integer,Double>(-189, 6.778210479684333E-4),
	        		new Pair<Integer,Double>(-188, 6.983290300118204E-4),
	        		new Pair<Integer,Double>(-187, 7.193414127212942E-4),
	        		new Pair<Integer,Double>(-186, 7.4085797724394E-4),
	        		new Pair<Integer,Double>(-185, 7.628785047268438E-4),
	        		new Pair<Integer,Double>(-184, 7.854027763170908E-4),
	        		new Pair<Integer,Double>(-183, 8.084305731617667E-4),
	        		new Pair<Integer,Double>(-182, 8.319616764079573E-4),
	        		new Pair<Integer,Double>(-181, 8.559958672027479E-4),
	        		new Pair<Integer,Double>(-180, 8.80532926693224E-4),
	        		new Pair<Integer,Double>(-179, 9.055726360264715E-4),
	        		new Pair<Integer,Double>(-178, 9.311147763495754E-4),
	        		new Pair<Integer,Double>(-177, 9.57159128809622E-4),
	        		new Pair<Integer,Double>(-176, 9.837054745536963E-4),
	        		new Pair<Integer,Double>(-175, 0.0010107535947288845),
	        		new Pair<Integer,Double>(-174, 0.0010383052420264377),
	        		new Pair<Integer,Double>(-173, 0.0010663700553142743),
	        		new Pair<Integer,Double>(-172, 0.0010949596450044764),
	        		new Pair<Integer,Double>(-171, 0.001124085621509129),
	        		new Pair<Integer,Double>(-170, 0.0011537595952403152),
	        		new Pair<Integer,Double>(-169, 0.0011839931766101186),
	        		new Pair<Integer,Double>(-168, 0.0012147979760306227),
	        		new Pair<Integer,Double>(-167, 0.0012461856039139116),
	        		new Pair<Integer,Double>(-166, 0.0012781676706720686),
	        		new Pair<Integer,Double>(-165, 0.0013107557867171776),
	        		new Pair<Integer,Double>(-164, 0.0013439615624613215),
	        		new Pair<Integer,Double>(-163, 0.0013777966083165852),
	        		new Pair<Integer,Double>(-162, 0.001412272534695051),
	        		new Pair<Integer,Double>(-161, 0.0014474009520088038),
	        		new Pair<Integer,Double>(-160, 0.001483193470669926),
	        		new Pair<Integer,Double>(-159, 0.001519661701090502),
	        		new Pair<Integer,Double>(-158, 0.0015568172536826151),
	        		new Pair<Integer,Double>(-157, 0.0015946717388583497),
	        		new Pair<Integer,Double>(-156, 0.0016332367670297881),
	        		new Pair<Integer,Double>(-155, 0.0016725239486090147),
	        		new Pair<Integer,Double>(-154, 0.0017125448940081137),
	        		new Pair<Integer,Double>(-153, 0.0017533112136391678),
	        		new Pair<Integer,Double>(-152, 0.0017948345179142613),
	        		new Pair<Integer,Double>(-151, 0.0018371264172454772),
	        		new Pair<Integer,Double>(-150, 0.001880198522044899),
	        		new Pair<Integer,Double>(-149, 0.0019240624427246116),
	        		new Pair<Integer,Double>(-148, 0.001968729789696697),
	        		new Pair<Integer,Double>(-147, 0.0020142121733732407),
	        		new Pair<Integer,Double>(-146, 0.002060521204166324),
	        		new Pair<Integer,Double>(-145, 0.0021076684924880326),
	        		new Pair<Integer,Double>(-144, 0.00215566144256165),
	        		new Pair<Integer,Double>(-143, 0.002204490633855266),
	        		new Pair<Integer,Double>(-142, 0.00225414243964817),
	        		new Pair<Integer,Double>(-141, 0.0023046032332196516),
	        		new Pair<Integer,Double>(-140, 0.0023558593878490013),
	        		new Pair<Integer,Double>(-139, 0.002407897276815508),
	        		new Pair<Integer,Double>(-138, 0.0024607032733984626),
	        		new Pair<Integer,Double>(-137, 0.002514263750877155),
	        		new Pair<Integer,Double>(-136, 0.002568565082530874),
	        		new Pair<Integer,Double>(-135, 0.002623593641638911),
	        		new Pair<Integer,Double>(-134, 0.002679335801480555),
	        		new Pair<Integer,Double>(-133, 0.0027357779353350967),
	        		new Pair<Integer,Double>(-132, 0.002792906416481825),
	        		new Pair<Integer,Double>(-131, 0.00285070761820003),
	        		new Pair<Integer,Double>(-130, 0.0029091679137690026),
	        		new Pair<Integer,Double>(-129, 0.002968273676468031),
	        		new Pair<Integer,Double>(-128, 0.003028011279576406),
	        		new Pair<Integer,Double>(-127, 0.0030883670963734186),
	        		new Pair<Integer,Double>(-126, 0.003149327500138357),
	        		new Pair<Integer,Double>(-125, 0.003210878864150512),
	        		new Pair<Integer,Double>(-124, 0.0032730075616891736),
	        		new Pair<Integer,Double>(-123, 0.0033356999660336304),
	        		new Pair<Integer,Double>(-122, 0.0033989424504631744),
	        		new Pair<Integer,Double>(-121, 0.003462721388257094),
	        		new Pair<Integer,Double>(-120, 0.0035270231526946796),
	        		new Pair<Integer,Double>(-119, 0.0035918341170552213),
	        		new Pair<Integer,Double>(-118, 0.0036571406546180084),
	        		new Pair<Integer,Double>(-117, 0.0037229291386623306),
	        		new Pair<Integer,Double>(-116, 0.003789185942467479),
	        		new Pair<Integer,Double>(-115, 0.003855897439312743),
	        		new Pair<Integer,Double>(-114, 0.003923043602622363),
	        		new Pair<Integer,Double>(-113, 0.003990578806400385),
	        		new Pair<Integer,Double>(-112, 0.004058451024795804),
	        		new Pair<Integer,Double>(-111, 0.0041266082319576175),
	        		new Pair<Integer,Double>(-110, 0.00419499840203482),
	        		new Pair<Integer,Double>(-109, 0.0042635695091764085),
	        		new Pair<Integer,Double>(-108, 0.004332269527531378),
	        		new Pair<Integer,Double>(-107, 0.004401046431248725),
	        		new Pair<Integer,Double>(-106, 0.004469848194477447),
	        		new Pair<Integer,Double>(-105, 0.004538622791366541),
	        		new Pair<Integer,Double>(-104, 0.004607318196064999),
	        		new Pair<Integer,Double>(-103, 0.0046758823827218205),
	        		new Pair<Integer,Double>(-102, 0.004744263325486001),
	        		new Pair<Integer,Double>(-101, 0.004812408998506535),
	        		new Pair<Integer,Double>(-100, 0.004880267375932421),
	        		new Pair<Integer,Double>(-99, 0.004947786431912652),
	        		new Pair<Integer,Double>(-98, 0.005014914140596227),
	        		new Pair<Integer,Double>(-97, 0.00508159847613214),
	        		new Pair<Integer,Double>(-96, 0.005147787412669388),
	        		new Pair<Integer,Double>(-95, 0.005213428924356968),
	        		new Pair<Integer,Double>(-94, 0.0052784709853438755),
	        		new Pair<Integer,Double>(-93, 0.005342861569779106),
	        		new Pair<Integer,Double>(-92, 0.005406548651811655),
	        		new Pair<Integer,Double>(-91, 0.00546948020559052),
	        		new Pair<Integer,Double>(-90, 0.0055316042052646975),
	        		new Pair<Integer,Double>(-89, 0.0055928686249831815),
	        		new Pair<Integer,Double>(-88, 0.00565322143889497),
	        		new Pair<Integer,Double>(-87, 0.005712610621149058),
	        		new Pair<Integer,Double>(-86, 0.005770984145894442),
	        		new Pair<Integer,Double>(-85, 0.005828289987280117),
	        		new Pair<Integer,Double>(-84, 0.005884484079524634),
	        		new Pair<Integer,Double>(-83, 0.0059395541971247485),
	        		new Pair<Integer,Double>(-82, 0.005993496074646772),
	        		new Pair<Integer,Double>(-81, 0.006046305446657014),
	        		new Pair<Integer,Double>(-80, 0.006097978047721785),
	        		new Pair<Integer,Double>(-79, 0.006148509612407397),
	        		new Pair<Integer,Double>(-78, 0.0061978958752801575),
	        		new Pair<Integer,Double>(-77, 0.0062461325709063775),
	        		new Pair<Integer,Double>(-76, 0.006293215433852369),
	        		new Pair<Integer,Double>(-75, 0.0063391401986844405),
	        		new Pair<Integer,Double>(-74, 0.006383902599968901),
	        		new Pair<Integer,Double>(-73, 0.006427498372272066),
	        		new Pair<Integer,Double>(-72, 0.00646992325016024),
	        		new Pair<Integer,Double>(-71, 0.0065111729681997365),
	        		new Pair<Integer,Double>(-70, 0.006551243260956864),
	        		new Pair<Integer,Double>(-69, 0.006590129862997933),
	        		new Pair<Integer,Double>(-68, 0.006627828508889254),
	        		new Pair<Integer,Double>(-67, 0.006664334933197137),
	        		new Pair<Integer,Double>(-66, 0.006699644870487895),
	        		new Pair<Integer,Double>(-65, 0.006733754055327833),
	        		new Pair<Integer,Double>(-64, 0.006766658222283267),
	        		new Pair<Integer,Double>(-63, 0.006798353105920504),
	        		new Pair<Integer,Double>(-62, 0.006828834440805855),
	        		new Pair<Integer,Double>(-61, 0.00685809796150563),
	        		new Pair<Integer,Double>(-60, 0.006886139402586138),
	        		new Pair<Integer,Double>(-59, 0.006912954498613692),
	        		new Pair<Integer,Double>(-58, 0.0069385389841546),
	        		new Pair<Integer,Double>(-57, 0.0069628885937751735),
	        		new Pair<Integer,Double>(-56, 0.006985999062041721),
	        		new Pair<Integer,Double>(-55, 0.0070078661235205555),
	        		new Pair<Integer,Double>(-54, 0.007028467255690983),
	        		new Pair<Integer,Double>(-53, 0.007047706907684303),
	        		new Pair<Integer,Double>(-52, 0.007065471271544807),
	        		new Pair<Integer,Double>(-51, 0.007081646539316792),
	        		new Pair<Integer,Double>(-50, 0.007096118903044555),
	        		new Pair<Integer,Double>(-49, 0.00710877455477239),
	        		new Pair<Integer,Double>(-48, 0.007119499686544591),
	        		new Pair<Integer,Double>(-47, 0.007128180490405454),
	        		new Pair<Integer,Double>(-46, 0.007134703158399274),
	        		new Pair<Integer,Double>(-45, 0.007138953882570347),
	        		new Pair<Integer,Double>(-44, 0.007140818854962967),
	        		new Pair<Integer,Double>(-43, 0.00714018426762143),
	        		new Pair<Integer,Double>(-42, 0.007136936312590032),
	        		new Pair<Integer,Double>(-41, 0.007130961181913066),
	        		new Pair<Integer,Double>(-40, 0.007122145067634828),
	        		new Pair<Integer,Double>(-39, 0.007110374161799614),
	        		new Pair<Integer,Double>(-38, 0.007095534656451718),
	        		new Pair<Integer,Double>(-37, 0.007077512743635436),
	        		new Pair<Integer,Double>(-36, 0.007056194615395063),
	        		new Pair<Integer,Double>(-35, 0.0070314664637748945),
	        		new Pair<Integer,Double>(-34, 0.0070032144808192245),
	        		new Pair<Integer,Double>(-33, 0.0069713248585723505),
	        		new Pair<Integer,Double>(-32, 0.006935683789078565),
	        		new Pair<Integer,Double>(-31, 0.006896177464382165),
	        		new Pair<Integer,Double>(-30, 0.006852692076527444),
	        		new Pair<Integer,Double>(-29, 0.0068051138175586985),
	        		new Pair<Integer,Double>(-28, 0.006753328879520223),
	        		new Pair<Integer,Double>(-27, 0.006697223454456314),
	        		new Pair<Integer,Double>(-26, 0.006636683734411266),
	        		new Pair<Integer,Double>(-25, 0.006571595911429373),
	        		new Pair<Integer,Double>(-24, 0.006501901865600835),
	        		new Pair<Integer,Double>(-23, 0.006427766229199461),
	        		new Pair<Integer,Double>(-22, 0.006349409322544967),
	        		new Pair<Integer,Double>(-21, 0.006267051465957067),
	        		new Pair<Integer,Double>(-20, 0.006180912979755476),
	        		new Pair<Integer,Double>(-19, 0.006091214184259908),
	        		new Pair<Integer,Double>(-18, 0.005998175399790077),
	        		new Pair<Integer,Double>(-17, 0.005902016946665696),
	        		new Pair<Integer,Double>(-16, 0.005802959145206481),
	        		new Pair<Integer,Double>(-15, 0.005701222315732148),
	        		new Pair<Integer,Double>(-14, 0.005597026778562408),
	        		new Pair<Integer,Double>(-13, 0.005490592854016977),
	        		new Pair<Integer,Double>(-12, 0.005382140862415569),
	        		new Pair<Integer,Double>(-11, 0.005271891124077898),
	        		new Pair<Integer,Double>(-10, 0.00516006395932368),
	        		new Pair<Integer,Double>(-9, 0.005046879688472627),
	        		new Pair<Integer,Double>(-8, 0.004932558631844455),
	        		new Pair<Integer,Double>(-7, 0.004817321109758877),
	        		new Pair<Integer,Double>(-6, 0.004701387442535608),
	        		new Pair<Integer,Double>(-5, 0.0045849779504943625),
	        		new Pair<Integer,Double>(-4, 0.004468312953954856),
	        		new Pair<Integer,Double>(-3, 0.0043516127732368005),
	        		new Pair<Integer,Double>(-2, 0.004235097728659913),
	        		new Pair<Integer,Double>(-1, 0.004118988140543906),
	        		new Pair<Integer,Double>(0, 0.004003504329208493),
	        		new Pair<Integer,Double>(1, 0.0038888666149733894),
	        		new Pair<Integer,Double>(2, 0.00377529531815831),
	        		new Pair<Integer,Double>(3, 0.0036630107590829688),
	        		new Pair<Integer,Double>(4, 0.00355223325806708),
	        		new Pair<Integer,Double>(5, 0.003443183135430358),
	        		new Pair<Integer,Double>(6, 0.003336044668723669),
	        		new Pair<Integer,Double>(7, 0.0032308579644224874),
	        		new Pair<Integer,Double>(8, 0.003127627086233437),
	        		new Pair<Integer,Double>(9, 0.0030263560978631447),
	        		new Pair<Integer,Double>(10, 0.0029270490630182357),
	        		new Pair<Integer,Double>(11, 0.002829710045405335),
	        		new Pair<Integer,Double>(12, 0.002734343108731069),
	        		new Pair<Integer,Double>(13, 0.0026409523167020615),
	        		new Pair<Integer,Double>(14, 0.0025495417330249396),
	        		new Pair<Integer,Double>(15, 0.002460115421406328),
	        		new Pair<Integer,Double>(16, 0.0023726774455528534),
	        		new Pair<Integer,Double>(17, 0.0022872318691711394),
	        		new Pair<Integer,Double>(18, 0.002203782755967813),
	        		new Pair<Integer,Double>(19, 0.0021223341696494998),
	        		new Pair<Integer,Double>(20, 0.002042890173922824),
	        		new Pair<Integer,Double>(21, 0.0019654548324944117),
	        		new Pair<Integer,Double>(22, 0.0018900322090708885),
	        		new Pair<Integer,Double>(23, 0.00181662636735888),
	        		new Pair<Integer,Double>(24, 0.0017452413710650112),
	        		new Pair<Integer,Double>(25, 0.001675881283895908),
	        		new Pair<Integer,Double>(26, 0.001608550169558196),
	        		new Pair<Integer,Double>(27, 0.0015432520917585007),
	        		new Pair<Integer,Double>(28, 0.0014799911142034472),
	        		new Pair<Integer,Double>(29, 0.0014187713005996614),
	        		new Pair<Integer,Double>(30, 0.001359596714653768),
	        		new Pair<Integer,Double>(31, 0.0013024714200723937),
	        		new Pair<Integer,Double>(32, 0.0012473994805621634),
	        		new Pair<Integer,Double>(33, 0.0011943849598297022),
	        		new Pair<Integer,Double>(34, 0.001143431921581636),
	        		new Pair<Integer,Double>(35, 0.0010945444295245902),
	        		new Pair<Integer,Double>(36, 0.0010477178032926953),
	        		new Pair<Integer,Double>(37, 0.0010029123862301005),
	        		new Pair<Integer,Double>(38, 9.600797776084589E-4),
	        		new Pair<Integer,Double>(39, 9.191715766994251E-4),
	        		new Pair<Integer,Double>(40, 8.801393827746526E-4),
	        		new Pair<Integer,Double>(41, 8.429347951057955E-4),
	        		new Pair<Integer,Double>(42, 8.075094129645075E-4),
	        		new Pair<Integer,Double>(43, 7.738148356224423E-4),
	        		new Pair<Integer,Double>(44, 7.418026623512542E-4),
	        		new Pair<Integer,Double>(45, 7.114244924225969E-4),
	        		new Pair<Integer,Double>(46, 6.826319251081241E-4),
	        		new Pair<Integer,Double>(47, 6.5537655967949E-4),
	        		new Pair<Integer,Double>(48, 6.296099954083481E-4),
	        		new Pair<Integer,Double>(49, 6.052838315663525E-4),
	        		new Pair<Integer,Double>(50, 5.823496674251568E-4),
	        		new Pair<Integer,Double>(51, 5.607591022564153E-4),
	        		new Pair<Integer,Double>(52, 5.404637353317815E-4),
	        		new Pair<Integer,Double>(53, 5.214151659229095E-4),
	        		new Pair<Integer,Double>(54, 5.035649933014529E-4),
	        		new Pair<Integer,Double>(55, 4.8686481673906575E-4),
	        		new Pair<Integer,Double>(56, 4.7126623550740196E-4),
	        		new Pair<Integer,Double>(57, 4.567208488781154E-4),
	        		new Pair<Integer,Double>(58, 4.431802561228598E-4),
	        		new Pair<Integer,Double>(59, 4.30596056513289E-4),
	        		new Pair<Integer,Double>(60, 4.189198493210571E-4),
	        		new Pair<Integer,Double>(61, 4.081032338178178E-4),
	        		new Pair<Integer,Double>(62, 3.9809780927522494E-4),
	        		new Pair<Integer,Double>(63, 3.888551749649324E-4),
	        		new Pair<Integer,Double>(64, 3.803269301585942E-4),
	        		new Pair<Integer,Double>(65, 3.7246467412786397E-4),
	        		new Pair<Integer,Double>(66, 3.6522384180415175E-4),
	        		new Pair<Integer,Double>(67, 3.585752107578902E-4),
	        		new Pair<Integer,Double>(68, 3.524933942192684E-4),
	        		new Pair<Integer,Double>(69, 3.4695300541847507E-4),
	        		new Pair<Integer,Double>(70, 3.4192865758569924E-4),
	        		new Pair<Integer,Double>(71, 3.3739496395112967E-4),
	        		new Pair<Integer,Double>(72, 3.3332653774495517E-4),
	        		new Pair<Integer,Double>(73, 3.296979921973647E-4),
	        		new Pair<Integer,Double>(74, 3.2648394053854714E-4),
	        		new Pair<Integer,Double>(75, 3.2365899599869134E-4),
	        		new Pair<Integer,Double>(76, 3.2119777180798613E-4),
	        		new Pair<Integer,Double>(77, 3.190748811966203E-4),
	        		new Pair<Integer,Double>(78, 3.1726493739478293E-4),
	        		new Pair<Integer,Double>(79, 3.1574255363266265E-4),
	        		new Pair<Integer,Double>(80, 3.1448234314044844E-4),
	        		new Pair<Integer,Double>(81, 3.1345891914832914E-4),
	        		new Pair<Integer,Double>(82, 3.126468948864936E-4),
	        		new Pair<Integer,Double>(83, 3.120208835851308E-4),
	        		new Pair<Integer,Double>(84, 3.115554984744293E-4),
	        		new Pair<Integer,Double>(85, 3.112253527845783E-4),
	        		new Pair<Integer,Double>(86, 3.1100505974576654E-4),
	        		new Pair<Integer,Double>(87, 3.1086923258818284E-4),
	        		new Pair<Integer,Double>(88, 3.1079248454201605E-4),
	        		new Pair<Integer,Double>(89, 3.1074942883745516E-4),
	        		new Pair<Integer,Double>(90, 3.1071467870468893E-4),
	        		new Pair<Integer,Double>(91, 3.1066284737390625E-4),
	        		new Pair<Integer,Double>(92, 3.105685480752959E-4),
	        		new Pair<Integer,Double>(93, 3.104063940390469E-4),
	        		new Pair<Integer,Double>(94, 3.1015099849534807E-4),
	        		new Pair<Integer,Double>(95, 3.097769746743882E-4),
	        		new Pair<Integer,Double>(96, 3.092643687909541E-4),
	        		new Pair<Integer,Double>(97, 3.086149589982245E-4),
	        		new Pair<Integer,Double>(98, 3.078359564339763E-4),
	        		new Pair<Integer,Double>(99, 3.069345722359859E-4),
	        		new Pair<Integer,Double>(100, 3.059180175420301E-4),
	        		new Pair<Integer,Double>(101, 3.047935034898855E-4),
	        		new Pair<Integer,Double>(102, 3.0356824121732887E-4),
	        		new Pair<Integer,Double>(103, 3.022494418621368E-4),
	        		new Pair<Integer,Double>(104, 3.008443165620861E-4),
	        		new Pair<Integer,Double>(105, 2.993600764549533E-4),
	        		new Pair<Integer,Double>(106, 2.9780393267851517E-4),
	        		new Pair<Integer,Double>(107, 2.9618309637054847E-4),
	        		new Pair<Integer,Double>(108, 2.945047786688296E-4),
	        		new Pair<Integer,Double>(109, 2.9277619071113554E-4),
	        		new Pair<Integer,Double>(110, 2.910045436352428E-4),
	        		new Pair<Integer,Double>(111, 2.8919704857892814E-4),
	        		new Pair<Integer,Double>(112, 2.8736091667996817E-4),
	        		new Pair<Integer,Double>(113, 2.8550335907613953E-4),
	        		new Pair<Integer,Double>(114, 2.83631586905219E-4),
	        		new Pair<Integer,Double>(115, 2.817528113049833E-4),
	        		new Pair<Integer,Double>(116, 2.7987424341320897E-4),
	        		new Pair<Integer,Double>(117, 2.780030943676727E-4),
	        		new Pair<Integer,Double>(118, 2.7614657530615136E-4),
	        		new Pair<Integer,Double>(119, 2.7431189736642136E-4),
	        		new Pair<Integer,Double>(120, 2.725062716862596E-4),
	        		new Pair<Integer,Double>(121, 2.707369094034427E-4),
	        		new Pair<Integer,Double>(122, 2.6901102165574717E-4),
	        		new Pair<Integer,Double>(123, 2.6733581958095E-4),
	        		new Pair<Integer,Double>(124, 2.6571851431682754E-4),
	        		new Pair<Integer,Double>(125, 2.641663170011567E-4),
	        		new Pair<Integer,Double>(126, 2.6268485293332744E-4),
	        		new Pair<Integer,Double>(127, 2.612734040591833E-4),
	        		new Pair<Integer,Double>(128, 2.5992966648618145E-4),
	        		new Pair<Integer,Double>(129, 2.586513363217787E-4),
	        		new Pair<Integer,Double>(130, 2.574361096734321E-4),
	        		new Pair<Integer,Double>(131, 2.5628168264859853E-4),
	        		new Pair<Integer,Double>(132, 2.551857513547352E-4),
	        		new Pair<Integer,Double>(133, 2.5414601189929884E-4),
	        		new Pair<Integer,Double>(134, 2.5316016038974656E-4),
	        		new Pair<Integer,Double>(135, 2.5222589293353524E-4),
	        		new Pair<Integer,Double>(136, 2.5134090563812197E-4),
	        		new Pair<Integer,Double>(137, 2.5050289461096365E-4),
	        		new Pair<Integer,Double>(138, 2.4970955595951733E-4),
	        		new Pair<Integer,Double>(139, 2.4895858579123987E-4),
	        		new Pair<Integer,Double>(140, 2.4824768021358836E-4),
	        		new Pair<Integer,Double>(141, 2.475745353340198E-4),
	        		new Pair<Integer,Double>(142, 2.4693684725999107E-4),
	        		new Pair<Integer,Double>(143, 2.463323120989591E-4),
	        		new Pair<Integer,Double>(144, 2.45758625958381E-4),
	        		new Pair<Integer,Double>(145, 2.452134849457137E-4),
	        		new Pair<Integer,Double>(146, 2.446945851684142E-4),
	        		new Pair<Integer,Double>(147, 2.441996227339394E-4),
	        		new Pair<Integer,Double>(148, 2.4372629374974647E-4),
	        		new Pair<Integer,Double>(149, 2.4327229432329214E-4),
	        		new Pair<Integer,Double>(150, 2.4283532056203352E-4),
	        		new Pair<Integer,Double>(151, 2.4241306857342756E-4),
	        		new Pair<Integer,Double>(152, 2.4200323446493127E-4),
	        		new Pair<Integer,Double>(153, 2.4160351434400154E-4),
	        		new Pair<Integer,Double>(154, 2.4121160431809548E-4),
	        		new Pair<Integer,Double>(155, 2.4082520049466992E-4),
	        		new Pair<Integer,Double>(156, 2.4044231832542033E-4),
	        		new Pair<Integer,Double>(157, 2.4006225063899538E-4),
	        		new Pair<Integer,Double>(158, 2.3968460960828208E-4),
	        		new Pair<Integer,Double>(159, 2.3930900740616758E-4),
	        		new Pair<Integer,Double>(160, 2.389350562055389E-4),
	        		new Pair<Integer,Double>(161, 2.3856236817928317E-4),
	        		new Pair<Integer,Double>(162, 2.3819055550028742E-4),
	        		new Pair<Integer,Double>(163, 2.378192303414387E-4),
	        		new Pair<Integer,Double>(164, 2.3744800487562412E-4),
	        		new Pair<Integer,Double>(165, 2.3707649127573078E-4),
	        		new Pair<Integer,Double>(166, 2.367043017146457E-4),
	        		new Pair<Integer,Double>(167, 2.3633104836525597E-4),
	        		new Pair<Integer,Double>(168, 2.3595634340044865E-4),
	        		new Pair<Integer,Double>(169, 2.3557979899311085E-4),
	        		new Pair<Integer,Double>(170, 2.3520102731612954E-4),
	        		new Pair<Integer,Double>(171, 2.3481964054239196E-4),
	        		new Pair<Integer,Double>(172, 2.3443525084478504E-4),
	        		new Pair<Integer,Double>(173, 2.3404747039619593E-4),
	        		new Pair<Integer,Double>(174, 2.3365591136951164E-4),
	        		new Pair<Integer,Double>(175, 2.3326018593761925E-4),
	        		new Pair<Integer,Double>(176, 2.328599062734059E-4),
	        		new Pair<Integer,Double>(177, 2.3245468454975861E-4),
	        		new Pair<Integer,Double>(178, 2.320441329395645E-4),
	        		new Pair<Integer,Double>(179, 2.3162786361571058E-4),
	        		new Pair<Integer,Double>(180, 2.3120548875108394E-4),
	        		new Pair<Integer,Double>(181, 2.3077662051857164E-4),
	        		new Pair<Integer,Double>(182, 2.303408710910608E-4),
	        		new Pair<Integer,Double>(183, 2.298978526414384E-4),
	        		new Pair<Integer,Double>(184, 2.2944717734259162E-4),
	        		new Pair<Integer,Double>(185, 2.289884573674075E-4),
	        		new Pair<Integer,Double>(186, 2.2852139114842897E-4),
	        		new Pair<Integer,Double>(187, 2.2804602215682244E-4),
	        		new Pair<Integer,Double>(188, 2.2756248012341024E-4),
	        		new Pair<Integer,Double>(189, 2.270708947790147E-4),
	        		new Pair<Integer,Double>(190, 2.265713958544582E-4),
	        		new Pair<Integer,Double>(191, 2.260641130805629E-4),
	        		new Pair<Integer,Double>(192, 2.2554917618815122E-4),
	        		new Pair<Integer,Double>(193, 2.250267149080454E-4),
	        		new Pair<Integer,Double>(194, 2.2449685897106786E-4),
	        		new Pair<Integer,Double>(195, 2.239597381080408E-4),
	        		new Pair<Integer,Double>(196, 2.234154820497866E-4),
	        		new Pair<Integer,Double>(197, 2.2286422052712754E-4),
	        		new Pair<Integer,Double>(198, 2.2230608327088595E-4),
	        		new Pair<Integer,Double>(199, 2.217412000118841E-4),
	        		new Pair<Integer,Double>(200, 2.2116970048094436E-4),
	        		new Pair<Integer,Double>(201, 2.2059171440888905E-4),
	        		new Pair<Integer,Double>(202, 2.2000737152654044E-4),
	        		new Pair<Integer,Double>(203, 2.1941680156472082E-4),
	        		new Pair<Integer,Double>(204, 2.1882013425425254E-4),
	        		new Pair<Integer,Double>(205, 2.1821749932595793E-4),
	        		new Pair<Integer,Double>(206, 2.1760902651065926E-4),
	        		new Pair<Integer,Double>(207, 2.1699484553917882E-4),
	        		new Pair<Integer,Double>(208, 2.1637508614233904E-4),
	        		new Pair<Integer,Double>(209, 2.157498780509621E-4),
	        		new Pair<Integer,Double>(210, 2.1511935099587041E-4),
	        		new Pair<Integer,Double>(211, 2.1448363470788622E-4),
	        		new Pair<Integer,Double>(212, 2.1384285891783188E-4),
	        		new Pair<Integer,Double>(213, 2.1319715335652964E-4),
	        		new Pair<Integer,Double>(214, 2.125466477548019E-4),
	        		new Pair<Integer,Double>(215, 2.1189147184347087E-4),
	        		new Pair<Integer,Double>(216, 2.1123175535335897E-4),
	        		new Pair<Integer,Double>(217, 2.1056762801528845E-4),
	        		new Pair<Integer,Double>(218, 2.098992195600816E-4),
	        		new Pair<Integer,Double>(219, 2.092266597185608E-4),
	        		new Pair<Integer,Double>(220, 2.085500782215483E-4),
	        		new Pair<Integer,Double>(221, 2.0786960479986648E-4),
	        		new Pair<Integer,Double>(222, 2.071853691843376E-4),
	        		new Pair<Integer,Double>(223, 2.0649750110578397E-4),
	        		new Pair<Integer,Double>(224, 2.0580613029502791E-4),
	        		new Pair<Integer,Double>(225, 2.0511138648289176E-4),
	        		new Pair<Integer,Double>(226, 2.0441339940019778E-4),
	        		new Pair<Integer,Double>(227, 2.0371229877776833E-4),
	        		new Pair<Integer,Double>(228, 2.0300821434642566E-4),
	        		new Pair<Integer,Double>(229, 2.023012758369922E-4),
	        		new Pair<Integer,Double>(230, 2.0159161298029013E-4),
	        		new Pair<Integer,Double>(231, 2.0087935550714187E-4),
	        		new Pair<Integer,Double>(232, 2.0016463314836963E-4),
	        		new Pair<Integer,Double>(233, 1.9944757563479577E-4),
	        		new Pair<Integer,Double>(234, 1.9872831269724265E-4),
	        		new Pair<Integer,Double>(235, 1.980069740665325E-4),
	        		new Pair<Integer,Double>(236, 1.972836894734877E-4),
	        		new Pair<Integer,Double>(237, 1.9655858864893051E-4),
	        		new Pair<Integer,Double>(238, 1.9583180132368323E-4),
	        		new Pair<Integer,Double>(239, 1.9510345722856829E-4),
	        		new Pair<Integer,Double>(240, 1.9437368609440787E-4),
	        		new Pair<Integer,Double>(241, 1.9364261765202433E-4),
	        		new Pair<Integer,Double>(242, 1.9291038163224E-4),
	        		new Pair<Integer,Double>(243, 1.9217710776587714E-4),
	        		new Pair<Integer,Double>(244, 1.914429257837581E-4),
	        		new Pair<Integer,Double>(245, 1.907079654167052E-4),
	        		new Pair<Integer,Double>(246, 1.8997235639554077E-4),
	        		new Pair<Integer,Double>(247, 1.8923622845108706E-4),
	        		new Pair<Integer,Double>(248, 1.8849971131416638E-4),
	        		new Pair<Integer,Double>(249, 1.8776293471560115E-4),
	        		new Pair<Integer,Double>(250, 1.8702602838621354E-4));
	
}
