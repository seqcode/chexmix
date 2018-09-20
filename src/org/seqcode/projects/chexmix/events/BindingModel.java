package org.seqcode.projects.chexmix.events;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;
import org.seqcode.math.stats.StatUtil;


/**
 * A BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * The probability is directional (i.e. stranded). 
 * Given a signed distance from the event, the BindingModel should return 
 * a relative probability of seeing a read at that distance.
 * 
 * @author shaunmahony
 *
 */
public class BindingModel {
	public final static int SMOOTHING_STEPSIZE = 30;
	public final static int SMOOTHING_AVG_PTS = 30;
	protected final double LOG2 = Math.log(2);
	protected int min, max;		// the start and end position
	protected int summit;		// the position of highest prob point
	protected double[] data;
	protected double[] probs;
	protected double[] logProbs;
	protected double bgProb, logBgProb;
	protected int influenceRange; //95% probability range (over both strands)
	protected String fileName;

	protected List<Pair<Integer, Double>> empiricalDistribution;
	
	protected BindingModel(){};
	
	public BindingModel(File f, int minDist, int maxDist){
		min=0; max=0;
		fileName = f.getName();
		try {
			empiricalDistribution = new LinkedList<Pair<Integer,Double>>(); 
			BufferedReader reader = new BufferedReader(new FileReader(f));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
	              Integer dist = new Integer(words[0]);
	              //make sure the current data point is within the specified range
	              if ((dist.intValue() >= minDist) && (dist.intValue() <= maxDist)) {
	                Pair<Integer,Double> p = new Pair<Integer,Double>(dist, new Double(words[1]));
	                if (p.cdr().doubleValue()>=0)	// should be non-negative value
	                	empiricalDistribution.add(p);
	                else {
	                	System.err.println("\nRead distribution file contains negative probability(count) value!"); 
	                	System.exit(1);
	                }
	              }
	            }
	        }
	        loadData(empiricalDistribution);
			makeProbabilities();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public BindingModel(File f) {
		this(f, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	
	public BindingModel(List<Pair<Integer, Double>> bindingDist){
		min=0; max=0;
		empiricalDistribution=bindingDist;
		loadData(bindingDist);
		makeProbabilities();
	}
	
	//Accessors
	public int getMin(){return min;}
	public int getMax(){return max;}
	public int getSummit(){return summit;}
	public int getInfluenceRange(){return influenceRange;}
	public double[] getProbabilities(){	return  probs.clone();}
	public double[] getLogProbabilities() { return logProbs.clone();}
	public String getFileName() {
		return fileName;
	}
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}	
	//Return a pair of distances corresponding to the central probability interval provided
	//Can be used to provide hit extension lengths
	protected Pair<Integer,Integer> probIntervalDistances(double prob){
		double ends=(1-prob)/2;
		double probSum=0;
		boolean firstFound=false, secondFound=false;
		int first=min, second=max;
		for(int i=min; i<=max; i++){
			probSum+=probability(i);
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
		influenceRange = longest*2;
	}
	//Load data
	protected void loadData(List<Pair<Integer, Double>> bindingDist){
		//Assumes the list is sorted//
		
		//Find max, min values first
		for(Pair<Integer, Double> p : bindingDist){
			if(p.car()<min)
				min=p.car();
			if(p.car()>max)
				max=p.car();
		}
		//Initialize arrays
		data = new double[(max-min)+1];
		probs = new double[(max-min)+1];
		logProbs = new double[(max-min)+1];
		for(int i=0; i<=(max-min); i++){
			data[i]=0; probs[i]=0; logProbs[i] = Double.NEGATIVE_INFINITY;
		}
		
		//Populate the data array (assumes sorted)
		int last=min-1;
		for(Pair<Integer, Double> p : bindingDist){
			int index = p.car();
			double val = p.cdr();
			//if list is not properly sorted (need to make this into an exception)
			if(index-last<0){
				System.err.println("Incorrectly sorted binding read distribution data!"); 
				System.exit(1);
			}
			//if unevenly spaced, smooth linearly between values
			if(index-last>1){
				double lastVal=dataVal(last);
				double step = (val-lastVal)/(double)(index-last);
				for(int i=1; i<(index-last); i++){
					data[(last+i)-min]=lastVal+(step*(double)i);
				}
			}
			data[index-min]=val;
			
			last = p.car();
		}
	}
	
	//Set a probability landscape according to the data. 
	protected void makeProbabilities(){
		double totalVal=0, minProb=Double.MAX_VALUE;
		for(int i=min; i<=max; i++){
			totalVal+=dataVal(i);
		}
		for(int i=min; i<=max; i++){
			probs[i-min] = dataVal(i)/totalVal; 
			logProbs[i-min] = Math.log(probs[i-min])/LOG2;
			if(probs[i-min]<minProb)
				minProb = probs[i-min];
		}
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
		
		// update empiricalDistribution with normalized probability
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>> ();
		for(int i=min; i<=max; i++){
			newDist.add(new Pair<Integer, Double>(i, probability(i)));
		}
		empiricalDistribution=newDist;
		bgProb = minProb/1000;
		logBgProb = Math.log(bgProb)/LOG2;

//		updateInfluenceRange();
	}
	
	protected void smooth(int splineStepSize, int avgStepSize){
		probs=StatUtil.cubicSpline(probs, splineStepSize, avgStepSize);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
	}
	protected void smoothGaussian (int kernelWidth){
		probs=StatUtil.gaussianSmoother(probs, kernelWidth);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
	}	
	//Look up the probability corresponding to a distance
	// Distance should be defined as (Read position - Peak position)
	public double probability(int distance){
		if(distance<min || distance>max){
			return(bgProb);
		}else{
			return(probs[distance-min]);
		}
	}
	
	public double logProbability(int distance) {
		if(distance<min || distance>max){
		  return(logBgProb);
		}else{
		  return(logProbs[distance-min]);
		}	  
	}
	
	//Look up the data corresponding to a distance
	protected double dataVal(int distance){
		if(distance<min || distance>max){
			return(0.0);
		}else{
			return(data[distance-min]);
		}
	}
		
	//Print probs to a file
	public void printToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int i=min; i<=max; i++){
				fout.write(i+"\t"+probability(i)+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Command-line interface to load a BindingModel from a file
	 * Example running: <br>
	 * <tt>
	 * --in oct4.shear.ext.txt --out out_oct4.txt --out_smooth out_oct4_smooth.txt
	 * </tt>
	 */
	public static void main(String[] args){
		if(Args.parseArgs(args).contains("in")){
			String infile = Args.parseString(args, "in", null);
			String outfile = Args.parseString(args, "out", "out.model");
			String outfile_smooth = Args.parseString(args, "out_smooth", "out_smooth.model");
			File pFile = new File(infile);
			if(!pFile.isFile()){
				System.err.println("Invalid file name");
				System.exit(1);
			}
	        //File loaded, make a BindingModel
	        BindingModel model = new BindingModel(pFile);
	        model.printToFile(outfile);
	        
	        model.smooth(SMOOTHING_STEPSIZE, SMOOTHING_AVG_PTS);
	        model.printToFile(outfile_smooth);
		}else{
			System.out.println("Usage: BindingModel --in GPSfileName --out outfile");
		}
	}
	
	public List<Pair<Integer, Double>> getEmpiricalDistribution() {
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>> ();
		for (Pair<Integer, Double> p: empiricalDistribution)
			newDist.add(p);
		return newDist;
	}
	
	/**
	 * Estimate a better read profile ranges
	 * @param minLeft minimum range on the left side
	 * @param minRight	minimum range on the right side
	 * @return The estimated profile ranges <left, right>
	 */
	public Pair<Integer, Integer> getNewEnds(int minLeft, int minRight){
		double halfHeight = probability(summit)*0.5;
		int leftHalfHeightEnd = 0;
		int rightHalfHeightEnd = 0;
		for (int i=min; i<=max; i++){
			if (probability(i)>=halfHeight){
				leftHalfHeightEnd = i;
				break;
			}
		}
		for (int i=max; i<=min; i--){
			if (probability(i)>=halfHeight){
				rightHalfHeightEnd = i;
				break;
			}
		}
		int left=Math.max(minLeft, Math.abs(summit-(summit-leftHalfHeightEnd)*4));
		int right=Math.max(minRight, Math.abs(summit+(rightHalfHeightEnd-summit)*3));

		return new Pair<Integer, Integer>(left, right);
	}
	
	// shift two array elements to give best KL-divergence
	// it will mutate the two input arrays.
	// assume a and b have same length
	public static int minKL_Shift (double[] a, double[] b){
		if (a.length!=b.length)
			return 0;
		int length = a.length;
		int range=20;
		double minKL=Double.MAX_VALUE;
		int shift=0;
		double[] a_new = new double[length];
		double[] b_new = new double[length];
		double[] a1 = new double[length];
		double[] b1 = new double[length];
		
		for (int i=-range;i<range;i++){
			if (i<0){				//a shift forward, b shift backward
				System.arraycopy(a, 0, a1, -i, length+i);
				System.arraycopy(a, length+i, a1, 0, -i);
				System.arraycopy(b, -i, b1, 0, length+i);
				System.arraycopy(b, 0, b1, length+i, -i);
			}
			else if (i>0){ 			//a shift backward, b shift forward
				System.arraycopy(a, i, a1, 0, length-i);
				System.arraycopy(a, 0, a1, length-i, i);
				System.arraycopy(b, 0, b1, i, length-i);
				System.arraycopy(b, length-i, b1, 0, i);
			}
			else{
				System.arraycopy(a, 0, a1, 0, length);
				System.arraycopy(b, 0, b1, 0, length);
			}

			double kl = StatUtil.log_KL_Divergence(a1, b1);

			if (kl<minKL){
				minKL=kl;
				shift=i;
				System.arraycopy(a1, 0, a_new, 0, length);
				System.arraycopy(b1, 0, b_new, 0, length);
			}
		}

		System.arraycopy(a_new, 0, a, 0, length);
		System.arraycopy(b_new, 0, b, 0, length);

		return shift;
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
	
	@SuppressWarnings("unchecked")
	public static final List<Pair<Integer, Double>> defaultChipExoEmpiricalDistribution =
			Arrays.asList(new Pair<Integer,Double>(-200, 0.00116031689418),	
					new Pair<Integer,Double>(-199, 0.00119151591863),
					new Pair<Integer,Double>(-198, 0.001219954197),
					new Pair<Integer,Double>(-197, 0.00124788861399),
					new Pair<Integer,Double>(-196, 0.00127436543682),
					new Pair<Integer,Double>(-195, 0.00129126146538),
					new Pair<Integer,Double>(-194, 0.0012917112882),
					new Pair<Integer,Double>(-193, 0.00128073300529),
					new Pair<Integer,Double>(-192, 0.00127301127302),
					new Pair<Integer,Double>(-191, 0.00127935938932),
					new Pair<Integer,Double>(-190, 0.00129847127724),
					new Pair<Integer,Double>(-189, 0.00132199724721),
					new Pair<Integer,Double>(-188, 0.00134344669845),
					new Pair<Integer,Double>(-187, 0.00136117684701),
					new Pair<Integer,Double>(-186, 0.00137589485105),
					new Pair<Integer,Double>(-185, 0.00138752222911),
					new Pair<Integer,Double>(-184, 0.00139530781234),
					new Pair<Integer,Double>(-183, 0.00140157273522),
					new Pair<Integer,Double>(-182, 0.00141305304162),
					new Pair<Integer,Double>(-181, 0.00143502366568),
					new Pair<Integer,Double>(-180, 0.00146420022137),
					new Pair<Integer,Double>(-179, 0.00148949279202),
					new Pair<Integer,Double>(-178, 0.00149937929142),
					new Pair<Integer,Double>(-177, 0.00148799194835),
					new Pair<Integer,Double>(-176, 0.00145667255413),
					new Pair<Integer,Double>(-175, 0.00141315862719),
					new Pair<Integer,Double>(-174, 0.00137009238531),
					new Pair<Integer,Double>(-173, 0.00134167350578),
					new Pair<Integer,Double>(-172, 0.0013377970071),
					new Pair<Integer,Double>(-171, 0.00136070657817),
					new Pair<Integer,Double>(-170, 0.00140756430311),
					new Pair<Integer,Double>(-169, 0.00147263643461),
					new Pair<Integer,Double>(-168, 0.00154379402499),
					new Pair<Integer,Double>(-167, 0.00160157638558),
					new Pair<Integer,Double>(-166, 0.00162884280434),
					new Pair<Integer,Double>(-165, 0.00162281804418),
					new Pair<Integer,Double>(-164, 0.00159606676918),
					new Pair<Integer,Double>(-163, 0.00156655624898),
					new Pair<Integer,Double>(-162, 0.00154694769488),
					new Pair<Integer,Double>(-161, 0.00153924065577),
					new Pair<Integer,Double>(-160, 0.00153610061521),
					new Pair<Integer,Double>(-159, 0.00152885500223),
					new Pair<Integer,Double>(-158, 0.00151712783066),
					new Pair<Integer,Double>(-157, 0.00151065758042),
					new Pair<Integer,Double>(-156, 0.00152088412944),
					new Pair<Integer,Double>(-155, 0.00155166824814),
					new Pair<Integer,Double>(-154, 0.00159788311171),
					new Pair<Integer,Double>(-153, 0.00164950879589),
					new Pair<Integer,Double>(-152, 0.00169515036787),
					new Pair<Integer,Double>(-151, 0.00172466580524),
					new Pair<Integer,Double>(-150, 0.00173199813793),
					new Pair<Integer,Double>(-149, 0.00171687370037),
					new Pair<Integer,Double>(-148, 0.00168511910616),
					new Pair<Integer,Double>(-147, 0.00164866237057),
					new Pair<Integer,Double>(-146, 0.00162277808231),
					new Pair<Integer,Double>(-145, 0.00161806449866),
					new Pair<Integer,Double>(-144, 0.00163281025886),
					new Pair<Integer,Double>(-143, 0.00165532997233),
					new Pair<Integer,Double>(-142, 0.00167434141666),
					new Pair<Integer,Double>(-141, 0.00168597592181),
					new Pair<Integer,Double>(-140, 0.00169254416726),
					new Pair<Integer,Double>(-139, 0.00169731636618),
					new Pair<Integer,Double>(-138, 0.00170061661825),
					new Pair<Integer,Double>(-137, 0.00169989881437),
					new Pair<Integer,Double>(-136, 0.00169443925576),
					new Pair<Integer,Double>(-135, 0.00169060068765),
					new Pair<Integer,Double>(-134, 0.00170017189377),
					new Pair<Integer,Double>(-133, 0.00172950471614),
					new Pair<Integer,Double>(-132, 0.00176965594705),
					new Pair<Integer,Double>(-131, 0.00180010249419),
					new Pair<Integer,Double>(-130, 0.00180348715929),
					new Pair<Integer,Double>(-129, 0.00177719533385),
					new Pair<Integer,Double>(-128, 0.00173404790006),
					new Pair<Integer,Double>(-127, 0.00169524293065),
					new Pair<Integer,Double>(-126, 0.00167952657457),
					new Pair<Integer,Double>(-125, 0.00169306782964),
					new Pair<Integer,Double>(-124, 0.00172752254996),
					new Pair<Integer,Double>(-123, 0.00176795914712),
					new Pair<Integer,Double>(-122, 0.00180187122305),
					new Pair<Integer,Double>(-121, 0.00182233310741),
					new Pair<Integer,Double>(-120, 0.00182778114568),
					new Pair<Integer,Double>(-119, 0.00182044442677),
					new Pair<Integer,Double>(-118, 0.00180300713201),
					new Pair<Integer,Double>(-117, 0.001777003049),
					new Pair<Integer,Double>(-116, 0.00174583822668),
					new Pair<Integer,Double>(-115, 0.00171682027254),
					new Pair<Integer,Double>(-114, 0.00169799403607),
					new Pair<Integer,Double>(-113, 0.00169494481843),
					new Pair<Integer,Double>(-112, 0.00171075252595),
					new Pair<Integer,Double>(-111, 0.00174455939304),
					new Pair<Integer,Double>(-110, 0.00178862858883),
					new Pair<Integer,Double>(-109, 0.00182985998943),
					new Pair<Integer,Double>(-108, 0.00185703488608),
					new Pair<Integer,Double>(-107, 0.00186793215748),
					new Pair<Integer,Double>(-106, 0.00186934691335),
					new Pair<Integer,Double>(-105, 0.0018698215687),
					new Pair<Integer,Double>(-104, 0.00187393887728),
					new Pair<Integer,Double>(-103, 0.00188370679424),
					new Pair<Integer,Double>(-102, 0.00190083554127),
					new Pair<Integer,Double>(-101, 0.00192384307933),
					new Pair<Integer,Double>(-100, 0.00194490921407),
					new Pair<Integer,Double>(-99, 0.00195412679664),
					new Pair<Integer,Double>(-98, 0.00194869817062),
					new Pair<Integer,Double>(-97, 0.00193685906307),
					new Pair<Integer,Double>(-96, 0.00193258265662),
					new Pair<Integer,Double>(-95, 0.00194616183108),
					new Pair<Integer,Double>(-94, 0.00197702622518),
					new Pair<Integer,Double>(-93, 0.0020133863262),
					new Pair<Integer,Double>(-92, 0.00204062107571),
					new Pair<Integer,Double>(-91, 0.00205272721716),
					new Pair<Integer,Double>(-90, 0.00205623339718),
					new Pair<Integer,Double>(-89, 0.00206240452274),
					new Pair<Integer,Double>(-88, 0.00207581619785),
					new Pair<Integer,Double>(-87, 0.00209128020082),
					new Pair<Integer,Double>(-86, 0.00210162180429),
					new Pair<Integer,Double>(-85, 0.00210694694509),
					new Pair<Integer,Double>(-84, 0.00211510158348),
					new Pair<Integer,Double>(-83, 0.00213328912338),
					new Pair<Integer,Double>(-82, 0.00216016473463),
					new Pair<Integer,Double>(-81, 0.00218637952393),
					new Pair<Integer,Double>(-80, 0.00220210750389),
					new Pair<Integer,Double>(-79, 0.00220357039255),
					new Pair<Integer,Double>(-78, 0.0021945325703),
					new Pair<Integer,Double>(-77, 0.0021851768327),
					new Pair<Integer,Double>(-76, 0.00218901625862),
					new Pair<Integer,Double>(-75, 0.0022164934691),
					new Pair<Integer,Double>(-74, 0.00226786831296),
					new Pair<Integer,Double>(-73, 0.00233202339507),
					new Pair<Integer,Double>(-72, 0.00239425655393),
					new Pair<Integer,Double>(-71, 0.00244684398628),
					new Pair<Integer,Double>(-70, 0.00249084545895),
					new Pair<Integer,Double>(-69, 0.00252643559345),
					new Pair<Integer,Double>(-68, 0.0025458399497),
					new Pair<Integer,Double>(-67, 0.00254168470468),
					new Pair<Integer,Double>(-66, 0.00251986896609),
					new Pair<Integer,Double>(-65, 0.00249693771471),
					new Pair<Integer,Double>(-64, 0.00248432541995),
					new Pair<Integer,Double>(-63, 0.00248151448345),
					new Pair<Integer,Double>(-62, 0.00248664781232),
					new Pair<Integer,Double>(-61, 0.00250617670694),
					new Pair<Integer,Double>(-60, 0.00254753107498),
					new Pair<Integer,Double>(-59, 0.00260638701136),
					new Pair<Integer,Double>(-58, 0.00266743916197),
					new Pair<Integer,Double>(-57, 0.00271730343174),
					new Pair<Integer,Double>(-56, 0.00275523681205),
					new Pair<Integer,Double>(-55, 0.00279376667421),
					new Pair<Integer,Double>(-54, 0.00284959654668),
					new Pair<Integer,Double>(-53, 0.00292996172514),
					new Pair<Integer,Double>(-52, 0.00302562433661),
					new Pair<Integer,Double>(-51, 0.00311951327157),
					new Pair<Integer,Double>(-50, 0.00320314327497),
					new Pair<Integer,Double>(-49, 0.00328411471184),
					new Pair<Integer,Double>(-48, 0.00337934830131),
					new Pair<Integer,Double>(-47, 0.00350124436305),
					new Pair<Integer,Double>(-46, 0.00364616264291),
					new Pair<Integer,Double>(-45, 0.00379225289039),
					new Pair<Integer,Double>(-44, 0.00390853237632),
					new Pair<Integer,Double>(-43, 0.00397164277527),
					new Pair<Integer,Double>(-42, 0.0039819438024),
					new Pair<Integer,Double>(-41, 0.00396821662805),
					new Pair<Integer,Double>(-40, 0.00397591961098),
					new Pair<Integer,Double>(-39, 0.00404489792006),
					new Pair<Integer,Double>(-38, 0.00418892975502),
					new Pair<Integer,Double>(-37, 0.00438727061344),
					new Pair<Integer,Double>(-36, 0.00459178995926),
					new Pair<Integer,Double>(-35, 0.00474805042582),
					new Pair<Integer,Double>(-34, 0.00482299909396),
					new Pair<Integer,Double>(-33, 0.00482197205735),
					new Pair<Integer,Double>(-32, 0.00477715522865),
					new Pair<Integer,Double>(-31, 0.00471598310237),
					new Pair<Integer,Double>(-30, 0.00464488514386),
					new Pair<Integer,Double>(-29, 0.00456320488895),
					new Pair<Integer,Double>(-28, 0.00447927395254),
					new Pair<Integer,Double>(-27, 0.00440292175835),
					new Pair<Integer,Double>(-26, 0.00433024963587),
					new Pair<Integer,Double>(-25, 0.00425012913431),
					new Pair<Integer,Double>(-24, 0.00416961789072),
					new Pair<Integer,Double>(-23, 0.00412845024798),
					new Pair<Integer,Double>(-22, 0.00418811326623),
					new Pair<Integer,Double>(-21, 0.00443517031822),
					new Pair<Integer,Double>(-20, 0.00503886752527),
					new Pair<Integer,Double>(-19, 0.00633163233456),
					new Pair<Integer,Double>(-18, 0.00873393948363),
					new Pair<Integer,Double>(-17, 0.0123386092788),
					new Pair<Integer,Double>(-16, 0.0163835237051),
					new Pair<Integer,Double>(-15, 0.0193440544473),
					new Pair<Integer,Double>(-14, 0.0199546431357),
					new Pair<Integer,Double>(-13, 0.01823172755),
					new Pair<Integer,Double>(-12, 0.0153520829615),
					new Pair<Integer,Double>(-11, 0.0126098117635),
					new Pair<Integer,Double>(-10, 0.0106240429004),
					new Pair<Integer,Double>(-9, 0.00933359384468),
					new Pair<Integer,Double>(-8, 0.00839300452524),
					new Pair<Integer,Double>(-7, 0.0075047478554),
					new Pair<Integer,Double>(-6, 0.00656499302384),
					new Pair<Integer,Double>(-5, 0.00565868664955),
					new Pair<Integer,Double>(-4, 0.0049377646129),
					new Pair<Integer,Double>(-3, 0.00447446059104),
					new Pair<Integer,Double>(-2, 0.00421546620153),
					new Pair<Integer,Double>(-1, 0.00407851902705),
					new Pair<Integer,Double>(0, 0.00407751680367),
					new Pair<Integer,Double>(1, 0.00431537925517),
					new Pair<Integer,Double>(2, 0.00482890217969),
					new Pair<Integer,Double>(3, 0.0054679327004),
					new Pair<Integer,Double>(4, 0.00597158540991),
					new Pair<Integer,Double>(5, 0.00616844942732),
					new Pair<Integer,Double>(6, 0.00608696734977),
					new Pair<Integer,Double>(7, 0.00588635712449),
					new Pair<Integer,Double>(8, 0.00571036125885),
					new Pair<Integer,Double>(9, 0.00560207064173),
					new Pair<Integer,Double>(10, 0.00552147717793),
					new Pair<Integer,Double>(11, 0.00541583015905),
					new Pair<Integer,Double>(12, 0.00526469328199),
					new Pair<Integer,Double>(13, 0.00506780859201),
					new Pair<Integer,Double>(14, 0.0048179089131),
					new Pair<Integer,Double>(15, 0.00451020383195),
					new Pair<Integer,Double>(16, 0.0041727525812),
					new Pair<Integer,Double>(17, 0.00386112795631),
					new Pair<Integer,Double>(18, 0.00361577272677),
					new Pair<Integer,Double>(19, 0.00343835120674),
					new Pair<Integer,Double>(20, 0.00331135026949),
					new Pair<Integer,Double>(21, 0.00322258353327),
					new Pair<Integer,Double>(22, 0.00316191721415),
					new Pair<Integer,Double>(23, 0.00311036159731),
					new Pair<Integer,Double>(24, 0.00304800702894),
					new Pair<Integer,Double>(25, 0.00296858143756),
					new Pair<Integer,Double>(26, 0.00287781656025),
					new Pair<Integer,Double>(27, 0.00278001487406),
					new Pair<Integer,Double>(28, 0.00267093352497),
					new Pair<Integer,Double>(29, 0.00254347654296),
					new Pair<Integer,Double>(30, 0.00239934497412),
					new Pair<Integer,Double>(31, 0.00225483339105),
					new Pair<Integer,Double>(32, 0.00213278885859),
					new Pair<Integer,Double>(33, 0.00204620740924),
					new Pair<Integer,Double>(34, 0.0019895468377),
					new Pair<Integer,Double>(35, 0.00194494324095),
					new Pair<Integer,Double>(36, 0.00189534017747),
					new Pair<Integer,Double>(37, 0.00183454875222),
					new Pair<Integer,Double>(38, 0.0017687928489),
					new Pair<Integer,Double>(39, 0.00170863302602),
					new Pair<Integer,Double>(40, 0.00165910958786),
					new Pair<Integer,Double>(41, 0.00161870497332),
					new Pair<Integer,Double>(42, 0.001585810396),
					new Pair<Integer,Double>(43, 0.00156258580514),
					new Pair<Integer,Double>(44, 0.0015515789543),
					new Pair<Integer,Double>(45, 0.00154961214787),
					new Pair<Integer,Double>(46, 0.0015458222008),
					new Pair<Integer,Double>(47, 0.00152807537877),
					new Pair<Integer,Double>(48, 0.00149363632243),
					new Pair<Integer,Double>(49, 0.00145245661458),
					new Pair<Integer,Double>(50, 0.00141895655599),
					new Pair<Integer,Double>(51, 0.00140183443011),
					new Pair<Integer,Double>(52, 0.00140018846931),
					new Pair<Integer,Double>(53, 0.00140521397247),
					new Pair<Integer,Double>(54, 0.00140490757726),
					new Pair<Integer,Double>(55, 0.00139028613963),
					new Pair<Integer,Double>(56, 0.00135992112851),
					new Pair<Integer,Double>(57, 0.00131977609019),
					new Pair<Integer,Double>(58, 0.00127926291664),
					new Pair<Integer,Double>(59, 0.00124679737381),
					new Pair<Integer,Double>(60, 0.00122762817726),
					new Pair<Integer,Double>(61, 0.00122404603908),
					new Pair<Integer,Double>(62, 0.00123622294076),
					new Pair<Integer,Double>(63, 0.00126229851951),
					new Pair<Integer,Double>(64, 0.00129787793024),
					new Pair<Integer,Double>(65, 0.00133591835365),
					new Pair<Integer,Double>(66, 0.00136832430709),
					new Pair<Integer,Double>(67, 0.00138941451742),
					new Pair<Integer,Double>(68, 0.00139813649626),
					new Pair<Integer,Double>(69, 0.00139718938285),
					new Pair<Integer,Double>(70, 0.00139202796565),
					new Pair<Integer,Double>(71, 0.00139054084133),
					new Pair<Integer,Double>(72, 0.00139920617572),
					new Pair<Integer,Double>(73, 0.00141646844133),
					new Pair<Integer,Double>(74, 0.00143175953433),
					new Pair<Integer,Double>(75, 0.00143315358276),
					new Pair<Integer,Double>(76, 0.00141531566658),
					new Pair<Integer,Double>(77, 0.00138013578264),
					new Pair<Integer,Double>(78, 0.00133283708376),
					new Pair<Integer,Double>(79, 0.00128049687629),
					new Pair<Integer,Double>(80, 0.0012328397454),
					new Pair<Integer,Double>(81, 0.00119995997),
					new Pair<Integer,Double>(82, 0.00118736177124),
					new Pair<Integer,Double>(83, 0.00119325461317),
					new Pair<Integer,Double>(84, 0.00120988286496),
					new Pair<Integer,Double>(85, 0.00122774745206),
					new Pair<Integer,Double>(86, 0.00124067962182),
					new Pair<Integer,Double>(87, 0.00124823513885),
					new Pair<Integer,Double>(88, 0.00125338324981),
					new Pair<Integer,Double>(89, 0.00125863218092),
					new Pair<Integer,Double>(90, 0.00126458458157),
					new Pair<Integer,Double>(91, 0.00126960875233),
					new Pair<Integer,Double>(92, 0.0012685220946),
					new Pair<Integer,Double>(93, 0.00125403873688),
					new Pair<Integer,Double>(94, 0.00122367038818),
					new Pair<Integer,Double>(95, 0.00118542486732),
					new Pair<Integer,Double>(96, 0.00115381843117),
					new Pair<Integer,Double>(97, 0.00113949451426),
					new Pair<Integer,Double>(98, 0.00114375006046),
					new Pair<Integer,Double>(99, 0.00116230362561),
					new Pair<Integer,Double>(100, 0.00119132969052),
					new Pair<Integer,Double>(101, 0.00122817613346),
					new Pair<Integer,Double>(102, 0.00126777686085),
					new Pair<Integer,Double>(103, 0.00130270248284),
					new Pair<Integer,Double>(104, 0.00132872433846),
					new Pair<Integer,Double>(105, 0.00134686272092),
					new Pair<Integer,Double>(106, 0.0013566775772),
					new Pair<Integer,Double>(107, 0.00135134873901),
					new Pair<Integer,Double>(108, 0.00132516217387),
					new Pair<Integer,Double>(109, 0.00128462688486),
					new Pair<Integer,Double>(110, 0.00124696027413),
					new Pair<Integer,Double>(111, 0.00122686531615),
					new Pair<Integer,Double>(112, 0.00122639543019),
					new Pair<Integer,Double>(113, 0.00123585058374),
					new Pair<Integer,Double>(114, 0.00124225350819),
					new Pair<Integer,Double>(115, 0.00123862316176),
					new Pair<Integer,Double>(116, 0.00122794325123),
					new Pair<Integer,Double>(117, 0.00121957926962),
					new Pair<Integer,Double>(118, 0.00122192141153),
					new Pair<Integer,Double>(119, 0.00123733636032),
					new Pair<Integer,Double>(120, 0.00126093916681),
					new Pair<Integer,Double>(121, 0.00128156787006),
					new Pair<Integer,Double>(122, 0.00128580294311),
					new Pair<Integer,Double>(123, 0.00126626948304),
					new Pair<Integer,Double>(124, 0.0012276025824),
					new Pair<Integer,Double>(125, 0.00118229133218),
					new Pair<Integer,Double>(126, 0.00114046879258),
					new Pair<Integer,Double>(127, 0.00110408707044),
					new Pair<Integer,Double>(128, 0.00106909881939),
					new Pair<Integer,Double>(129, 0.00103149482918),
					new Pair<Integer,Double>(130, 0.000991435673673),
					new Pair<Integer,Double>(131, 0.000953091072761),
					new Pair<Integer,Double>(132, 0.00092229921941),
					new Pair<Integer,Double>(133, 0.000904176770421),
					new Pair<Integer,Double>(134, 0.000899801174527),
					new Pair<Integer,Double>(135, 0.00090327548629),
					new Pair<Integer,Double>(136, 0.000903990327694),
					new Pair<Integer,Double>(137, 0.000894176777689),
					new Pair<Integer,Double>(138, 0.000874570190291),
					new Pair<Integer,Double>(139, 0.000851963811615),
					new Pair<Integer,Double>(140, 0.000831496315642),
					new Pair<Integer,Double>(141, 0.000813066028895),
					new Pair<Integer,Double>(142, 0.000795356041162),
					new Pair<Integer,Double>(143, 0.000780076183475),
					new Pair<Integer,Double>(144, 0.000769870621891),
					new Pair<Integer,Double>(145, 0.000763787273733),
					new Pair<Integer,Double>(146, 0.000757377996392),
					new Pair<Integer,Double>(147, 0.000747143009951),
					new Pair<Integer,Double>(148, 0.000733611759392),
					new Pair<Integer,Double>(149, 0.00072046442412),
					new Pair<Integer,Double>(150, 0.000711439294768),
					new Pair<Integer,Double>(151, 0.000707632642357),
					new Pair<Integer,Double>(152, 0.000707236704062),
					new Pair<Integer,Double>(153, 0.000707488921501),
					new Pair<Integer,Double>(154, 0.000706738308035),
					new Pair<Integer,Double>(155, 0.000705180337753),
					new Pair<Integer,Double>(156, 0.000704129682437),
					new Pair<Integer,Double>(157, 0.000704947601082),
					new Pair<Integer,Double>(158, 0.000708822610843),
					new Pair<Integer,Double>(159, 0.000716161071127),
					new Pair<Integer,Double>(160, 0.00072446072585),
					new Pair<Integer,Double>(161, 0.00072798393636),
					new Pair<Integer,Double>(162, 0.000721912835448),
					new Pair<Integer,Double>(163, 0.000706644737427),
					new Pair<Integer,Double>(164, 0.000686969716968),
					new Pair<Integer,Double>(165, 0.000668897288883),
					new Pair<Integer,Double>(166, 0.000659095227796),
					new Pair<Integer,Double>(167, 0.000664646363553),
					new Pair<Integer,Double>(168, 0.000689354225131),
					new Pair<Integer,Double>(169, 0.000729203575363),
					new Pair<Integer,Double>(170, 0.000772238967628),
					new Pair<Integer,Double>(171, 0.00080427885623),
					new Pair<Integer,Double>(172, 0.0008171841262),
					new Pair<Integer,Double>(173, 0.000813611636715),
					new Pair<Integer,Double>(174, 0.000803631185639),
					new Pair<Integer,Double>(175, 0.000795655406089),
					new Pair<Integer,Double>(176, 0.000790036227775),
					new Pair<Integer,Double>(177, 0.000780640554249),
					new Pair<Integer,Double>(178, 0.000761287970143),
					new Pair<Integer,Double>(179, 0.000731078359128),
					new Pair<Integer,Double>(180, 0.000695847909194),
					new Pair<Integer,Double>(181, 0.000665115417746),
					new Pair<Integer,Double>(182, 0.000646048243478),
					new Pair<Integer,Double>(183, 0.000639375123329),
					new Pair<Integer,Double>(184, 0.000640785846385),
					new Pair<Integer,Double>(185, 0.000645445564149),
					new Pair<Integer,Double>(186, 0.000651417483736),
					new Pair<Integer,Double>(187, 0.000659712573423),
					new Pair<Integer,Double>(188, 0.000671087716982),
					new Pair<Integer,Double>(189, 0.000683208018992),
					new Pair<Integer,Double>(190, 0.000691804732176),
					new Pair<Integer,Double>(191, 0.000694296817562),
					new Pair<Integer,Double>(192, 0.000691527553615),
					new Pair<Integer,Double>(193, 0.000686403633218),
					new Pair<Integer,Double>(194, 0.000681460298531),
					new Pair<Integer,Double>(195, 0.000677553819079),
					new Pair<Integer,Double>(196, 0.000674095330932),
					new Pair<Integer,Double>(197, 0.000669432986286),
					new Pair<Integer,Double>(198, 0.000660936221991),
					new Pair<Integer,Double>(199, 0.000646672152263),
					new Pair<Integer,Double>(200, 0.000628546124895));
}