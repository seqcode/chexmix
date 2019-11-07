package org.seqcode.projects.chexmix.framework;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TimeZone;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.io.IOUtil;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.motifs.FreqMatrixImport;

/**
 * ChExMixConfig: 
 * 		Maintains all constants needed by ChExMix. 
 *     
 * @author Naomi Yamada
 * @version	%I%, %G%
 */
public class ChExMixConfig {

	public static String version = "0.42";
	public boolean isGPS=true;
	protected GenomeConfig gconfig;
	protected Genome gen=null;
	protected String outName="chexmix", outBase="chexmix";
	protected File outDir=null, interDir=null, imagesDir=null;
	protected boolean printHelp=false;
	protected double sigLogConf=-7; 
	protected double prLogConf=-6; 
	protected int minModelUpdateRounds=1; //Minimum number of outer EM training rounds
	protected int maxModelUpdateRounds=3; //Maximum number of outer EM training rounds
	protected int posPriorScaling=10;
	protected double modelConvergenceKL=-25; //KL-divergence threshold for convergence 
	protected int maxThreads=1;				//Number of threads to use. Default is 1 for single processor machines. 
	protected double alphaScalingFactor = 1.0; //Scale the condition-specific alpha value by this factor
	protected double fixedAlpha = 0.0; //Fixed alpha value if above 0
	protected double betaScalingFactor = 0.05; //Scale the condition and component-specific beta value by this factor
	protected double epsilonScalingFactor = 0.2; //Scale the condition and component-specific epsilon value by this factor
	protected double motifMinROC = 0.70; //Motif prior is used only if the ROC is greater than this .	
	protected double extendWindow = 500; //Range extension around gff points
	protected int bmAnalysisWindowMax=2000;
	protected int minComponentsForBMUpdate = 50;
	protected int minRefsForBMUpdate = 25;
	protected double minSubtypeFraction = 0.05; // A subtype needs to be associated with at least this fraction of binding events to be supported 
	protected double minComponentReadFactorForBM = 3; //Components must have (this factor times the condition alpha) number of reads assigned before being included in BM update
	protected boolean updateBM=true; //Set to false to turn off binding model update
    protected double gauss_smooth = 1; //Variance for Gaussian smoothing 
    protected int addFlankingComponentSpacing=20; //In non-first rounds of EM, the components are initialized using the positions from the last round with additional flanking components added at this spacing
	protected boolean addFlankingComponent=true;
    protected List<Region> regionsToPlot = new ArrayList<Region>(); //List of regions that will be printed during EM training (for debugging/demonstration)
	protected List<Region> regionsToIgnore = new ArrayList<Region>(); //List of regions that will be ignored during EM training (i.e. known towers, etc)
	protected List<Point> initialPos=null;	//List of points loaded from peak file and used to place binding components.
	protected List<StrandedRegion> motifRegions = new ArrayList<StrandedRegion>(); //List of regions to construct cross-linking point histograms (for testing)
	protected boolean findMotifs = true; //Run motif-finding for motif prior
	protected boolean motif_posprior=true; //You can have motif-finding without using the motif-prior
	protected boolean clusterReadDistributions = true; //Cluster read distributions
	protected boolean doReadFilter=false;	// Turn on per base read filter in case of highly duplicated experiment
	protected String MEMEpath="";
	protected String MEMEargs=" -dna -mod zoops -revcomp -nostatus ";    //Do not try using -p here; it leads to MEME runtime errors
	protected boolean MEMEnonparallel=false; //flag to enforce use of non-parallel version
	public int MEMEminw=7;
	public int MEMEmaxw=18;
	protected MarkovBackgroundModel markovBackMode; // Markov background model file
	protected boolean verbose = false; //Print extra output
	protected List<List<StrandedPoint>> initialClustPoints = null; // Initial cluster points
	protected String MetaMakerArgs="--win 250 --bins 250 --batch --nocolorbar --linemax 3 --linethick 1 --noborder";
	protected int initComponentSpacing=30;	//Initial component spacing
	protected String distribA=null;	//Stranded distribution A
	protected String distribB=null;	//Stranded distribution B
	protected List<WeightMatrix> initMotifs=null; //Initial motifs
	protected double preferenceValue = -0.1; // Preference value for read distribution clustering
	protected int clusteringWindow = 150;
	protected int numClusteringComps = 500;	// Number of components to perform AP clustering
	protected double MarkovBackSeqRmThres = 0.1; // Markov background threshold for removing sequences
	protected int modelRange = 50;	// Window size to extract tag counts for potential region filter, Poisson significance test, and other purposes
	protected boolean galaxyhtml=false; // Output simpler html file for galaxy 
	protected boolean shareSubtypes=true;	// Share subtypes across experiments
	protected boolean standardMode=true; // Mode in which events are reported if summed tag counts significant over background in condition as a whole.
	protected boolean lenientMode=false; // Mode in which events are reported if significant over background in >=1 replicate *or* the condition as a whole.  
	protected boolean lenientPlusMode=false; // Mode in which events are reported if summed tag counts significant over background in condition as a whole OR (significant over background in >=1 replicate AND no significant difference between replicates).  
	
	
	//Constants
	public final double LOG2 = Math.log(2);
	public final int POTREG_BIN_STEP = 100; //Sliding window step in potential region scanner
	public final int MAXSECTION = 50000000;
	public final double INIT_CS_TO_XL_RATIO=0.05; 	//Initial ratio of CS component pi values to sum of XO pi values.
	public final double MIN_CS_PI = 0.05; //Minimum pi value for CS component
	public final double MIN_ALPHA = 0.01; //Minimum alpha 
	public final boolean FIXED_XL_OFFSET=true; //Estimate the XL component offset (mean)?
	public final boolean XL_DISTRIB_SYMMETRIC=true; //Is the sigma associated with XL components symmetric?
//    public final int INIT_COMPONENT_SPACING=30;  //Initial component spacing
	public final int MAX_EM_ITER=2000;
    public final int EM_ML_ITER=50;     				//Run EM up until <tt>ML_ITER</tt> without using sparse prior	**Changed from 100 to 50
    public final int ML_ML_ITER=50;     				//Run ML up until <tt>ML_ITER</tt> without using sparse prior	**Changed from 100 to 50
    public final int ALPHA_ANNEALING_ITER=100;     //Run EM up until <tt>ALPHA_ANNEALING_ITER</tt> with smaller alpha based on the current iteration 
    public final int POSPRIOR_ITER=150;     //Run EM up until <tt>ALPHA_ANNEALING_ITER</tt> with uniform positional prior and then up until at least <tt>POSPRIOR_ANNEALING_ITER</tt> with activated positional prior
    public final int EM_MU_UPDATE_WIN=30; //Half the window size in which to look for mu maximization (i.e. component position) during EM.
    public final double EM_CONVERGENCE = 1e-10; //EM convergence between the likelihood of the current and the previous step
    public final double EM_STATE_EQUIV_THRES = 1e-10; //EM state equivalence threshold 
    public final int EM_STATE_EQUIV_ROUNDS = 3; //Number of training rounds where the EM states have to be equivalent
    public final double NOISE_EMISSION_MIN = 0.05; //Arbitrary floor on the emission probability of noise (must be non-zero to mop up noise reads)
    public final double NOISE_EMISSION_MAX = 0.75; //Arbitrary ceiling on the emission probability of noise
    public final int NOISE_DISTRIB_SMOOTHING_WIN = 50; //Smoothing window for the noise distribution used in the BindingMixture
    public final int MAX_BINDINGMODEL_WIDTH=500; //Maximum width for binding models (affects how large the read profiles are for binding components)
	public final boolean CALC_LL=false; //Calculate the log-likelihood during EM.
	public final boolean CALC_COMP_LL=false; //Calculate component-wise log-likelihoods during ML
    public final int MOTIF_FINDING_SEQWINDOW=60; //Bases to extract around components for motif-finding			public final boolean CALC_LL=false; //Calculate the log-likelihood during EM. 
    public final int MOTIF_FINDING_LOCAL_SEQWINDOW=20; //Bases to extract around components for focused motif-finding
    public final int MOTIF_FINDING_TOPSEQS=1000; //Number of top components to analyze			public final boolean CALC_COMP_LL=false; //Calculate component-wise log-likelihoods during ML 
    public final double MOTIF_FINDING_ALLOWED_REPETITIVE = 0.2; //Percentage of the examined sequence window allowed to be lowercase or N			
    public final int MOTIF_FINDING_NEGSEQ=5000; //Number of negative sequences for motif significance tests		
    public final double MARKOV_BACK_MODEL_THRES = 0.05; // Markov background threshold for making models
    public final int SLIDING_WINDOW=60; // Sliding window range in computing KL divergence 
    public final double MOTIF_PCC_THRES = 0.95; //Motif length adjusted similarity threshold for selecting one motif
    public final int MARKOV_NUM_TEST=100000;
    public final int KMEANS_TRAIN_REPEATS=10;
    public final int KMEANS_MAX_ITER = 100;
    public final int KMEANS_K = 2;
    public final double KMEANS_CONVERGENCE_THRES = 0.01;
    public final double KL_DIVERGENCE_BM_THRES = -10;
    public final double CORR_THRES =100; //100 is no threshold
    public final long RANDOMSEED = 1000; //setting the random seed for the sake of reproducibility
	
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public ChExMixConfig(GenomeConfig gcon, String[] arguments){this(gcon, arguments, true);}
	public ChExMixConfig(GenomeConfig gcon, String[] arguments, boolean isGPS){
		System.setProperty("java.awt.headless", "true");
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args=arguments; 
		this.isGPS=isGPS;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			printHelp=true;			
		}else{
			try{
				//Test for a config file... if there is concatenate the contents into the args
				if(ap.hasKey("config")){
					ArrayList<String> confArgs = new ArrayList<String>();
					String confName = ap.getKeyValue("config");
					File confFile = new File(confName);
					if(!confFile.isFile())
						System.err.println("\nCannot find configuration file: "+confName);
					BufferedReader reader = new BufferedReader(new FileReader(confFile));
				    String line;
			        while ((line = reader.readLine()) != null) {
			        	line = line.trim();
			        	String[] words = line.split("\\s+");
			        	if(!words[0].startsWith("--"))
			        		words[0] = new String("--"+words[0]);
			        	confArgs.add(words[0]); 
			        	if(words.length>1){
				        	String rest=words[1];
				        	for(int w=2; w<words.length; w++)
				        		rest = rest+" "+words[w];
				        	confArgs.add(rest);
			        	}
			        }
			        String [] confArgsArr = confArgs.toArray(new String[confArgs.size()]);
			        String [] newargs =new String[args.length + confArgsArr.length];
			        System.arraycopy(args, 0, newargs, 0, args.length);
			        System.arraycopy(confArgsArr, 0, newargs, args.length, confArgsArr.length);
			        args = newargs;
			        ap = new ArgParser(args);
				}
				
				/****Miscellaneous arguments****/
				//Maximum number of model update rounds
				maxModelUpdateRounds = Args.parseInteger(args,"round", maxModelUpdateRounds);
				//Turn off binding model updates
				updateBM = Args.parseFlags(args).contains("nomodelupdate") ? false : true;
				//Minimum number of components to support a binding model update		
				minComponentsForBMUpdate = Args.parseInteger(args,"minmodelupdateevents",minComponentsForBMUpdate);
				//Minimum number of motif references  to support a binding model update		
				minRefsForBMUpdate = Args.parseInteger(args,"minmodelupdaterefs",minRefsForBMUpdate);
				//Parameter for Gaussian smoothing (std. dev.)
				gauss_smooth = Args.parseDouble(args,"gausssmoothparam",gauss_smooth);
				//Output path
				DateFormat df = new SimpleDateFormat("yyyy-MM-dd-hh-mm-ss");  
			    df.setTimeZone(TimeZone.getTimeZone("EST"));
				outName = Args.parseString(args, "out", outName+"_"+df.format(new Date()));
				outDir =  new File(outName); //Output directory
				outBase = outDir.getName(); //Last part of name
				//Background model parameters		
				sigLogConf = Args.parseDouble(args,"highlogconf",sigLogConf);		
				prLogConf = Args.parseDouble(args,"prlogconf",prLogConf);
				//Threads
				maxThreads = Args.parseInteger(args,"threads",maxThreads);
				maxThreads = Math.min(maxThreads, java.lang.Runtime.getRuntime().availableProcessors());
				//Alpha scaling factor
				alphaScalingFactor = Args.parseDouble(args,"alphascale",alphaScalingFactor);
				//Fixed alpha value
				fixedAlpha = Args.parseDouble(args,"fixedalpha",fixedAlpha);
				//Beta scaling factor
				betaScalingFactor = Args.parseDouble(args,"betascale",betaScalingFactor);
				//Epsilon scaling factor
				epsilonScalingFactor = Args.parseDouble(args,"epsilonscale",epsilonScalingFactor);
				//Motif prior is used only if the ROC is greater than this .
				motifMinROC = Args.parseDouble(args, "minroc", motifMinROC);
				//Number of base pair to extend around gff
				extendWindow = Args.parseDouble(args, "extwin", extendWindow);
				//Initial component spacing
				initComponentSpacing = Args.parseInteger(args,"compspacing",initComponentSpacing);
				//Preference value for AP clustering
				preferenceValue = Args.parseDouble(args, "pref", preferenceValue);
				//Number of components to perform AP clustering
				numClusteringComps = Args.parseInteger(args,"numcomps",numClusteringComps);
				//Window size for extracting tag counts
				modelRange = Args.parseInteger(args,"mrange",modelRange);
				//Max window size for running a mixture model over binding events
				System.out.println("before: "+bmAnalysisWindowMax);
				bmAnalysisWindowMax = Args.parseInteger(args,"bmwindowmax",bmAnalysisWindowMax);
				System.out.println("before: "+bmAnalysisWindowMax);
				
				if(ap.hasKey("plotregions"))
					regionsToPlot = RegionFileUtilities.loadRegionsFromFile(Args.parseString(args, "plotregions", null), gen, -1);
				//Regions to ignore during EM training
				if(ap.hasKey("exclude"))
					regionsToIgnore = RegionFileUtilities.loadRegionsFromFile(Args.parseString(args, "exclude", null), gen, -1);
				else if (ap.hasKey("excludebed"))
					regionsToIgnore = RegionFileUtilities.loadRegionsFromBEDFile(gen, Args.parseString(args, "excludebed", null), -1);
				//Initial peak file
				if (ap.hasKey("peakf"))
					initialPos = RegionFileUtilities.loadPeaksFromPeakFile(gen, Args.parseString(args, "peakf", null));
				//Motif for plotting components
				if (ap.hasKey("motifregions"))
					motifRegions = RegionFileUtilities.loadStrandedRegionsFromMotifFile(gen, Args.parseString(args, "motifregions", null), -1);

				//Turn off motif-finding 
				findMotifs = Args.parseFlags(args).contains("nomotifs") ? false : true;
				//Turn off motif prior only
				motif_posprior = (findMotifs && Args.parseFlags(args).contains("nomotifprior")) ? false : true;				
				//Check whether sequence is available (affects motif-finding behavior)
				if(isGPS && !gconfig.getSequenceGenerator().usingLocalFiles()){
					findMotifs=false;
					motif_posprior=false;
					System.err.println("No genome sequence data was provided with --seq, so motif-finding and the motif prior are switched off.");
				}
				//Turn off read distribution clustering
				clusterReadDistributions = (Args.parseFlags(args).contains("noclustering")||!updateBM) ? false : true; 
				
				//Turn off adding franking components
				addFlankingComponent = Args.parseFlags(args).contains("noflanking") ? false : true; 
				
				// Positional prior weights
				posPriorScaling = Args.parseInteger(args,"pospriorscale",posPriorScaling);
				// Turn on per base read filtering
				doReadFilter = Args.parseFlags(args).contains("readfilter") ? true : false;	
				// Markov background threshold to remove sequences with a discovered motif
				MarkovBackSeqRmThres = Args.parseDouble(args, "seqrmthres", MarkovBackSeqRmThres);
				
				//MEME path
				MEMEpath = Args.parseString(args, "memepath", MEMEpath);
				if(!MEMEpath.equals("") && !MEMEpath.endsWith("/")){ MEMEpath= MEMEpath+"/";}
				//MEME args
				MEMEargs = Args.parseString(args, "memeargs", MEMEargs);
				//MEME minw
				MEMEminw = Args.parseInteger(args, "mememinw", MEMEminw);
				//MEME maxw
				MEMEmaxw = Args.parseInteger(args, "mememaxw", MEMEmaxw);
				//MEME nmotifs option
				int MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
				MEMEargs = MEMEargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
				//Enforce non-parallel MEME
				MEMEnonparallel = Args.parseFlags(args).contains("meme1proc");
				
				// Markov background model
				String backfile = Args.parseString(args, "back", null);				
				//Load the background model or make background model
				if (findMotifs){
					try{       
						if(backfile == null)
							markovBackMode = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gen)); // this doesn't seem working for sacCer3
						else
							markovBackMode = BackgroundModelIO.parseMarkovBackgroundModel(backfile, gen);
					} catch (IOException | ParseException e) {
						e.printStackTrace();
					}
				}
				
				//Extra output
				verbose = Args.parseFlags(args).contains("verbose") ? true : false;
				
				//Event reporting mode - determines which events to print to .events file
				if(Args.parseFlags(args).contains("standard") || Args.parseFlags(args).contains("lenient") || Args.parseFlags(args).contains("lenientplus")){
					if(Args.parseFlags(args).contains("standard")){
						standardMode = true; lenientMode=false; lenientPlusMode=false;
					}else if(Args.parseFlags(args).contains("lenient")){
						standardMode = false; lenientMode=true; lenientPlusMode=false;
					}else if(Args.parseFlags(args).contains("lenientplus")){
						standardMode = false; lenientMode=false; lenientPlusMode=true;
					}
				}
				
				//Galaxy html output
				galaxyhtml = Args.parseFlags(args).contains("galaxyhtml") ? true : false;
				
				//Not share subtype motifs across experiments
				shareSubtypes = Args.parseFlags(args).contains("subtypenotshared") ? false : true;
								
				//Initial clustering points
				String fname=null;
				if (ap.hasKey("plist")){
					initialClustPoints = new ArrayList<List<StrandedPoint>>();
					fname=Args.parseString(args, "plist", null);
				}
				if (fname!=null){
					String[] lines= IOUtil.readFile2Array(fname);
					for (int i=0; i <lines.length; i++)
						initialClustPoints.add(RegionFileUtilities.loadStrandedPointsFromFile(gen, lines[i]));
				}
				
				// Stranded read distribution file
				if (ap.hasKey("distrA"))
					distribA=Args.parseString(args, "distrA", null);
				if (ap.hasKey("distrB"))
					distribB=Args.parseString(args, "distrB", null);	
				
				// Motif file list
				if (ap.hasKey("motfile")){
					String motfile=ap.getKeyValue("motfile");
					initMotifs= new ArrayList<WeightMatrix>();
					FreqMatrixImport motifImport = new FreqMatrixImport();
					motifImport.setBackground(markovBackMode);
					initMotifs.addAll(motifImport.readTransfacMatrices(motfile));
				}

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	/**
	 * Merge a set of estimated genomes 
	 * @param estGenomes
	 * @return
	 */
	public Genome mergeGenomes(List<Genome> estGenomes){
		//Combine the chromosome information
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		for(Genome e : estGenomes){
			Map<String, Integer> currMap = e.getChromLengthMap();
			for(String s: currMap.keySet()){
				if(!chrLenMap.containsKey(s) || chrLenMap.get(s)<currMap.get(s))
					chrLenMap.put(s, currMap.get(s));
			}
		}
		gen =new Genome("Genome", chrLenMap);
		return gen;		
	}
	
	//Accessors
	public Genome getGenome(){return gen;}
	public boolean helpWanted(){return printHelp;}
	public double getSigLogConf(){return sigLogConf;}
	public double getPRLogConf(){return prLogConf;}
	public int getMaxThreads(){return maxThreads;}
	public double getAlphaScalingFactor(){return alphaScalingFactor;}
	public double getBetaScalingFactor(){return betaScalingFactor;}
	public double getEpsilonScalingFactor(){return epsilonScalingFactor;}
	public double getFixedAlpha(){return fixedAlpha;}
	public double getWindowExtension(){return extendWindow;}
	public double getMotifMinROC(){return motifMinROC;}
	public int getBMAnalysisWindowMax(){return bmAnalysisWindowMax;}
	public int getAddFlankingComponentSpacing(){return addFlankingComponentSpacing;}
	public boolean getAddFlankingComponents(){return addFlankingComponent;}
	public List<Region> getRegionsToPlot(){return regionsToPlot;}
	public List<Region> getRegionsToIgnore(){return regionsToIgnore;}
	public List<Point> getInitialPos(){return initialPos;}
	public List<StrandedRegion> getMotifRegions(){return motifRegions;}
	public boolean doBMUpdate(){return updateBM;}
	public int getMinComponentsForBMUpdate(){return minComponentsForBMUpdate;}
	public int getMinRefsForBMUpdate(){return minRefsForBMUpdate;}
	public double getMinSubtypeFraction(){return minSubtypeFraction;}
	public double getMinComponentReadFactorForBM(){return minComponentReadFactorForBM;}
	public double getGaussSmoothParam(){return gauss_smooth;}
	public int getMinModelUpdateRounds(){return minModelUpdateRounds;}
	public int getMaxModelUpdateRounds(){return maxModelUpdateRounds;}
	public double getModelConvergenceKL(){return modelConvergenceKL;}
	public boolean getFindingMotifs(){return findMotifs;}
	public boolean useMotifPrior(){return motif_posprior;}
	public boolean getClusteringReads(){return clusterReadDistributions;}
	public boolean useReadFilter(){return doReadFilter;}
	public String getMEMEpath(){return MEMEpath;}
	public String getMEMEargs(){return MEMEargs;}
	public String getMetaMakerArgs(){return MetaMakerArgs;}
	public MarkovBackgroundModel getBackModel(){return markovBackMode;}
	public double getPosPriorScaling(){return posPriorScaling;}
	public int getMinMotifLength(){return MEMEminw;}
	public boolean isVerbose(){return verbose;}
	public List<List<StrandedPoint>> getInitialClustPoints(){return initialClustPoints;}
	public int getInitialCompSpacing(){return initComponentSpacing;}
	public String getDistribA(){return distribA;}
	public String getDistribB(){return distribB;}
	public List<WeightMatrix> getInitialMotifs(){return initMotifs;}
	public double getPreferenceValue(){return preferenceValue;}
	public int getClusterWindowSize(){return clusteringWindow;}
	public int getNumClusteringComps(){return numClusteringComps;}
	public double getMarkovBackSeqRmThres(){return MarkovBackSeqRmThres;}
	public boolean getMEMEnonparallel(){return MEMEnonparallel;}
	public int getModelRange(){return modelRange;}
	public boolean isStandardMode(){return standardMode;}
	public boolean isLenientMode(){return lenientMode;}
	public boolean isLenientPlusMode(){return lenientPlusMode;}
	public boolean useGalaxyhtml(){return galaxyhtml;}
	public boolean getShareSubtypes(){return shareSubtypes;}
	
	
	/**
	 * Make some output directories used by ChExMix
	 */
	public void makeChExMixOutputDirs(boolean makeInterAndImageDirs){
		//Test if output directory already exists. If it does,  recursively delete contents
		outDir =  new File(outName);
//		if(outDir.exists())
//			deleteDirectory(outDir);
		outBase = outDir.getName();
		//(re)make the output directory
		outDir.mkdirs();
		if(makeInterAndImageDirs){
			//Make the gps intermediate results output directory
			interDir = new File(outDir.getAbsolutePath()+File.separator+"intermediate-results");
			interDir.mkdirs();
			//Make the image results output directory
			imagesDir = new File(outDir.getAbsolutePath()+File.separator+"images");
			imagesDir.mkdirs();
		}
	}
	public String getOutName(){return outName;}
	public String getOutBase(){return outBase;}
	public File getOutputParentDir(){return outDir;}
	public File getOutputIntermediateDir(){return interDir;}
	public File getOutputImagesDir(){return imagesDir;}
	
	/**
	 * Delete a directory
	 */
	public boolean deleteDirectory(File path) {
	    if( path.exists() ) {
	      File[] files = path.listFiles();
	      for(int i=0; i<files.length; i++) {
	         if(files[i].isDirectory()) {
	           deleteDirectory(files[i]);
	         }
	         else {
	           files[i].delete();
	         }
	      }
	    }
	    return( path.delete() );
	}
	
	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				" Running ChExMix:\n" +
				"\t--round <max. model update rounds (default=3)>\n" +
				"\t--nomodelupdate [flag to turn off binding model updates]\n" +
				"\t--minmodelupdateevents <minimum number of events to support an update (default=100)>\n" +
				"\t--prlogconf <Poisson log threshold for potential region scanning (default=-6)>\n" +
				"\t--fixedalpha <binding events must have at least this number of reads (default: set automatically)>\n" +
				"\t--alphascale <alpha scaling factor; increase for stricter event calls (default=1.0)>\n" +
				"\t--betascale <beta scaling factor; prior on subtype assignment (default=0.05)>\n" +
				"\t--epsilonscale <epsilon scaling factor; increase for more weight on motif in subtype assignment (default=0.2)>\n" +
				"\t--peakf <file of peaks to initialize component positions>\n" +
				"\t--exclude <file of regions to ignore> OR --excludebed <file of regions to ignore in bed format>\n" +
				"\t--galaxyhtml [flag to produce a html output appropreate for galaxy]\n" +
				" Binding event reporting mode (which events to put into .events file):\n" +
				"\t--standard [report events that pass significance threshold in condition as a whole (default mode)]\n" +
				"\t--lenient [report events that pass significance in >=1 replicate *or* the condition as a whole.]\n" +
				"\t--lenientplus [report events that pass significance in condition OR (>=1 replicate AND no signif diff between replicates)]\n" +
				" Finding ChExMix subtypes using motif:\n"+
				"\t--motfile <file of motifs in transfac format to initialize subtype motifs>\n" +
				"\t--memepath <path to the meme bin dir (default: meme is in $PATH)>\n" +
				"\t--nomotifs [flag to turn off motif-finding & motif priors]\n" +
				"\t--nomotifprior [flag to turn off motif priors only]\n" +
				"\t--memenmotifs <number of motifs MEME should find for each condition (default=3)>\n" +
				"\t--mememinw <minw arg for MEME (default=6)>\n"+
				"\t--mememaxw <maxw arg for MEME (default=18)>\n"+
				"\t--memeargs <additional args for MEME (default=  -dna -mod zoops -revcomp -nostatus)>\n"+
				"\t--minroc <minimum motif ROC value (default=0.7)>\n"+
				"\t--minmodelupdaterefs <minimum number of motif reference to support an subtype distribution update (default=50)>\n"+
				"\t--seqrmthres <Filter out sequences with motifs below this threshold for recursively finding motifs (default=0.1)>\n" +
				" Finding ChExMix subtypes using read distribution clustering:\n"+
				"\t--noclustering [flag to turn off read distribution clustering]\n" +
				"\t--pref <preference value for read distribution clustering (default=-0.1)>\n"+
				"\t--numcomps <number of components to cluster (default=500)>\n"+
				"\t--win <window size of read profiles (default=150)>\n"+
				" Reporting binding events:\n" +
				"\t--q <Q-value minimum (default=0.01)>\n" +
				"\t--minfold <minimum event fold-change vs scaled control (default=1.5)>\n" +
				""));
	}
}

