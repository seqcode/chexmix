package org.seqcode.projects.chexmix.shapealign.alignment;

import java.io.File;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * @author naomi yamada
 */

public class ShapeAlignConfig {
	
	protected GenomeConfig gconf;
	protected SequenceGenerator<Region> seqgen = null;
	protected String[] args;
	protected String MEMEpath="";
	protected String MEMEargs=" -dna -mod zoops -revcomp -nostatus ";    
	public int MEMEminw=6;
	public int MEMEmaxw=18;
	protected int MEMEnmotifs = 3;	
	/** The base name of the output directory */
	protected String outbase;
	/** The output directory file */
	protected File outdir;
	
	protected int win; // window size to get reads
	protected int K; // number of K in K means clustering
	protected List<StrandedPoint> spoints;
	protected List<StrandedRegion> strandedRegions;
	protected double percentile; // percentile to retain only confident alignment
	
	protected boolean weightedPC = false;
	protected boolean isSmithWaterman = false;
	protected boolean printError = false;
	protected boolean printSimScore = false;
	protected boolean printThresDistances = false;
	protected boolean printOffsetArray = false;	
	
	/** Similarity metrics for Smith-Waterman Alignment    **/
	protected boolean euclidean = true; //defualt for Smith-Waterman is euclidean
	protected boolean squaredChi = false;
	protected boolean linear = false;
	protected boolean sorensen = false;
	protected boolean soergel = false;
	protected boolean lorentzian = false;
	protected boolean pce = false;
	protected boolean divergence = false;
	protected boolean clark = false;
	protected boolean AKL = false;
	protected boolean ALLR = false;
	
	// accessors	
	public List<StrandedPoint> getStrandedPoints(){return spoints;}
	public List<StrandedRegion> getStrandedRegions(){return strandedRegions;}
	public String getMemePath(){return MEMEpath;}	
	public String getMEMEargs(){return MEMEargs;}
	public File getOutDir(){return outdir;}
	public String getOutBase(){return outbase;}
	public int getWindowSize(){return win;}
	public int getK(){return K;}	
	public double getPercentile(){return percentile;}
	public boolean isSmithWaterman(){return isSmithWaterman;} // Do Smith-Waterman Alignment
	public boolean weightedPC(){return weightedPC;} // Do weighted person correlation
	public boolean printError(){return printError;}
	public boolean printSimScore(){return printSimScore;}
	public boolean printThresDistances(){return printThresDistances;}
	public boolean printOffsetArray(){return printOffsetArray;}
	public boolean euclidean(){return euclidean;}
	public boolean squaredChi(){return squaredChi;}
	public boolean linear(){return linear;}
	public boolean sorensen(){return sorensen;}
	public boolean soergel(){return soergel;}
	public boolean lorentzian(){return lorentzian;}
	public boolean pce(){return pce;}
	public boolean divergence(){return divergence;}
	public boolean clark(){return clark;}
	public boolean AKL(){return AKL;}
	public boolean ALLR(){return ALLR;}
	
	// setters
	public void setStrandedRegions(List<StrandedRegion> strandedRegs){strandedRegions = strandedRegs;}
	public void setStrandedPoints(List<StrandedPoint> spos){spoints=spos;}
	public void setWindowSize(int w){win=w;}
	
	public ShapeAlignConfig(String[] arguments) {
		
		args = arguments;
		ArgParser ap = new ArgParser(args);
//		if(ap.hasKey("h") || ap.hasKey("help") || args.length == 0){
//			System.err.println(ShapeAlignConfig.getShapAlignArgsList());
//			System.exit(1);
//		}
		
		gconf = new GenomeConfig(args);
		seqgen = new SequenceGenerator<Region>(gconf.getGenome());
		
		win = Args.parseInteger(args, "win", 200);
		K = Args.parseInteger(args, "K", 5);
		percentile = Args.parseInteger(args, "percentile", 10);
		if (ap.hasKey("peakf"))
			spoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peakf"), win);
		
		// Get outdir and outbase and make them; delete dirs that exist with the same
		outbase = Args.parseString(args, "out", "shapealign");
		outdir = new File(outbase);
		
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
		
		// Smith-Waterman similarity metrics
		if (Args.parseFlags(args).contains("linear")){ linear = true;}
		if (Args.parseFlags(args).contains("euclidean")){euclidean = true;}
		if (Args.parseFlags(args).contains("sorensen")){ sorensen = true; ;}
		if (Args.parseFlags(args).contains("soergel")){ soergel = true;}
		if (Args.parseFlags(args).contains("lorentzian")){ lorentzian = true;}
		if (Args.parseFlags(args).contains("pce")){ pce = true;}
		if (Args.parseFlags(args).contains("squaredChi")){ squaredChi = true;}
		if (Args.parseFlags(args).contains("divergence")){ divergence = true;}
		if (Args.parseFlags(args).contains("clark")){ clark = true;}
		if (Args.parseFlags(args).contains("AKL")){ AKL = true;}
		if (Args.parseFlags(args).contains("ALLR")){ ALLR = true;}
		
		if (Args.parseFlags(args).contains("SW")){ isSmithWaterman = true;}
		if (Args.parseFlags(args).contains("weighted")){weightedPC = true;}
		if (Args.parseFlags(args).contains("printError")){printError = true;}
		if (Args.parseFlags(args).contains("printSimScore")){printSimScore = true;}
		if (Args.parseFlags(args).contains("printThresDistances")){printThresDistances = true;}
		if (Args.parseFlags(args).contains("printOffset")){printOffsetArray = true;}		
	}
	
	public void makeOutputDir(){
		//Test if output directory already exists. If it does,  recursively delete contents
		if(outdir.exists())
			deleteDirectory(outdir);
		outbase = outdir.getName();
		//(re)make the output directory
		outdir.mkdirs();
	}
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
	
	public static String getShapAlignArgsList(){
		return(new String("" +
                "ShapeAlignment\n" +
                "--species <organism;genome> OR\n" +
                "--geninfo <genome info> AND --seq <path to seqs>\n" +
                "--peaks <file containing coordinates of peaks> \n" +
                "--win <window of regions to take around peaks> \n" +
                "\nOPTIONS:General\n" +
                "--printSimScore <print similarity scores> \n" +
                "--printThresDistances <print threshould distance matrix> \n" +
                "--printOffset <print offset distributions> \n" +
                "\nOPTIONS:Pearson Correlation\n" +
                "--weighted <weighted PC> \n" +
                "\nOPTIONS:Smith-Waterman\n" +
                "--SW <Smith-Waterman overlap alignment> \n" +
                "Similarity metrics\n" +
                "--euclidean \n" +
                "--squaredChi \n" +
                "--linear \n" +
                "--sorensen \n" +
                "--soergel \n" +
                "--lorentzian \n" +
                "--squaredChi \n" +
                "--pce \n" +
                "--divergence \n" +
                "--clark \n" +
                ""));
	}

}
