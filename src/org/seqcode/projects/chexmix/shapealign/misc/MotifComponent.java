package org.seqcode.projects.chexmix.shapealign.misc;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.io.parsing.PWMParser;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;

/**
 * EM algorithm to find a motif given a prior, meant to be integrated with ChIP-exo shape component. 
 * Reference: Lawrence and Reilly (1990) An Expectation Maximization (EM) Algorithm for the Identification and Characterization of Common Sites in Unaligned Biopolymer Sequences. 
 * @author naomi yamada
 */

public class MotifComponent {
	protected GenomeConfig gconfig;	
	protected List<StrandedPoint> strandedPoints;
	protected int window;
	protected List<String> sequences = new ArrayList<String>();
	
	boolean converged = false;
	final static double epsilon = 0.01;  // parameter for convergence	
	final static int  q = 100 ; // maximum number of EM iterations
	int N; // number of sequences
	int W = 10; // length of the motif    /** hard coded
	int L; // length of each sequences /** hard coded
	double[][][] z; // estimate after q iterations of EM of the probabilities that the site begins at position j in sequence i given the model and the data.  
	double[][][] p; // estimate after q iterations of EM of the probabilities of letter l appearing in position k of the motif
	double[][] po; // estimate after q iterations of EM of the base frequencies of outside of motif
	
	public MotifComponent(GenomeConfig gcon, List<StrandedPoint> spoints, int win, double[][] pwm){
		gconfig = gcon;	
		strandedPoints = spoints;
		window = win;	
		N = strandedPoints.size();
		L = window+1;
		
		z = new double[N][L-W+1][q+1];
		p = new double[4][W][q+1];
		po = new double [4][q+1];		
		
		// initialize all the matrix
		for (int i = 0; i <N ; i++)
			for (int j = 0; j <= L-W; j ++)
				for (int itr = 0 ; itr <= q ; itr++)
					z[i][j][itr] = 0;		
		for (int base = 0 ; base < 4; base++){ // copy PWM to the p matrix
			for (int w = 0; w <W ; w++){
				p[base][w][0] = pwm[base][w];
				for (int itr = 1 ; itr <= q ; itr++)
					p[base][w][itr] = 0;
			}
		}
		for (int base = 0; base <4; base++)
			for (int itr = 0; itr <= q ; itr++)
				po[base][itr] = 0;
	}
	
	public void loadSequencesFromRegions(){
		// converting stranded points to stranded regions
		List<StrandedRegion> regionList = new ArrayList<StrandedRegion>();
		for(Point p: strandedPoints){		
			int start = Math.max(1, p.getLocation() - window/2 );
			int end = Math.min(p.getLocation() + window/2, p.getGenome().getChromLength(p.getChrom()));				
			StrandedRegion strandedReg = new StrandedRegion(p.getGenome(), p.getChrom(), start, end, p.getStrand());					
			regionList.add(strandedReg);
		}
		
		// get sequences from regions
		SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>(gconfig.getGenome());		
		for (StrandedRegion reg : regionList){
			String seq = seqgen.execute(reg);
			sequences.add(seq.toUpperCase());
		}	
	}
	
	// This is created to test the EM with very simple example
	public void loadTestSequence(){
		String seq1 = "TATTTG";
		String seq2 = "TTATTC";
		String seq3 = "GATTTG";
		sequences.add(seq1.toUpperCase());
		sequences.add(seq2.toUpperCase());
		sequences.add(seq3.toUpperCase());
		
		System.out.println("printing test sequences");
		for (String seq : sequences){
			System.out.println(seq);
		}
		
		System.out.println("checking how pwm 2d matrix is arranged");
		for (int base = 0 ; base < p.length; base++ ){
			for (int w = 0; w <W ; w++)
				System.out.print(p[base][w][0]+"\t");
			System.out.println();
		}	
	}
	
	public void runMotifEM(){		
		int round = 0;
        while (!converged && round <q){
        	
        	computeLogLikelihood(round);
        	
        	updatePositions(round);
        	
        	updateMotifFrequencies(round);
        	
        	updateBackgroundBaseFrequencies(round);
        	
        	round ++;
        }	
        	
	}
	
	public double computeLogLikelihood(int round){
		double logLikelihood = 0;
		double motifLikelihood = 0;
		double backgroundLikelihood = 0;
		for (int w = 0 ; w < W ; w++){ // compute likelihood from motif
			for (int b = 0; b <4; b++){
				motifLikelihood += (p[b][w][round]*Math.log(p[b][w][round]));
			}
		}
		for (int b = 0; b <4 ; b++){ // compute likelihood from background
			backgroundLikelihood += (po[b][round]*Math.log(po[b][round]));
		}
		logLikelihood = N*(motifLikelihood + backgroundLikelihood*(L-W));
		
		//printing for test
		System.out.println("log likelihood is "+logLikelihood);
		
		return logLikelihood;		
	}
	
    /**
     * Update positions z based on the current estimate of motif frequencies.
     */
	public void updatePositions(int round){	
		int n = 0;
		for (String seq : sequences){ // for each sequence
			double z_d = 0; // initialize denominator
			double[] z_n = new double[L-W+1];
			for (int j = 0 ; j <= L-W; j ++) //initialize numerator
				z_n[j] = 1;		
			
			for (int j = 0; j <= L-W; j++){
				for(int w = 0; w < W; w++){
					z_n[j] *= p[getBaseIndex(seq,j+w)][w][round];  
				}
				z_d += z_n[j];
			}
			for (int j = 0 ; j <= L-W; j++)
				z[n][j][round] = z_n[j]/z_d;
			n++;
		}
		
		// printing for test
		System.out.println("printing z: iteration # "+round);
		for (int i = 0; i < 10; i++){
			for (int j = 0; j < z[0].length; j++){
				System.out.print(z[i][j][round]+"\t");
			}
			System.out.println();
		}		
	}
	
    /**
     * Update motif frequencies p based on the current estimate of positions z
     */
	public void updateMotifFrequencies(int round){
		double[][] expectedMotifFreq = new double [4][W];
		for (int base = 0; base <4; base++)
			for (int w = 0 ; w <W ; w++)
				expectedMotifFreq[base][w] = 0;
		int n = 0;
		for (String seq : sequences){ // for each sequence
			for (int j = 0; j <= L-W; j++){
				for(int w = 0; w < W; w++){
					expectedMotifFreq[getBaseIndex(seq,j+w)][w] += z[n][j][round]; // make sure that this is correct
				}
			}			
			n++;
		}
		for (int base = 0; base < 4; base++){
			for (int w = 0; w <W ; w++)
				p[base][w][round+1] = expectedMotifFreq[base][w]/sequences.size();
		}
		
		//printing for test
		System.out.println("printing p : iteration # "+round);
		for (int i = 0; i <p.length; i++){
			for (int j = 0; j < p[0].length; j++){
				System.out.print(p[i][j][round+1]+"\t");
			}
			System.out.println();
		}
		
		// check for convergence based on values in p
		converged = true;
		for (int i = 0; i <p.length; i++){
			for (int j = 0; j < p[0].length; j++){
				if (Math.abs(p[i][j][round+1] - p[i][j][round]) > epsilon){
					converged = false;
					break;
				}
			}
		}
	}
	
    /**
     * Update background frequencies based on the current estimate of motif position z , and motif frequencies p
     */
	public void updateBackgroundBaseFrequencies(int round){
		double[] background = new double[4];
		for (int base = 0; base <4; base++)
			background[base] = 0;
		int n = 0;
		for (String seq : sequences){
			double[] totalBaseCounts = new double[4];
			for (int base = 0; base <4 ; base++)
				totalBaseCounts[base] = 0;
			for (int l = 0; l < L ; l++)
				totalBaseCounts[getBaseIndex(seq,l)] ++;  //total base counts for one sequence
			for (int j = 0; j <= L-W ; j ++){ // subtract bases for the window size W each time for j times
				double[] seqBaseCounts = new double[4];
				System.arraycopy(totalBaseCounts, 0, seqBaseCounts, 0, totalBaseCounts.length);
				for (int w = 0; w < W; w++){
					seqBaseCounts[getBaseIndex(seq,w+j)]--;
				}
				for (int base = 0; base <4; base++){
					background[base] += seqBaseCounts[base]*z[n][j][round];
				}
			}
			n++;			
		}
		for (int base = 0; base <4 ; base++)
			po[base][round+1] = background[base]/(sequences.size()*(L-W)); // po total should add up to 1 !!
		
		//printing for test
		System.out.println("printing po : iteration # "+round);
		for (int i = 0; i <po.length; i++){
			System.out.println(po[i][round+1]);
		}
	}	
	
	// get the base index based on the observed sequences
	public int getBaseIndex(String seq, int j){
		int basePos = -1;
		if (seq.charAt(j) == 'A'){
			basePos = 0;
		}else if (seq.charAt(j) == 'C'){
			basePos = 1;
		}else if (seq.charAt(j) == 'G'){
			basePos = 2;
		}else if (seq.charAt(j) == 'T'){
			basePos = 3;
		}else{
			System.err.println("only include ACGT bases");
			System.exit(0);
		}
		return basePos;		
	}	

	// method test
	public static void main(String[] args) throws IOException, ParseException{
		ArgParser ap = new ArgParser(args);		
        if((!ap.hasKey("peaks") && !ap.hasKey("regions")) ) { 
            System.err.println("Usage:\n " +
                               "MotifEM " +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info> AND --seq <path to seqs>\n" +
                               "--peaks <file containing coordinates of peaks> \n" +
                               "--win <window of sequence to take around peaks> OR\n" +
                               "--regions \n" +
                               "");
            System.exit(0);
        }
        
		GenomeConfig gconf = new GenomeConfig(args);
		
		// parsing command line arguments	
		int win = Args.parseInteger(args, "win", 100);
		List<StrandedPoint> strandedPoints = RegionFileUtilities.loadStrandedPointsFromMotifFile(gconf.getGenome(), ap.getKeyValue("peaks"), win);
		
		// This matrix goes with loadTestSequence()
//	    String aLine = " a  0.167 0.500 0.166"; 
//	    String cLine = " c  0.167 0.167 0.167"; 
//	    String gLine = " g  0.166 0.167 0.167";
//	    String tLine = " t  0.500 0.166 0.500"; 
		 
		// PWM for FoxA2 with noise
		// pos			   	1	    2	  3		4	  5		6	  7		8	  9		10
	    String aLine = " a  0.300 0.167 0.300 0.166 0.167 0.200 0.300 0.200 0.300 0.250"; 
	    String cLine = " c  0.200 0.167 0.200 0.167 0.167 0.200 0.200 0.300 0.200 0.250"; 
	    String gLine = " g  0.300 0.166 0.300 0.167 0.166 0.300 0.300 0.200 0.200 0.250";
	    String tLine = " t  0.200 0.500 0.200 0.500 0.500 0.300 0.200 0.300 0.300 0.250"; 
		
	    DenseDoubleMatrix2D pwm  = PWMParser.parsePWM(10, aLine, cLine, gLine, tLine); ///** 3 is hard coded
	    System.out.println(pwm.toString());
	    
	    MotifComponent motifEM = new MotifComponent(gconf, strandedPoints, win, pwm.toArray());
//		motifEM.loadTestSequence();
		motifEM.loadSequencesFromRegions();
		motifEM.runMotifEM();
	}
}
