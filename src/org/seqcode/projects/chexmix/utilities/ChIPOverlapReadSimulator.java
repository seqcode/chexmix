package org.seqcode.projects.chexmix.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.seqcode.deepseq.ReadHit;
import org.seqcode.deepseq.events.BindingModel;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.deepseq.utils.simulation.CountDataSimulator;
import org.seqcode.deepseq.utils.simulation.CountDataSimulator.SimCounts;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.Pair;


/**
 * Simulates single or double condition fragments and reads using BindingModels. <br> 
 * 
 * In any given simulated replicate experiment, reads are simulated according to a BindingModel.
 * Read counts for each binding event in each condition are determined from a CountDataSimulator,
 * which allows for differential binding between two conditions, and negative-binomial sampled
 * read counts in each replicate. 
 * Noise reads are distributed uniformly across the genome or sampled from an input tag file.
 * Events can have less than 200bp in distance to model closely spaced events.
 * 
 * 
 * @author Naomi Yamada
 *
 * Extended from org.seqcode.deepseq.utils.simulation.ChIPReadSimulator
 */
public class ChIPOverlapReadSimulator {

	private BindingModel model;
	private int numConditions = 2;
	private int numReplicates = 1;
	private String outPath;
	private FileWriter[][] writers;
	
	private Genome fakeGen;
	private long[] chromLens;
	private HashMap<String, Long> chromOffsets = new HashMap<String, Long>();
	private int rLen=32;
	private long genomeLength=-1;
	private double noiseProbabilities[][];
	private int numSigFrags[][], numTotalFrags[][];
	private int numReads =1000000;
	private List<SimCounts> simCounts;
	
	private int numEvents=0;
	private int eventSpacing = 2000;
	private double jointEventRate = 0.0;
	private int jointEventSpacing = 200;
	private int fragLenMean = 140;
	private int fragLenStdDev = 40;
	private List<Pair<Point, SimCounts>> events = new ArrayList<Pair<Point, SimCounts>>();
	private HashMap<Point, Boolean> eventIsJoint = new HashMap<Point, Boolean>();
	private List<ReadHit> noiseSource=null;
	private boolean subsampleControl=false;
	private boolean paired=false;
	private int overlapEvent=1000;
	
	public ChIPOverlapReadSimulator(BindingModel m, Genome g, List<SimCounts> counts, int numCond, int numRep, double noiseProb, double jointRate, int jointSpacing, int overlapEvent, String outPath){
		model=m;
		numConditions = numCond;
		numReplicates = numRep;
		jointEventRate = jointRate;
		jointEventSpacing = jointSpacing;
		this.outPath = outPath;
		fakeGen = g;
		genomeLength = (long)fakeGen.getGenomeLength();
		chromLens = new long[fakeGen.getChromList().size()];
		this.overlapEvent=overlapEvent;
		
		int c=0; long offset=0;
		for(String chr : fakeGen.getChromList()){
			chromLens[c]=fakeGen.getChromLength(chr); 
			chromOffsets.put(chr, offset);
			offset+=chromLens[c];
			c++;
		}
		
		noiseProbabilities = new double[numConditions][numReplicates];
		numSigFrags = new int[numConditions][numReplicates];
		numTotalFrags = new int[numConditions][numReplicates];
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
				noiseProbabilities[co][r]=noiseProb;
				numSigFrags[co][r]=0;
				numTotalFrags[co][r]=0;
			}
		
		try {
			writers = new FileWriter[numConditions][numReplicates];
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					FileWriter fout = new FileWriter(outPath+"_reads_C"+co+"_R"+r+".bed");
					writers[co][r]=fout;
				}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if(noiseProb<1.0){
			simCounts = counts;
			setBindingPositions();
			numEvents = events.size();
			
			//Count signal reads and infer total reads
			for(Pair<Point,SimCounts> event : events){
				int index=0;
				for(int co=0; co<numConditions; co++)
					for(int r=0; r<numReplicates; r++){
						if(!eventIsJoint.get(event.car())){
							numSigFrags[co][r]+=event.cdr().counts[index];
						}else{
							numSigFrags[co][r]+=event.cdr().backup[index];
						}
						numTotalFrags[co][r] = (int)((double)numSigFrags[co][r]/(1-noiseProbabilities[co][r]));
						index++;
					}
			}
		}
	}
	
	/**
	 * Make binding positions corresponding to all simulated counts.
	 * By default, shared events are in chr1, and diff events are in chrX
	 * If there is no chrX (e.g. yeast), shared and diff events are interspersed along all chroms
	 */
	private void setBindingPositions(){
		Random jointDice = new Random();
		long sharedOffset=eventSpacing, diffOffset=eventSpacing/2;
		if(chromOffsets.containsKey("X")){
			sharedOffset=chromOffsets.get("1"); diffOffset=chromOffsets.get("X");
		}		
		
		int eventC=0;
		for(SimCounts s : simCounts){
			long curroff =0;
			if(s.isDiff){
				curroff =diffOffset;
				diffOffset+=eventSpacing;
			}else{				
				if (overlapEvent>0 && eventC < overlapEvent && (eventC%2)==1){
					curroff = sharedOffset+jointDice.nextInt(jointEventSpacing-1)+1;					
				}else{				
					curroff =sharedOffset;
					sharedOffset+=eventSpacing;
				}
			}
			
			//Translate from offsets to chromosome name and start
			int c=0; long offset=0;
			String chr = fakeGen.getChromList().get(0);
			while(curroff>(offset+chromLens[c]) && c<fakeGen.getChromList().size()-1){
				c++;
				chr = fakeGen.getChromList().get(c);
				offset = chromOffsets.get(chr);
			}
			long start = curroff-offset;
			Point p = new Point(fakeGen, chr, (int)start);
			events.add(new Pair<Point, SimCounts>(p,s));
			eventIsJoint.put(p, false);
			
			//Simulate a joint event?
			double jointrand = jointDice.nextDouble();
			if(jointrand<jointEventRate){
				Point jp = new Point(fakeGen, chr, (int)start+jointEventSpacing);
				events.add(new Pair<Point, SimCounts>(jp,s));
				eventIsJoint.put(jp, true);
			}
			eventC++;
		}
	}
	
	
	/**
	 * Simulate a set of binding event reads for each replicate in each condition. 
	 * Number of reads, number of events, and strengths of events are all pre-determined. 
	 */
	private void simulateReads(){
		//Initialize the probability landscape
		int eventWidth=1000; int evoff = eventWidth/2;
		double[] forProbLand=new double[eventWidth]; double[] revProbLand=new double[eventWidth];
		double[] forProbCumul=new double[eventWidth]; double[] revProbCumul=new double[eventWidth];
		for(int i=0; i<eventWidth; i++){
			forProbLand[i]=0; revProbLand[i]=0;
			forProbCumul[i]=0; revProbCumul[i]=0;
		}
		int modelRange = Math.max(Math.abs(model.getMin()), Math.abs(model.getMax()));
		int winStart = Math.max(0, evoff-modelRange);
		int winStop = Math.min(evoff+modelRange, eventWidth-1);
		for(int i=winStart; i<winStop; i++){
			int forDist  = i-evoff;
			int revDist  = evoff-i;
			forProbLand[i]+=model.probability(forDist);
			revProbLand[i]+=model.probability(revDist);
		}
		//Set the cumulative scores
		double fTotal=0, rTotal=0;
		for(int i=0; i<eventWidth; i++){
			fTotal+=forProbLand[i];
			rTotal+=revProbLand[i];
			forProbCumul[i]=fTotal;
			revProbCumul[i]=rTotal;
		}
		//Normalize
		for(int i=0; i<eventWidth; i++){
			forProbLand[i]=forProbLand[i]/fTotal;
			revProbLand[i]=revProbLand[i]/rTotal;
			forProbCumul[i]=forProbCumul[i]/fTotal;
			revProbCumul[i]=revProbCumul[i]/rTotal;
		}
		
		
		//Generate the reads
		Random sigGenerator = new Random();
		Random noiseGenerator = new Random();
		Random strandGenerator = new Random();
		NormalDistribution fragLengthDistrib = new NormalDistribution(fragLenMean, fragLenStdDev);
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
				List<ReadHit> frags = new ArrayList<ReadHit>();
				List<ReadHit> fragPairs = new ArrayList<ReadHit>();
				
				// Generate event reads
				if(noiseProbabilities[co][r]<1){
					for(Pair<Point, SimCounts> ps : events){
						Point coord = ps.car();
						String chr = coord.getChrom();
						int chrLen = fakeGen.getChromLength(chr);
						SimCounts sc = ps.cdr();
						int sample = co*numReplicates+r;
						double readCount = eventIsJoint.get(coord) ? sc.backup[sample] : sc.counts[sample];
						
						ReadHit rh = null; ReadHit rhp=null;
						for(int x=0; x<readCount; x++){
							boolean forwardStrand = strandGenerator.nextDouble()<0.5 ? true:false;
							double rand = sigGenerator.nextDouble();
							int fivePrimeEnd=coord.getLocation()-evoff;
							//Find the probability interval
							if(forwardStrand){
								for(int j=0; j<eventWidth; j++){
									if(forProbCumul[j] > rand){
										fivePrimeEnd+=j;
										break;
									}
								}
								//Make the ReadHit
								if((fivePrimeEnd+rLen-1)<chrLen){
									rh = new ReadHit(chr, fivePrimeEnd, fivePrimeEnd+rLen-1, '+');
									if(paired){
										int len =Math.max(rLen, (int)fragLengthDistrib.sample()); 
										rhp = new ReadHit(chr, fivePrimeEnd+len-rLen+1, fivePrimeEnd+len, '-');
									}
								}
							}else{
								for(int j=eventWidth-1; j>=0; j--){
									if(revProbCumul[j] < rand){
										fivePrimeEnd+=j+1;
										break;
									}
								}
								//Make the ReadHit
								if(fivePrimeEnd<chrLen){
									rh = new ReadHit(chr, Math.max(1, fivePrimeEnd-rLen+1), fivePrimeEnd, '-');
									if(paired){
										int len =Math.max(rLen, (int)fragLengthDistrib.sample()); 
										rhp = new ReadHit(chr, Math.max(1, fivePrimeEnd-len), Math.max(1, fivePrimeEnd-len+rLen-1), '+');
									}
								}
							}
							if(rh!=null && (!paired || rhp!=null)){
								frags.add(rh);
								fragPairs.add(rhp);
							}
						}
					}
				}
				
				//Generate noise reads
				Random readSampler = new Random();
				int noiseFrags = (int)((double)numTotalFrags[co][r]*noiseProbabilities[co][r]);
				if(noiseSource==null) //Poisson
					for(int i=0; i<noiseFrags; i++){
						ReadHit rh=null, rhp=null;
						double noiserand = noiseGenerator.nextDouble();
						double strandrand = strandGenerator.nextDouble();
	
						long pos = (long)(noiserand*(genomeLength));
						
						//Translate from pos to chromosome name and start
						int c=0; long offset=0;
						String chr = fakeGen.getChromList().get(0);
						while(pos>(offset+chromLens[c]) && c<fakeGen.getChromList().size()-1){
							c++;
							chr = fakeGen.getChromList().get(c);
							offset = chromOffsets.get(chr);
						}
						long start = pos-offset;
						
						//Add the ReadHit
						if (strandrand<0.5){
							rh = new ReadHit(chr, (int)start, (int)start+rLen-1, '+');
							if(paired){
								int len =Math.max(rLen, (int)fragLengthDistrib.sample()); 
								rhp = new ReadHit(chr, (int)start+len-rLen+1, (int)start+len, '-');
							}
						}else{
							rh = new ReadHit(chr, Math.max(1, (int)start-rLen+1), (int)start, '-');
							if(paired){
								int len =Math.max(rLen, (int)fragLengthDistrib.sample()); 
								rhp = new ReadHit(chr, Math.max(1, (int)start-len), Math.max(1, (int)start-len+rLen-1), '+');
							}
						}
						frags.add(rh);
						if(paired)
							fragPairs.add(rhp);
					}
				else{
					//If the noise source is a control experiment, we shouldn't sample fragments with replacement.
					//Rather, we treat each control experiment read as a distinct fragment in the initial library 
					int noiseSourceSize = noiseSource.size();
					if(noiseSourceSize<noiseFrags){
						System.err.println("Provided control has fewer reads than the requested noise fragments");
						System.exit(1);
					}
					//fragindex is a list of non-repeating random numbers
					Integer[] fragindex = new Integer[noiseSourceSize];
				    for (int i = 0; i < noiseSourceSize; i++)
				    	fragindex[i] = i;
				    Collections.shuffle(Arrays.asList(fragindex));
					//Extract the fragments
					for(int i=0; i<noiseFrags; i++){
						ReadHit rh = noiseSource.get(fragindex[i]);
						frags.add(rh);
					}
				}
				
				//Count unique fragments
				HashMap<String, Integer> uniqueFragLabels = new HashMap<String, Integer>();
				for(ReadHit hit : frags){
					String label = hit.toString();
					int count =0;
					if(uniqueFragLabels.containsKey(label))
						count = uniqueFragLabels.get(label);
					uniqueFragLabels.put(label, count+1);
				}
				int uniqueFrags = uniqueFragLabels.keySet().size();
				System.out.println(uniqueFrags+" unique positions generated.");
				
				//Sample reads from the fragments & print
				try {
					int numFrags = frags.size();
					for(int x=0; x<numReads; x++){
						ReadHit rh=null, rhp=null;
						if(noiseSource==null || !subsampleControl){
							double rand = readSampler.nextDouble();
							int index = (int)((double)numFrags*rand);
							rh = frags.get(index);
							if(paired)
								rhp=fragPairs.get(index);
						}else{
							//Here, the options tell us to just *subsample* the control reads directly without replacement. 
							//This is easy here, as frags should contain a randomly ordered subset of control reads already
							rh = frags.get(x);
						}
						writers[co][r].write(rh.getChrom()+"\t"+rh.getStart()+"\t"+rh.getEnd()+"\tU\t0\t"+rh.getStrand()+"\n");
						if(paired)
							writers[co][r].write(rhp.getChrom()+"\t"+rhp.getStart()+"\t"+rhp.getEnd()+"\tU\t0\t"+rhp.getStrand()+"\n");
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		
		
	}
		
	//Accessors
	public void setTotalFrags(int totFrags){
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
					numTotalFrags[co][r] = totFrags;
			}
	}
	public void setTotalReads(int totReads){
		for(int co=0; co<numConditions; co++)
			for(int r=0; r<numReplicates; r++){
					numReads = totReads;
			}
	}
	public void setJointEventRate(double jointRate){
		this.jointEventRate = jointRate;
	}
	public void setReadLength(int r){
		this.rLen=r;
	}
	public void setPaired(boolean p){ this.paired=p;}
	public void setNoiseSource(List<Sample> controls, boolean subsample){
		subsampleControl = subsample;
		noiseSource = new ArrayList<ReadHit>();
		for(Sample s : controls){
			//Add each read hit as its own fragment
			for(ReadHit rh : s.exportReadHits(rLen)){
				for(float x=0; x<rh.getWeight(); x+=1.0){
					ReadHit newRH = new ReadHit(rh.getChrom(), rh.getStart(), rh.getEnd(), rh.getStrand());
					noiseSource.add(newRH);
				}
			}
		}
		System.err.println(noiseSource.size()+" control reads sourced as distinct fragments");
	}
	

	// clean up
	public void close(){
		try {
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++){
					writers[co][r].close();
				}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Write an output file of event positions & expected read counts per replicate
	 */
	public void printEvents(){
		try {
			String filename = outPath+".events";
			FileWriter fout = new FileWriter(filename);
			
			fout.write("Coord\tDiff");
			for(int co=0; co<numConditions; co++)
				for(int r=0; r<numReplicates; r++)
					fout.write("\tC"+co+"_R"+r);
			fout.write("\n");
			
			for(Pair<Point, SimCounts> ps : events){
				Point coord = ps.car();
				SimCounts sc = ps.cdr();
				int diff = sc.isDiff && !eventIsJoint.get(coord) ? 1 : 0;
				fout.write(coord.getLocationString()+"\t"+diff);
				for(int x=0; x<sc.counts.length; x++){
					if(!eventIsJoint.get(coord))
						fout.write("\t"+sc.counts[x]);
					else
						fout.write("\t"+sc.backup[x]);
				}fout.write("\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**

	 */
	public static void main(String[] args) {
		String empFile, outFile="out";
		int c=2, r=2, numdata, jointSpacing=200, rlen=32;
		int overlapEvent=1000;
		double frags=1000000, reads=1000000, a, up, down, diff, jointRate=0.0;
		String bmfile;
		boolean printEvents=true, subsampleControl=false, isPaired=false;
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h") || !ap.hasKey("emp")){
			System.err.println("ChIPReadSimulator:\n" +
					"\t--geninfo <genome info file>\n" +
					"\t--emp <empirical data file>\n" +
					"\t--numdata <number of events to simulate>\n" +
					"\t--c <num conditions>\n" +
					"\t--r <num replicates per condition>\n" +
					"\t--a <over-dispersion param>\n" +
					"\t--frags <avg total fragments>\n" +
					"\t--reads <avg total reads>\n" +
					"\t--rlen <length of each read>\n" +
					"\t--up <up-regulated fraction>\n" +
					"\t--down <down-regulated fraction>\n" +
					"\t--diff <basis of differential expression>\n" +
					"\t--model <binding model file>\n" +
					"\t--noise <noise probability per replicate>\n" +
					"\t--ctrl <control experiment to replace Poisson for noise data>\n" +
					"\t--subsample <subsample the control experiment non-redundantly>\n" +
					"\t--format <SAM/IDX>\n" +
					"\t--jointrate <proportion of peaks that are joint binding events>\n" +
					"\t--jointspacing <spacing between joint events>\n" +
					"\t--noevents [flag to turn off making some files]\n" +
					"\t--paired [flag to generate paired reads (not working for control-sampling yet)]\n" +
					"\t--out <output file>\n" +
					"");
		}else{
			GenomeConfig gcon = new GenomeConfig(args);
			ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
			CountDataSimulator cdsim = new CountDataSimulator();
			
			//////////////////////////////////////////////////
			// Read in parameters 
			if(ap.hasKey("emp")){
				empFile = ap.getKeyValue("emp");
				cdsim.loadEmpiricalFromFile(empFile);
			}if(ap.hasKey("numdata")){
				numdata = new Integer(ap.getKeyValue("numdata"));
				cdsim.setDataPoints(numdata);
			}if(ap.hasKey("c")){
				c = new Integer(ap.getKeyValue("c"));
				if(c!=1 && c!=2){
					System.err.println("Number of conditions must be either 1 or 2");
					System.exit(1);
				}
				cdsim.setConditions(c);
			}if(ap.hasKey("r")){
				r = new Integer(ap.getKeyValue("r"));
				cdsim.setReplicates(r);
			}if(ap.hasKey("a")){
				a = new Double(ap.getKeyValue("a"));
				cdsim.setAlpha(a);
			}if(ap.hasKey("frags")){
				frags = new Double(ap.getKeyValue("frags"));
				cdsim.setReadsA(frags);
				cdsim.setReadsB(frags);
			}if(ap.hasKey("reads")){
				reads = new Double(ap.getKeyValue("reads"));
			}if(ap.hasKey("rlen")){
				rlen = new Integer(ap.getKeyValue("rlen"));
			}if(ap.hasKey("up")){
				up = new Double(ap.getKeyValue("up"));
				cdsim.setUpRegFrac(up);
			}if(ap.hasKey("down")){
				down = new Double(ap.getKeyValue("down"));
				cdsim.setDownRegFrac(down);
			}if(ap.hasKey("diff")){
				diff = new Double(ap.getKeyValue("diff"));
				cdsim.setDiffExpLevel(diff);
			}if(ap.hasKey("jointrate")){
				jointRate = new Double(ap.getKeyValue("jointrate"));
			}if(ap.hasKey("jointspacing")){
				jointSpacing = new Integer(ap.getKeyValue("jointspacing"));
			}if(ap.hasKey("out")){
				outFile = ap.getKeyValue("out");
			}if(ap.hasKey("noevents")){
				printEvents=false;
			}if(ap.hasKey("subsample")){
				subsampleControl=true;
				if(reads>frags){
					System.err.println("Can't subsample in cases where frags is less than reads!");
					System.exit(1);
				}
			}if(ap.hasKey("paired")){
				isPaired=true;
			}if(ap.hasKey("overlap")){
				overlapEvent = new Integer(ap.getKeyValue("overlap"));
			}
			
			
			double noiseProb   = Args.parseDouble(args, "noise", 0.9);
			double [][] noiseProbs = new double[c][r];
			for(int x=0; x<c; x++)
				for(int y=0; y<noiseProbs[x].length; y++)
					noiseProbs[x][y]=noiseProb;
			cdsim.setReads(frags*(1-noiseProb));
			
			bmfile  = Args.parseString(args, "model", null);
			File mFile = new File(bmfile);
			if(!mFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
			
			//////////////////////////////////////////////////
			// Simulate counts 
			List<SimCounts> counts = null;
			if(noiseProb<1.0){
				counts = cdsim.simulate();
				if(printEvents)
					cdsim.printOutput(outFile, true);
			}
			
			//////////////////////////////////////////////////
			// Simulate reads according to counts and binding model
			BindingModel bm = new BindingModel(mFile);
	        //Initialize the MultiConditionReadSimulator
			ChIPOverlapReadSimulator sim = new ChIPOverlapReadSimulator(bm, gcon.getGenome(), counts, c, r, noiseProb, jointRate, jointSpacing, overlapEvent, outFile);
	        if(noiseProb==1.0)
	        	sim.setTotalFrags((int) frags);

	        ExperimentManager manager=null;
	        if(ap.hasKey("ctrl")){
				manager = new ExperimentManager(econ);
				sim.setNoiseSource(manager.getSamples(), subsampleControl);
			}
			
	        sim.setTotalReads((int) reads);
	        sim.setReadLength(rlen);
	        sim.setPaired(isPaired);
	        
	        if(noiseProb<1 && printEvents)
	        	sim.printEvents();	      
			sim.simulateReads();
			sim.close();
			
			if(manager!=null)
				manager.close();
		}
	}
}
