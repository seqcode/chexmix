package org.seqcode.projects.chexmix.events;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.AnnotationLoader;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;


/**
 * EventsConfig: 
 * 		Maintains settings for events and binding models. 
 *     
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class EventsConfig {
	protected GenomeConfig gconfig;
	protected Genome gen=null;
	protected boolean printHelp=false;
	protected BindingModel defaultModel=null;
	protected boolean addAnnotations=false;
	protected boolean addSequences=true;
	protected List<AnnotationLoader> geneAnnotations = new ArrayList<AnnotationLoader>();
	protected int maxAnnotDistance=50000;
	protected boolean annotOverlapOnly=false;
	protected boolean calcEventBaseCompositions=false; //Calculate base compositions around events and tags belonging to events. Useful for analyzing permanganate ChIP-seq
	protected double qMinThres=0.005;		//Minimum  Q-value for reported binding events
	protected double minEventFoldChange=1.5;
	protected double multiGPSqMinThres=0.01;	// Enrichment test after multiGPS style peak calling to remove non-significant binding events
	protected double multiGPSMinEventFoldChange= 1.2;
	protected double differentialSignificanceP = 0.01;
	protected boolean runDiffTests = true; //Run differential enrichment testing
	protected String Rpath="";
	protected double edger_overdispersion = 0.15; //Overdispersion used by EdgeR differential enrichment tests
	protected boolean eventsFileTXTExtension=false;
	protected boolean printBED=true;
	
	
	//Constants
	public final double LOG2 = Math.log(2);
	public final double LOG_FC_LIMIT = 10; //Maximum absolute log fold-change reported
	public final boolean CALC_EVENTS_LL=false; //Calculate component-wise log-likelihoods during ML

	
	protected String[] args;
	public String getArgs(){
		String a="";
		for(int i=0; i<args.length; i++)
			a = a+" "+args[i];
		return a;
	}
	
	public EventsConfig(GenomeConfig gcon, String[] arguments){
		gconfig = gcon;
		gen = gconfig.getGenome();
		this.args=arguments; 
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
				
				//Read distribution file
				String modelFile = Args.parseString(args, "d", null);	// read distribution file
				if (modelFile != null){
					File pFile = new File(modelFile);
					if(!pFile.isFile()){
						System.err.println("\nCannot find read distribution file: "+modelFile);
						System.exit(1);
					}
					defaultModel = new BindingModel(pFile);
				}
				
				
				/****Gene Annotation****/
				//Gene Annotations
				Collection<String> tfiles = Args.parseStrings(args,"transcripts");
				Collection<String> dbgenes = Args.parseStrings(args,"dbgenes");
				for(String s:dbgenes)
		        	geneAnnotations.add(new AnnotationLoader(gen, s, "refGene", maxAnnotDistance, annotOverlapOnly));
		        for(String s:tfiles)
		        	geneAnnotations.add(new AnnotationLoader(gen, s, "file", maxAnnotDistance, annotOverlapOnly));
				if(geneAnnotations.size()>0)
					addAnnotations=true;
				
				/****Miscellaneous arguments****/
				//Record event base compositions
				calcEventBaseCompositions = Args.parseFlags(args).contains("eventbasecomp") ? true : false;
					
				//Q-value threshold
				qMinThres = Args.parseDouble(args,"q",qMinThres);
				//differential p-value threshold
				differentialSignificanceP = Args.parseDouble(args,"diffp",differentialSignificanceP);
				//Event Fold-change minimum
				minEventFoldChange = Args.parseDouble(args,"minfold",minEventFoldChange);
				
				//Turn off DE testing
				runDiffTests = Args.parseFlags(args).contains("nodifftests") ? false : true;
				//R path
				Rpath = Args.parseString(args, "rpath", Rpath);
				if(Rpath.length()>0 && !Rpath.endsWith(File.separator))
					Rpath = Rpath+File.separator;
				//EdgeR overdispersion parameter
				edger_overdispersion = Args.parseDouble(args,"edgerod",edger_overdispersion);
				
				//Add .txt extension to events files
				eventsFileTXTExtension= Args.parseFlags(args).contains("eventsaretxt");
				//No BED output
				printBED = !(Args.parseFlags(args).contains("nobed"));
				
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
	public boolean isAddingAnnotations(){return addAnnotations;}
	public boolean isAddingSequences(){return addSequences;}
	public List<AnnotationLoader> getGeneAnnotations(){return geneAnnotations;}
	public int getMaxAnnotDistance(){return maxAnnotDistance;}
	public boolean getAnnotOverlapOnly(){return annotOverlapOnly;}
	public BindingModel getDefaultBindingModel(){return defaultModel;}
	public boolean getCalcEventBaseCompositions(){return calcEventBaseCompositions;}
	public double getQMinThres(){return qMinThres;}
	public double getMultiGPSQMinThres(){return multiGPSqMinThres;}
	public double getMinEventFoldChange(){return minEventFoldChange;}
	public double getMultiGPSMinEventFoldChange(){return multiGPSMinEventFoldChange;}
	public double getDiffPMinThres(){return differentialSignificanceP;}
	public boolean getRunDiffTests(){return runDiffTests;}
	public String getRpath(){return Rpath;}
	public double getEdgeROverDisp(){return edger_overdispersion;}
	public boolean getEventsFileTXTExtension(){return eventsFileTXTExtension;}
	public boolean getPrintBED(){return printBED;}
	
	/**
	 * returns a string describing the arguments handled by this parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String("" +
				"BindingModels:\n" +
				"\t--d <read distribution model file>\n" +
				"\t--eventbasecomp [flag to record event base compositions]\n"+
				"EventTesting:\n"+
				"\t--q <Q-value minimum, i.e corrected p-value(default="+qMinThres+")>\n" +
				"\t--minfold <min event fold-change (default="+minEventFoldChange+")>\n" +
				"\t--nodifftests [flag to turn off DE tests]\n" +
				"\t--rpath <path to the R bin dir (default: R is in $PATH). Note that you need to install edgeR separately>\n" +
				"\t--edgerod <EdgeR overdispersion (default="+edger_overdispersion+")>\n" +
				"\t--diffp <minimum p-value for differential enrichment (default="+differentialSignificanceP+")>\n" +
				"\t--eventsaretxt [add .txt to events file extension]\n"+
				"\t--nobed [do not print BED files]" +				
				"Annotations:\n" +
				"\t--transcripts <transcripts file>\n" +
				"\t--dbgenes refGene\n" +
				""));
	}
}
