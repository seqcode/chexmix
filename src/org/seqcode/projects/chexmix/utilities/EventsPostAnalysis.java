package org.seqcode.projects.chexmix.utilities;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.RealValuedHistogram;
import org.seqcode.motifs.DrawMotifs;
import org.seqcode.projects.chexmix.events.BindingEvent;
import org.seqcode.projects.chexmix.events.BindingManager;
import org.seqcode.projects.chexmix.events.BindingSubtype;
import org.seqcode.projects.chexmix.events.EventsConfig;
import org.seqcode.projects.chexmix.framework.ChExMixConfig;
import org.seqcode.projects.chexmix.motifs.MotifPlatform;
import org.seqcode.viz.metaprofile.MetaConfig;
import org.seqcode.viz.metaprofile.SequenceAlignmentFigure;


public class EventsPostAnalysis {
	
	protected GenomeConfig gconfig;
	protected EventsConfig evconfig;
	protected ChExMixConfig config;
	protected ExperimentManager manager;
	protected BindingManager bindingManager;
	protected MotifPlatform motifFinder = null;
	protected List<BindingEvent> events;
	protected double motifThres = 0.6;  // fraction of max score threshold
	
	public EventsPostAnalysis(GenomeConfig gcon,EventsConfig ec, ChExMixConfig c, ExperimentManager man, BindingManager bMan, List<BindingEvent> ev, MotifPlatform mp){
		evconfig = ec;
		config = c;
		manager = man;
		bindingManager = bMan;
		events = ev;
		motifFinder = mp;
		gconfig= gcon;
	}
	
	/**
	 * Run post-analysis of peaks.
	 * 	1) Histograms of peak-closestMotif distances
	 * 	2) Histograms of peak-peak distances (same condition)
	 * 	3) Histograms of peak-peak distances (inter-condition) 
	 */
	public void execute(int histoWin){
		System.err.println("Events post-analysis");
		String pcmfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".peaks2motifs.histo.txt";
		String ppdscfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".intraCondPeakDistances.histo.txt";
		String ppdicfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".interCondPeakDistances.histo.txt";
		String htmlfilename = config.getOutputParentDir()+File.separator+"ChExMix_"+config.getOutBase()+"_results.html";
		
		//0) Set up hash map structure for events by chromosome
		List<HashMap<String,List<Integer>>> eventStruct = new ArrayList<HashMap<String,List<Integer>>>();
		for(int c=0; c<manager.getNumConditions(); c++){
			ExperimentCondition cond = manager.getIndexedCondition(c);
			eventStruct.add(new HashMap<String, List<Integer>>());
			for(String chr : config.getGenome().getChromList())
				eventStruct.get(c).put(chr, new ArrayList<Integer>());
			for(BindingEvent ev : events){
				double Q = ev.getCondSigVCtrlP(cond);
	    		if(ev.isFoundInCondition(cond) && Q <=evconfig.getQMinThres()){
					String chr = ev.getPoint().getChrom();
					int loc = ev.getPoint().getLocation();
					eventStruct.get(c).get(chr).add(loc);
				}
			}
		}
		
		
		//1) Histograms of peak-closestMotif distances
		try {
			if(config.getFindingMotifs()){
				System.err.println("\tPeak-motif distance histograms");	    		
	    		FileWriter fout = new FileWriter(pcmfilename);
	    		fout.write("#Peaks to closest motifs distance histograms\n\n");
				for(ExperimentCondition cond : manager.getConditions()){
					List<BindingSubtype> subtypes = bindingManager.getBindingSubtype(cond);
					fout.write("#Condition:"+cond.getName()+"\n");
					RealValuedHistogram peakMotifHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
					for(BindingEvent ev : events){
						double Q = ev.getCondSigVCtrlP(cond);
						if(ev.isFoundInCondition(cond) && Q <=evconfig.getQMinThres()){			
							int loc = ev.getPoint().getLocation();
							if(ev.getContainingRegion()!=null){	
								if(loc - ev.getContainingRegion().getStart() > histoWin && ev.getContainingRegion().getEnd()-loc >histoWin){
									double[][] scores = motifFinder.scanRegionWithMotif(ev.getContainingRegion(), cond);
									int index = loc - ev.getContainingRegion().getStart();
									int closestMatch = Integer.MAX_VALUE;
									for (int m=0; m < subtypes.size(); m++){
										if (subtypes.get(m).hasMotif()){
											double currThreshold = subtypes.get(m).getMotif().getMaxScore() * motifThres;
											for(int x=0; x<scores[m].length; x++)
												if(scores[m][x]>=currThreshold && Math.abs(x-index)<closestMatch)
													closestMatch = Math.abs(x-index);
										}
									}
									peakMotifHisto.addValue(closestMatch);										
								}
							}
						}
					}
					fout.write(peakMotifHisto.contentsToString()+"\n");
				}
				fout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//2) Histograms of peak-peak distances (same condition)
		try {
			System.err.println("\tPeak-peak distance histograms (same condition)");    		
    		FileWriter fout = new FileWriter(ppdscfilename);
    		fout.write("#Peaks to other peaks in same condition distance histograms\n\n");
			for(ExperimentCondition cond : manager.getConditions()){
				RealValuedHistogram peakPeakHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
				fout.write("#Condition: "+cond.getName()+"\n");
				for(String chr : config.getGenome().getChromList()){
					List<Integer> currCondChrLocs = eventStruct.get(cond.getIndex()).get(chr);
					for(int x=0; x<currCondChrLocs.size(); x++){
						int xLoc = currCondChrLocs.get(x);
						int closestPeak = Integer.MAX_VALUE;
						for(int y=0; y<currCondChrLocs.size(); y++){ if(x!=y){
							int yLoc = currCondChrLocs.get(y);
							int dist = Math.abs(xLoc-yLoc);
							if(dist<closestPeak)
								closestPeak = dist;
						}}
						peakPeakHisto.addValue(closestPeak);	
					}
				}
				fout.write(peakPeakHisto.contentsToString()+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//3) Histograms of peak-peak distances (inter-condition)
		try {
			if(manager.getNumConditions()>1){
				System.err.println("\tPeak-peak distance histograms (different conditions)");	    		
	    		FileWriter fout = new FileWriter(ppdicfilename);
	    		fout.write("#Peaks to peaks in other conditions distance histograms\n\n");
				for(ExperimentCondition condA : manager.getConditions()){
					for(ExperimentCondition condB : manager.getConditions()){if(condA != condB){
						RealValuedHistogram peakPeakHisto = new RealValuedHistogram(0, histoWin, histoWin/5);
						fout.write("#Condition: "+condA.getName()+" vs "+condB.getName()+"\n");
						for(String chr : config.getGenome().getChromList()){
							List<Integer> currCondChrLocsA = eventStruct.get(condA.getIndex()).get(chr);
							List<Integer> currCondChrLocsB = eventStruct.get(condB.getIndex()).get(chr);
							for(int x=0; x<currCondChrLocsA.size(); x++){
								int xLoc = currCondChrLocsA.get(x);
								int closestPeak = Integer.MAX_VALUE;
								for(int y=0; y<currCondChrLocsB.size(); y++){ if(x!=y){
									int yLoc = currCondChrLocsB.get(y);
									int dist = Math.abs(xLoc-yLoc);
									if(dist<closestPeak)
										closestPeak = dist;
								}}
								peakPeakHisto.addValue(closestPeak);	
							}
						}
						fout.write(peakPeakHisto.contentsToString()+"\n");
					}}
				}
				fout.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//4) Print motif transfac file
		try {
			if(config.getFindingMotifs()){
				for(ExperimentCondition cond : manager.getConditions()){
					String motfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"."+cond.getName()+".transfac";
					FileWriter fout = new FileWriter(motfilename);
					String out = "";
					for (BindingSubtype sub : bindingManager.getBindingSubtype(cond))
						if (sub.hasMotif())
							out = out+WeightMatrix.printTransfacMatrix(sub.getFreqMatrix(),sub.getFreqMatrix().getName());
					fout.write(out);
					fout.close();
				}
			}
		}catch (IOException e) {
			e.printStackTrace();
		}
				
		//5) Print aligned sequences and 4 color sequence plot  
		if (gconfig.getSequenceGenerator().usingLocalFiles()){
			for(ExperimentCondition cond : manager.getConditions()){		
				try {
					String outFile = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_seq.png";
					String seqOutFile = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"."+cond.getName()+".seq";
					SequenceGenerator seqgen = gconfig.getSequenceGenerator();
		    	
					List<StrandedRegion> regions = new ArrayList<StrandedRegion>();
					for (List<StrandedPoint> pl : bindingManager.getAlignedEventPoints(cond))            
						for (StrandedPoint p : pl)
							regions.add(p.expand(evconfig.SEQPLOTWIN, evconfig.SEQPLOTWIN));
		                
					List<String> seqs = RegionFileUtilities.getSequencesForStrandedRegions(regions, seqgen);
		        	
					if(seqs !=null){
						SequenceAlignmentFigure fig = new SequenceAlignmentFigure();
						fig.setColors(Color.RED, Color.BLUE, Color.ORANGE, Color.GREEN);
						fig.visualizeSequences(seqs, 3, 1, new File(outFile));
				    	
						if(seqOutFile != null){
							FileWriter fout = new FileWriter(seqOutFile);
							for(String s: seqs){
								fout.write(String.format("%s\n", s));
							}
							fout.close();
						}
					}
				}catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
				
		//6) Make heatmap
		for (ExperimentCondition cond : manager.getConditions()){
			String pointArgs = " --peaks "+config.getOutputParentDir()+File.separator+config.getOutBase()+"_"+cond.getName()+".subtype.aligned.events";
			// Run for each strand
			System.out.println(config.getMetaMakerArgs()+pointArgs+" --strand + --color blue --noborder --out "+cond.getName()+".events");
			runMetaMaker(config.getMetaMakerArgs()+pointArgs+" --strand + --color blue --noborder --out "+cond.getName()+".events");
			runMetaMaker(config.getMetaMakerArgs()+pointArgs+" --strand - --color red --noborder --out "+cond.getName()+".events");
					
			// Combine plots
			for (ExperimentCondition pcond : manager.getConditions()){				
				String pngPath=config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+pcond.getName()+".events_"+cond.getName()+"_";
				try {
					Process proc = Runtime.getRuntime().exec("composite -dissolve 60,100 -transparent-color white "+pngPath+"+_lines.png "+pngPath+"-_lines.png "+pngPath+"heatmap.png");
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}				
		
		//7) HTML report
		try {
			System.err.println("Writing results report to: "+htmlfilename);
			
			//Write motif images
			HashMap<ExperimentCondition, List<String>> motifImageNames = new HashMap<ExperimentCondition, List<String>>();
			HashMap<ExperimentCondition, List<String>> motifRCImageNames = new HashMap<ExperimentCondition, List<String>>();
			if(config.getFindingMotifs()){
				for(ExperimentCondition cond : manager.getConditions()){
					List<String> mImageNames = new ArrayList<String>();
					List<String> mRCImageNames = new ArrayList<String>();
					for (BindingSubtype subtype : bindingManager.getBindingSubtype(cond)){
						if (subtype.hasMotif()){
							WeightMatrix currMotif = subtype.getMotif();
							String imName = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_"+WeightMatrix.getConsensus(currMotif)+".png";
							String imName2 = "images/"+config.getOutBase()+"_"+cond.getName()+"_"+WeightMatrix.getConsensus(currMotif)+".png";
							String motifLabel = cond.getName()+" "+WeightMatrix.getConsensus(currMotif)+", MEME";
							DrawMotifs.printMotifLogo(currMotif, new File(imName), 75, motifLabel);
							mImageNames.add(imName2);
							WeightMatrix wm_rc = WeightMatrix.reverseComplement(currMotif);
							imName = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_"+WeightMatrix.getConsensus(currMotif)+"_rc.png";
							imName2 = "images/"+config.getOutBase()+"_"+cond.getName()+"_"+WeightMatrix.getConsensus(currMotif)+"_rc.png";
							motifLabel = cond.getName()+" revcomp "+WeightMatrix.getConsensus(currMotif)+", MEME";
							DrawMotifs.printMotifLogo(wm_rc, new File(imName), 75, motifLabel);
							mRCImageNames.add(imName2);
						}else{
							mImageNames.add("");
							mRCImageNames.add("");
						}
						motifImageNames.put(cond, mImageNames);
						motifRCImageNames.put(cond, mRCImageNames);
					}
				}
			}
			
			//Build up the HTML file
			
			//Header and run information 
	    	FileWriter fout = new FileWriter(htmlfilename);
	    	fout.write("<html>\n" +
	    			"\t<head><title>ChExMix results ("+config.getOutBase()+")</title></head>\n" +
	    			"\t<style type='text/css'>/* <![CDATA[ */ table, th{border-color: #600;border-style: solid;} td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>\n" +
	    			"\t<script language='javascript' type='text/javascript'><!--\nfunction motifpopitup(url) {	newwindow=window.open(url,'name','height=75');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
	    			"\t<script language='javascript' type='text/javascript'><!--\nfunction fullpopitup(url) {	newwindow=window.open(url,'name');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
	    			"\t<body>\n" +
	    			"\t<h1>ChExMix results ("+config.getOutBase()+")</h1>\n" +
	    			"");
	    	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    	Date date = new Date();
	    	fout.write("\t<p>ChExMix version "+config.version+" run completed on: "+dateFormat.format(date));
	    	fout.write(" with arguments:\n "+config.getArgs()+"\n</p>\n");
	    	
	    	
	    	//Binding event information (per condition)
	    	fout.write("\t<h2>Binding events</h2>\n" +
	    			"\t<table>\n");
	    	fout.write("\t\t<tr>" +
	    			"\t\t<th>Condition</th>\n" +
	    			"\t\t<th>Events</th>\n" +
	    			"\t\t<th>File</th>\n");
	    	if(config.getFindingMotifs())
	    		fout.write("\t\t<th>Positional Prior Motif</th>\n" +
	    				"\t\t<th>Motif Relative Offset</th>\n");
	    	fout.write("\t\t</tr>\n");
	    	for(ExperimentCondition cond : manager.getConditions()){
	    		String eventFileName=config.getOutBase()+"_"+cond.getName()+".events";
	    		if(evconfig.getEventsFileTXTExtension())
	    			eventFileName = eventFileName+".txt";
	    		fout.write("\t\t<tr>" +
		    			"\t\t<td>"+cond.getName()+"</td>\n" +
	    				"\t\t<td>"+bindingManager.countEventsInCondition(cond, evconfig.getQMinThres())+"</td>\n" +
		    			"\t\t<td><a href='"+eventFileName+"'>"+eventFileName+"</a></td>\n");
		    	if(config.getFindingMotifs()){
		    		if(motifImageNames.get(cond)!=null){
		    			for (int i=0; i< motifImageNames.get(cond).size();i++){
		    				fout.write("\t\t<td><img src='"+motifImageNames.get(cond).get(i)+"'><a href='#' onclick='return motifpopitup(\""+motifRCImageNames.get(cond).get(i)+"\")'>rc</a></td>\n" +
		    					"\t\t<td>"+bindingManager.getBindingSubtype(cond).get(i).getMotifOffset()+"</td>\n");
		    			}
		    		}
		    		else
		    			fout.write("\t\t<td>No motif found</td>\n" +
		    					"\t\t<td>NA</td>\n");
		    	}
		    	fout.write("\t\t</tr>\n");
			}fout.write("\t</table>\n");
			
	    	
	    	//Input data read counts and read distribs (per replicate)
	    	fout.write("\t<h2>Input data</h2>\n" +
	    			"\t<table>\n");
	    	fout.write("\t\t<tr>" +
	    			"\t\t<th>Replicate</th>\n" +
	    			"\t\t<th>ReadCount</th>\n" +
	    			"\t\t<th>CtrlScaling</th>\n" +
	    			"\t\t<th>SignalFraction</th>\n");
	    	fout.write("\t\t</tr>\n");
	    	for(ControlledExperiment rep : manager.getReplicates()){
				String tmpscale = rep.hasControl()?String.format("%.3f",rep.getControlScaling()):"NA";
	    		fout.write("\t\t<tr>" +
		    			"\t\t<td>"+rep.getCondName()+" "+rep.getRepName()+"</td>\n" +
	    				"\t\t<td>"+rep.getSignal().getHitCount()+"</td>\n" +
	    				"\t\t<td>"+tmpscale+"</td>\n" +
	    				"\t\t<td>"+String.format("%.3f",rep.getSignalVsNoiseFraction())+"</td>\n");
	    		fout.write("\t\t</tr>\n");
			}fout.write("\t</table>\n");
			
			//Binding subtype information (per condition)
			int maxNumSubtypes=0;
			for(ExperimentCondition cond : manager.getConditions()){
				int condNumSubtype = bindingManager.getNumBindingType(cond);
				if (condNumSubtype > maxNumSubtypes){maxNumSubtypes=condNumSubtype;}	
			}
			fout.write("\t<h2>Binding event subtypes</h2>\n" +
	    			"\t<table>\n");
			fout.write("\t\t<tr>" +
	    			"\t\t<th>Condition</th>\n" +
	    			"\t\t<th>Events</th>\n" +
	    			"\t\t<th>File</th>\n" +
	    			"\t\t<th>Heatmap</th>\n");
			if(config.getFindingMotifs())
				fout.write("\t\t<th>Sequence color plot</th>\n");
			for (int i=0; i < maxNumSubtypes; i++)
				fout.write("\t\t<th>Subtype "+i+"</th>\n");
			fout.write("\t\t</tr>\n");
			
			for(ExperimentCondition cond : manager.getConditions()){
				String subtypeEventFileName = config.getOutBase()+"_"+cond.getName()+".subtype.events";
				String heatmapFileName = "images/"+config.getOutBase()+"_"+cond.getName()+".events_"+cond.getName()+"_"+"heatmap.png";
				String seqcolorplot = "images/"+config.getOutBase()+"_"+cond.getName()+"_seq.png";
	    		fout.write("\t\t<tr>" +
		    			"\t\t<td rowspan=3>"+cond.getName()+"</td>\n" +
		    			"\t\t<td rowspan=3>"+bindingManager.countEventsInCondition(cond, evconfig.getQMinThres())+"</td>\n" +
		    			"\t\t<td rowspan=3><a href='"+subtypeEventFileName+"'>"+subtypeEventFileName+"</a></td>\n" +
		    			"\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+heatmapFileName+"\")'><img src='"+heatmapFileName+"' height='500' width='100'></a></td>\n");
				if(config.getFindingMotifs()){
					fout.write("\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+seqcolorplot+"\")'><img src='"+seqcolorplot+"' height='500' width='100'></a></td>\n");
				}
	    		String replicateName = cond.getName()+"-"+cond.getReplicates().get(0).getRepName();
	    		for (int i=0; i < maxNumSubtypes; i++){
	    			if (i < bindingManager.getNumBindingType(cond)){
	    				String distribFilename = "images/"+config.getOutBase()+"_"+replicateName+"_"+i+"_Read_Distributions.png";
	    				fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\""+distribFilename+"\")'><img src='"+distribFilename+"' height='200'></a></td>\n");
	    			}else{
	    				fout.write("\t\t<td>NA</td>\n");
	    			}					
	    		}fout.write("\t\t</tr>\n");
	    		fout.write("\t\t<tr>");
	    		
	    		if(config.getFindingMotifs()){
	    			int mc=0;
	    			int colc=0;
	    			if(!motifImageNames.get(cond).isEmpty()){
	    				for (BindingSubtype subtype :bindingManager.getBindingSubtype(cond)){
	    					if (subtype.hasMotif()){
	    						fout.write("\t\t<td><img src='"+motifImageNames.get(cond).get(mc)+"'height='100'><a href='#' onclick='return motifpopitup(\""+motifRCImageNames.get(cond).get(mc)+"\")'>rc</a></td>\n");
	    						mc++;
	    					}else{
	    						fout.write("\t\t<td>NA</td>\n");
	    					}
	    					colc++;
	    				}
	    				for (int j=colc; j<maxNumSubtypes;j++)
    						fout.write("\t\t<td>NA</td>\n");
	    			}else{
	    				for (int i=0; i < maxNumSubtypes; i++)
			    			fout.write("\t\t<td>NA</td>\n");
	    			}
	    		}else{
	    			for (int i=0; i < maxNumSubtypes; i++)
		    			fout.write("\t\t<td>NA</td>\n");
	    		}fout.write("\t\t</tr>\n");
	    		
	    		// add number of subtype specific sites
	    		for (BindingSubtype subtype :bindingManager.getBindingSubtype(cond)){
	    			fout.write("\t\t<td>NA</td>\n");
	    		}fout.write("\t\t</tr>\n");
	    		
			}fout.write("\t</table>\n");
			
			
			//File list of extras (histograms, etc)
			fout.write("\t<h2>Miscellaneous files</h2>\n");
			if(config.getFindingMotifs())
				if(evconfig.getEventsFileTXTExtension())
					fout.write("\t<p><a href='"+config.getOutBase()+".motifs.txt'>Positional prior motifs.</a> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation.</p>\n");
				else
					fout.write("\t<p><a href='"+config.getOutBase()+".motifs'>Positional prior motifs.</a> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation.</p>\n");
			fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".intraCondPeakDistances.histo.txt'>Peak-peak distance histograms (same condition)</a></p>\n");
			if(manager.getNumConditions()>1)
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".interCondPeakDistances.histo.txt'>Peak-peak distance histograms (between conditions)</a></p>\n");
			if(config.getFindingMotifs())
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".peaks2motifs.histo.txt'>Peak-motif distance histograms</a></p>\n");
	    	
	    	
	    	fout.write("\t</body>\n</html>\n");
	    	fout.close();

	    	
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void runMetaMaker(String args){
		String[] words = args.split("\\s+");
		MetaConfig mconfig = new MetaConfig(words);
		MetaMaker maker = new MetaMaker(config, gconfig, mconfig, manager);		
		maker.run();	
	}
	
}
