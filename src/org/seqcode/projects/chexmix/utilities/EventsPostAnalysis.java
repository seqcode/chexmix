package org.seqcode.projects.chexmix.utilities;

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;

import javax.imageio.ImageIO;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.Pair;
import org.seqcode.gseutils.RealValuedHistogram;
import org.seqcode.motifs.DrawMotifs;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;
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
					String outFilename = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+cond.getName()+"_seq.png";
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
						File outFile = new File(outFilename);
						fig.visualizeSequences(seqs, 3, 1, outFile);
						
						// resize image
						BufferedImage seqImage = ImageIO.read(outFile);
						Image resizedImage = seqImage.getScaledInstance(Math.min(seqImage.getWidth(), 250),  Math.min(seqImage.getHeight(), 1000), Image.SCALE_DEFAULT);						
						// convert image back to buffered image
						BufferedImage bimg = new BufferedImage(Math.min(seqImage.getWidth(), 250),  Math.min(seqImage.getHeight(), 1000), BufferedImage.TYPE_INT_ARGB);
						bimg.getGraphics().drawImage(resizedImage,0,0, null);	
						ImageIO.write(bimg, "PNG", new File(outFilename));
						
						//Sequence file
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
			String pointArgs = " --peaks "+config.getOutputIntermediateDir()+File.separator+config.getOutBase()+"_"+cond.getName()+".subtype.aligned.events.txt";
						
			if(events.size()>0){			
				// Run for each strand
				System.out.println(config.getMetaMakerArgs()+pointArgs+" --strand + --color blue --noborder --out "+cond.getName()+".events");			
			
				runMetaMaker(config.getMetaMakerArgs()+pointArgs+" --strand + --color blue --noborder --out "+cond.getName()+".events");
				runMetaMaker(config.getMetaMakerArgs()+pointArgs+" --strand - --color red --noborder --out "+cond.getName()+".events");
					
				// Combine plots
				for (ExperimentCondition pcond : manager.getConditions()){				
					String pngPath=config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+pcond.getName()+".events_"+cond.getName()+"_";
					try {
										
						// load source images
						File posHeatmap = new File(pngPath+"+_lines.png");
						File negHeatmap = new File(pngPath+"-_lines.png");
						BufferedImage posImage = ImageIO.read(posHeatmap);
						BufferedImage negImage = ImageIO.read(negHeatmap);

						// create the new image, canvas size is at most 250x1000
						//BufferedImage combined = new BufferedImage(Math.min(posImage.getWidth(), 250), Math.min(posImage.getHeight(), 1000), BufferedImage.TYPE_INT_ARGB);
						// also create a new full image, canvas size unchanged
						BufferedImage combinedFull = new BufferedImage(posImage.getWidth(), posImage.getHeight(), BufferedImage.TYPE_INT_ARGB);
						
						// paint both images, preserving the alpha channels
						//Graphics g = combined.getGraphics();
						//g.drawImage(negImage, 0, 0, null);
						//g.drawImage(posImage, 0, 0, null);
						Graphics gfull = combinedFull.getGraphics();
						gfull.drawImage(negImage, 0, 0, null);
						gfull.drawImage(posImage, 0, 0, null);
					
						//((Graphics2D) g).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float) 0.6));
						((Graphics2D) gfull).setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, (float) 0.6));

						// Save as new image
						//ImageIO.write(combined, "PNG", new File(pngPath+"heatmap.png"));
						// resize image to 250 x 1000
						Image resizedImage = combinedFull.getScaledInstance(Math.min(combinedFull.getWidth(), 250),  Math.min(combinedFull.getHeight(), 1000), Image.SCALE_DEFAULT);						
						// convert image back to buffered image
						BufferedImage bimg = new BufferedImage(Math.min(combinedFull.getWidth(), 250),  Math.min(combinedFull.getHeight(), 1000), BufferedImage.TYPE_INT_ARGB);
						bimg.getGraphics().drawImage(resizedImage,0,0, null);	
						ImageIO.write(bimg, "PNG", new File(pngPath+"heatmap.png"));
					
						// delete source images
						posHeatmap.delete();
						negHeatmap.delete();
					
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}	
		//7) Report subtype model width
		try{
			String modelwidthfilename = config.getOutputIntermediateDir()+File.separator+config.getOutBase()+".subtypeModelWidth.txt";
			FileWriter fout = new FileWriter(modelwidthfilename);
			fout.write("#Subtype model width (i.e. peak width): reporting 90% probability interval distance in a model\n\n");
			for(ExperimentCondition cond : manager.getConditions()){
				fout.write("#Condition: "+cond.getName()+"\n");
				String firstline="Group\t";
				String secondline="ModelRange(bp)\t";
				double aveWidth=0;
				int[] subtypeCounts=bindingManager.countSubtypeEventsInCondition(cond, evconfig.getQMinThres());
				
				int subIndex=0;
				for (BindingSubtype sub : bindingManager.getBindingSubtype(cond)){
					TagProbabilityDensity model =sub.getBindingModel(0);
					Pair<Integer,Integer> intervals = model.probIntervalDistances(0.90);
					int longest = Math.max(Math.abs(intervals.car()), Math.abs(intervals.cdr()));
					secondline+=(longest+"\t");
					double weightedModelSize = (double) longest*(double) subtypeCounts[subIndex];
					aveWidth+=weightedModelSize;
					
					firstline+="Subtype"+subIndex+"(n="+subtypeCounts[subIndex]+")\t";
					
					subIndex++;
				}
				aveWidth/=bindingManager.countEventsInCondition(cond, evconfig.getQMinThres());
				firstline+="average\n";
				secondline+=((int)aveWidth+"\n");
								
				fout.write(firstline);
				fout.write(secondline);
			}
			fout.close();
		}catch (IOException e) {
			e.printStackTrace();
		}
		
		//8) HTML report
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
	    	
	    	if (config.useGalaxyhtml()){
		    	fout.write("<html>\n" +
		    			"\t<head><title>ChExMix results ("+config.getOutBase()+")</title></head>\n" +
		    			"\t<style type='text/css'>/* <![CDATA[ */ table, th{border-color: #600;border-style: solid;} td{border-color: #600;border-style: solid;} table{border-width: 0 0 1px 1px; border-spacing: 0;border-collapse: collapse;} th{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} td{margin: 0;padding: 4px;border-width: 1px 1px 0 0;} /* ]]> */</style>\n" +
		    			"\t<script language='javascript' type='text/javascript'><!--\nfunction motifpopitup(url) {	newwindow=window.open(url,'name',' height=75');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
		    			"\t<script language='javascript' type='text/javascript'><!--\nfunction fullpopitup(url) {	newwindow=window.open(url,'name');	if (window.focus) {newwindow.focus()}	return false;}// --></script>\n" +
		    			"\t<body>\n" +
		    			"\t<h1>ChExMix results ("+config.getOutBase()+")</h1>\n" +
		    			"");
		    	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		    	Date date = new Date();
		    	fout.write("\t<p>ChExMix version "+config.version+" run completed on: "+dateFormat.format(date));
		    	fout.write(" with arguments:\n "+config.getArgs()+"\n</p>\n");
		   
		    	
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
		    	fout.write("\t<h2>Binding event subtypes</h2>\n" +
		    			"\t<table>\n");
		    	fout.write("\t\t<tr>" +
		    			"\t\t<th>Condition</th>\n" +
		    			"\t\t<th>Events</th>\n" +
		    			"\t\t<th>File</th>\n");
		    	fout.write("\t\t</tr>\n");
		    	for(ExperimentCondition cond : manager.getConditions()){
		    		String subtypeEventFileName = config.getOutBase()+"_"+cond.getName()+".subtype.events";
		    		fout.write("\t\t<tr>" +
			    			"\t\t<td>"+cond.getName()+"</td>\n" +
		    				"\t\t<td>"+bindingManager.countEventsInCondition(cond, evconfig.getQMinThres())+"</td>\n" +
			    			"\t\t<td><a href='"+subtypeEventFileName+"'>"+subtypeEventFileName+"</a></td>\n");
			    	fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
				int maxNumSubtypes=0;
				for(ExperimentCondition cond : manager.getConditions()){
					int condNumSubtype = bindingManager.getNumBindingType(cond);
					if (condNumSubtype > maxNumSubtypes){maxNumSubtypes=condNumSubtype;}	
				}
//				fout.write("\t<h2>Binding event subtype images</h2>\n" +
//		    			"\t<table>\n");
				fout.write("\t<p></p>\n"+
						"\t<div style='overflow-x:auto;'>\n" +
						"\t<table>\n");
				fout.write("\t\t<tr>" +
		    			"\t\t<th>Condition</th>\n" +
		    			"\t\t<th>Heatmap (+/-250bp)</th>\n");
				if(config.getFindingMotifs())
					fout.write("\t\t<th>Sequence plot (+/-"+evconfig.SEQPLOTWIN+"bp)</th>\n");
				for (int i=0; i < maxNumSubtypes; i++)
					fout.write("\t\t<th>Subtype "+i+"</th>\n");
				fout.write("\t\t</tr>\n");
				
				for(ExperimentCondition cond : manager.getConditions()){
					int numEvents = bindingManager.countEventsInCondition(cond, evconfig.getQMinThres());
	    			String heatmapFileName = "images/"+config.getOutBase()+"_"+cond.getName()+".events_"+cond.getName()+"_"+"heatmap.png";
					String seqcolorplot = "images/"+config.getOutBase()+"_"+cond.getName()+"_seq.png";
					fout.write("\t\t<tr>" +
			    			"\t\t<td rowspan=3>"+cond.getName()+"</td>\n");
		    		if(numEvents>0)
		    			fout.write("\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+heatmapFileName+"\")'><img src='"+heatmapFileName+"' height='400' width='150'></a></td>\n");
		    		else
		    			fout.write("\t\t\t\t\t\t\t<td rowspan=3>No events</td>\n");
	    			if(config.getFindingMotifs()){
	    				if(numEvents>0)
	    					fout.write("\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+seqcolorplot+"\")'><img src='"+seqcolorplot+"' height='400' width='150'></a></td>\n");
	    				else
			    			fout.write("\t\t\t\t\t\t\t<td rowspan=3>No events</td>\n");
	    			}
					
		    		String replicateName = cond.getName()+"-"+cond.getReplicates().get(0).getRepName();
		    		for (int i=0; i < maxNumSubtypes; i++){
		    			if (i < bindingManager.getNumBindingType(cond)){
		    				String distribFilename = "images/"+config.getOutBase()+"_"+replicateName+"_"+i+"_Read_Distributions.png";
		    				fout.write("\t\t<td><a href='#' onclick='return fullpopitup(\""+distribFilename+"\")'><img src='"+distribFilename+"' height='300'></a></td>\n");
		    			}else{
		    				fout.write("\t\t<td>NA</td>\n");
		    			}					
		    		}fout.write("\t\t</tr>\n");
		    		fout.write("\t\t<tr>\n");
		    		
		    		if(config.getFindingMotifs()){
		    			int mc=0;
		    			int colc=0;
		    			if(!motifImageNames.get(cond).isEmpty()){
		    				for (BindingSubtype subtype :bindingManager.getBindingSubtype(cond)){
		    					if (subtype.hasMotif()){
		    						fout.write("\t\t<td><img src='"+motifImageNames.get(cond).get(mc)+"' height='70' width='250'><a href='#' onclick='return motifpopitup(\""+motifRCImageNames.get(cond).get(mc)+"\")'>rc</a></td>\n");
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
		    		fout.write("\t\t<tr>\n");
		    		// add number of subtype specific sites
		    		int[] subtypeCounts=bindingManager.countSubtypeEventsInCondition(cond, evconfig.getQMinThres());
		    		int colc=0;
		    		for (int i=0; i< bindingManager.getNumBindingType(cond);i++){
		    			fout.write("\t\t<td>"+subtypeCounts[i]+" events</td>\n");
		    			colc++;
		    		}
		    		for (int j=colc; j<maxNumSubtypes;j++)
						fout.write("\t\t<td>NA</td>\n");
		    		fout.write("\t\t</tr>\n");
		    		
				}fout.write("\t</table>\n\t</div>\n");				
				
				//Replicate consistency info	
				fout.write("\t<h2>Replicate information</h2>\n" +
		    			"\t<table>\n");
				// Replication code
		    	fout.write("\t\t<tr>" +
		    			"\t\t<th>Condition</th>\n" +
		    			"\t\t<th>Replication code file</th>\n");
		    	fout.write("\t\t</tr>\n");
		    	for(ExperimentCondition cond : manager.getConditions()){
		    		String repCodeFileName = config.getOutBase()+".all.replicationcodes.table";
		    		fout.write("\t\t<tr>" +
			    			"\t\t<td>"+cond.getName()+"</td>\n" +
			    			"\t\t<td><a href='"+repCodeFileName+"'>"+repCodeFileName+"</a></td>\n");
			    	fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
				//Replicate events
				fout.write("\t<p></p>\n"+ 
						"\t<table>\n");
				// Replication code
		    	fout.write("\t\t<tr>" +
		    			"\t\t<th>Replicate</th>\n" +
		    			"\t\t<th>Event file</th>\n");
				for (ControlledExperiment rep : manager.getReplicates()){
					String repEventFileName = config.getOutBase()+"_"+rep.getCondName()+"_"+rep.getRepName()+".repevents.txt";
					fout.write("\t\t<tr>" +
			    			"\t\t<td>"+rep.getCondName()+" "+rep.getRepName()+"</td>\n" +
			    			"\t\t<td><a href='"+repEventFileName+"'>"+repEventFileName+"</a></td>\n");
					fout.write("\t\t</tr>\n");
				}fout.write("\t</table>\n");
				
				//File list of extras (histograms, etc)
				fout.write("\t<h2>Miscellaneous files</h2>\n");
				if(config.getFindingMotifs()){
					for(ExperimentCondition cond : manager.getConditions())
						fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+"."+cond.getName()+".transfac'>"+cond.getName()+" subtype motifs</a></p>\n");
					fout.write("\t<p> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation</p>\n");
				}
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".subtypeModelWidth.txt'>Subtype model width</a></p>\n");
				fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".intraCondPeakDistances.histo.txt'>Peak-peak distance histograms (same condition)</a></p>\n");
				if(manager.getNumConditions()>1)
					fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".interCondPeakDistances.histo.txt'>Peak-peak distance histograms (between conditions)</a></p>\n");
				if(config.getFindingMotifs())
					fout.write("\t<p><a href='intermediate-results/"+config.getOutBase()+".peaks2motifs.histo.txt'>Peak-motif distance histograms</a></p>\n");
		    	
		    	
		    	fout.write("\t</body>\n</html>\n");
		    	fout.close();
	    		
	    	}else{ // bootstrap html version
	    	
	    		fout.write("<!doctype html>\n<html lang='en'>\n\t<head>\n");
	    		//Required meta tags
	    		fout.write("\t\t<meta charset='utf-8'>\n"+
	    				"\t\t<meta name='viewport' content='width=device-width, initial-scale=1, shrink-to-fit=no'>\n");
	    		//Bootstrap CSS
	    		fout.write("\t\t<link rel='stylesheet' href='https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css' integrity='sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO' crossorigin='anonymous'>\n");
	    		//google fonts
	    		fout.write("\t\t<link href='https://fonts.googleapis.com/css?family=Roboto' rel='stylesheet'>\n");
	    		//font awesome
	    		fout.write("\t\t<link rel='stylesheet' href='https://use.fontawesome.com/releases/v5.3.1/css/all.css' integrity='sha384-mzrmE5qonljUremFsqc01SB46JvROS7bZs3IO2EmfFsd15uHvIt+Y8vEf7N7fWAU' crossorigin='anonymous'>\n");
	    		//Custom styles
	    		fout.write("\t\t<style>\n"+
	    				"\t\tbody{\n\t\t\tfont-family: 'Roboto', Arial;\t\t\n\t\t}\n"+
	    				"\t\t.myimg{\n\t\t\tmax-height:400px;\n\t\t\tmax-width: 250px;\n\t\t}\n"+
	    				"\t\ti{\n\t\t\tfont-size:20px;\n\t\t}\n"+    	  
	    				"\t\t.card-header{\n\t\t\tbackground-color: #EBF5FB;\n\t\t}\n"+
	    				"\t\t</style>\n"+
	    				"\t\t<title>ChExMix</title>\n\t</head>\n\t<body>\n"+
	    				"\t\t<div class='container'>\n");
	    		//heading and command section
	    		fout.write("\t\t\t<br><h2>ChExMix Results</h2><br>\n"+
	    				"\t\t\t<div class='command' id='commandRun'>\n\t\t\t<div class='card'>\n"+
	    				"\t\t\t\t<div class='card-header' id='headingOne'>\n"+
	    				"\t\t\t\t\t<button class='btn btn-link collapsed text-dark' type='button' data-toggle='collapse' data-target='#collapseTwo' aria-expanded='false' aria-controls='collapseTwo'>\n"+
	    				"\t\t\t\t\t\t<h5 class='mb-0'> Command</h5>\n\t\t\t\t\t</button>\n\t\t\t\t\t<i class='fas fa-sort-down float-right'></i>\n"+
	    				"\t\t\t\t\t</div>\n\t\t\t\t\t<div id='collapseTwo' class='collapse' aria-labelledby='headingOne' data-parent='#commandRun'>\n");
	    		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	    		Date date = new Date();
	    		fout.write("\t\t\t\t<div class='card-body'>\n\t\t\t\t\tChExMix version "+config.version+" run completed on: "+dateFormat.format(date));
	    		fout.write(" with arguments:\n "+config.getArgs()+"\n\t\t\t\t</div>\n\t\t\t\t</div>\n\t\t\t</div>\n\t\t\t</div>\n\t\t\t<br>\n");	   
	    	
	    		//Input data read counts and read distribs (per replicate)
	    		fout.write("\t\t\t<div class='card'>\n" +
	    				"\t\t\t\t<h5 class='card-header'>Input data</h5>\n"+
	    				"\t\t\t\t<div class='card-body'>\n"+
	    				"\t\t\t\t\t<div style='overflow-x:auto;'>\n" +
	    				"\t\t\t\t\t<table class='table table-bordered'>\n");
	    		fout.write("\t\t\t\t\t\t<tr>\n"+
	    				"\t\t\t\t\t\t\t<th>Replicate</th>\n" +
	    				"\t\t\t\t\t\t\t<th>ReadCount</th>\n" +
	    				"\t\t\t\t\t\t\t<th>CtrlScaling</th>\n" +
	    				"\t\t\t\t\t\t\t<th>SignalFraction</th>\n");
	    		fout.write("\t\t\t\t\t\t</tr>\n");
	    		for(ControlledExperiment rep : manager.getReplicates()){
	    			String tmpscale = rep.hasControl()?String.format("%.3f",rep.getControlScaling()):"NA";
	    			fout.write("\t\t\t\t\t\t<tr>\n" +
	    					"\t\t\t\t\t\t\t<td>"+rep.getCondName()+" "+rep.getRepName()+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td>"+rep.getSignal().getHitCount()+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td>"+tmpscale+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td>"+String.format("%.3f",rep.getSignalVsNoiseFraction())+"</td>\n");
	    			fout.write("\t\t\t\t\t\t</tr>\n");
	    		}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n");
	    		fout.write("\t\t\t\t</div>\n\t\t\t</div>\n\t\t\t<br>\n");
			
	    		//Binding subtype information (per condition)
	    		fout.write("\t\t\t<div class='card'>\n"+
	    				"\t\t\t\t<h5 class='card-header'>Binding event subtypes</h5>\n"+
	    				"\t\t\t\t<div class='card-body'>\n");
	    		fout.write("\t\t\t\t\t<div style='overflow-x:auto;'>\n" +
	    				"\t\t\t\t\t<table class='table table-bordered'>\n"+
	    				"\t\t\t\t\t\t<tr>\n" +
	    				"\t\t\t\t\t\t\t<th>Condition</th>\n" +
	    				"\t\t\t\t\t\t\t<th>Events</th>\n" +
	    				"\t\t\t\t\t\t\t<th>File</th>\n");
	    		fout.write("\t\t\t\t\t\t</tr>\n");
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			String subtypeEventFileName = config.getOutBase()+"_"+cond.getName()+".events";
	    			fout.write("\t\t\t\t\t\t<tr>\n" +
	    					"\t\t\t\t\t\t\t<td>"+cond.getName()+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td>"+bindingManager.countEventsInCondition(cond, evconfig.getQMinThres())+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td><a href='"+subtypeEventFileName+"'>"+subtypeEventFileName+"</a></td>\n");
	    			fout.write("\t\t\t\t\t\t</tr>\n");
//				}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n\t\t\t</div>\n\t\t\t<br>\n");
	    		}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n\t\t\t\t\t<br>\n");
			
	    		// Heatmap and Motif information
	    		int maxNumSubtypes=0;
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			int condNumSubtype = bindingManager.getNumBindingType(cond);
	    			if (condNumSubtype > maxNumSubtypes){maxNumSubtypes=condNumSubtype;}	
	    		}
			
//			fout.write("\t\t\t<div class='card'>\n"+
//					"\t\t\t\t<h5 class='card-header'>Heatmap & Motif Information</h5>\n"+
//					"\t\t\t\t<div class='card-body'>\n");			
	    		fout.write("\t\t\t\t\t<div style='overflow-x:auto;'>\n" +
	    				"\t\t\t\t\t<table class='table table-bordered'>\n"+
	    				"\t\t\t\t\t\t<tr>\n"+
	    				"\t\t\t\t\t\t\t<th>Condition</th>\n" +
	    				"\t\t<th>Heatmap (+/-250bp)</th>\n");
				if(config.getFindingMotifs())
					fout.write("\t\t<th>Sequence plot (+/-"+evconfig.SEQPLOTWIN+"bp)</th>\n");
	    		for (int i=0; i < maxNumSubtypes; i++)
	    			fout.write("\t\t\t\t\t\t\t<th>Subtype "+i+"</th>\n");
	    		fout.write("\t\t\t\t\t\t</tr>\n");
			
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			int numEvents = bindingManager.countEventsInCondition(cond, evconfig.getQMinThres());
	    			String heatmapFileName = "images/"+config.getOutBase()+"_"+cond.getName()+".events_"+cond.getName()+"_"+"heatmap.png";
					String seqcolorplot = "images/"+config.getOutBase()+"_"+cond.getName()+"_seq.png";
		    		fout.write("\t\t\t\t\t\t<tr>" +
	    					"\t\t\t\t\t\t\t<td rowspan=3>"+cond.getName()+"</td>\n");
		    		if(numEvents>0)
		    			fout.write("\t\t\t\t\t\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+heatmapFileName+"\")'><img class='myimg mx-auto d-block' src='"+heatmapFileName+"' height='400' width='150'></a></td>\n");
		    		else
		    			fout.write("\t\t\t\t\t\t\t<td rowspan=3>No events</td>\n");
	    			if(config.getFindingMotifs()){
	    				if(numEvents>0)
	    					fout.write("\t\t\t\t\t\t\t<td rowspan=3><a href='#' onclick='return fullpopitup(\""+seqcolorplot+"\")'><img class='myimg mx-auto d-block' src='"+seqcolorplot+"' height='400' width='150'></a></td>\n");
	    				else
			    			fout.write("\t\t\t\t\t\t\t<td rowspan=3>No events</td>\n");
	    			}
	    			String replicateName = cond.getName()+"-"+cond.getReplicates().get(0).getRepName();
	    			for (int i=0; i < maxNumSubtypes; i++){
	    				if (i < bindingManager.getNumBindingType(cond)){
	    					String distribFilename = "images/"+config.getOutBase()+"_"+replicateName+"_"+i+"_Read_Distributions.png";
	    					fout.write("\t\t\t\t\t\t\t<td><a href='#' onclick='return fullpopitup(\""+distribFilename+"\")'><img class='myimg mx-auto d-block' src='"+distribFilename+"' height='300'></a></td>\n");
	    				}else{
	    					fout.write("\t\t\t\t\t\t\t<td>NA</td>\n");
	    				}					
	    			}fout.write("\t\t\t\t\t\t</tr>\n");
	    			fout.write("\t\t\t\t\t\t<tr>\n");
	    		
	    			if(config.getFindingMotifs()){
	    				int mc=0;
	    				int colc=0;
	    				if(!motifImageNames.get(cond).isEmpty()){
	    					for (BindingSubtype subtype :bindingManager.getBindingSubtype(cond)){
	    						if (subtype.hasMotif()){
	    							fout.write("\t\t\t\t\t\t<td><img src='"+motifImageNames.get(cond).get(mc)+"' height='70' width='250'><a href='#' onclick='return motifpopitup(\""+motifRCImageNames.get(cond).get(mc)+"\")'>rc</a></td>\n");
	    							mc++;
	    						}else{
	    							fout.write("\t\t\t\t\t\t<td>NA</td>\n");
	    						}
	    						colc++;
	    					}
	    					for (int j=colc; j<maxNumSubtypes;j++)
	    						fout.write("\t\t\t\t\t\t<td>NA</td>\n");
	    				}else{
	    					for (int i=0; i < maxNumSubtypes; i++)
	    						fout.write("\t\t\t\t\t\t<td>NA</td>\n");
	    				}
	    			}else{
	    				for (int i=0; i < maxNumSubtypes; i++)
	    					fout.write("\t\t\t\t\t\t<td>NA</td>\n");
	    			}fout.write("\t\t\t\t\t\t</tr>\n");
	    			fout.write("\t\t\t\t\t\t<tr>\n");
	    			
	    			// add number of subtype specific sites
	    			int[] subtypeCounts=bindingManager.countSubtypeEventsInCondition(cond, evconfig.getQMinThres());
	    			int colc=0;
	    			for (int i=0; i< bindingManager.getNumBindingType(cond);i++){
	    				fout.write("\t\t\t\t\t\t<td>"+subtypeCounts[i]+" events</td>\n");
	    				colc++;
	    			}
	    			for (int j=colc; j<maxNumSubtypes;j++)
	    				fout.write("\t\t\t\t\t\t<td>NA</td>\n");
	    			fout.write("\t\t\t\t\t\t</tr>\n");
	    		
	    		}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n\t\t\t\t</div>\n\t\t\t</div>\n\t\t\t<br>\n");
	    		
	    		//Replicate consistency info	    		
	    		fout.write("\t\t\t<div class='card'>\n"+
	    				"\t\t\t\t<h5 class='card-header'>Replicate information</h5>\n"+
	    				"\t\t\t\t<div class='card-body'>\n");
	    		// Replication code
	    		fout.write("\t\t\t\t\t<div style='overflow-x:auto;'>\n" +
	    				"\t\t\t\t\t<table class='table table-bordered'>\n"+
	    				"\t\t\t\t\t\t<tr>\n" +
	    				"\t\t\t\t\t\t\t<th>Condition</th>\n" +
	    				"\t\t\t\t\t\t\t<th>Replication code file</th>\n");
	    		fout.write("\t\t\t\t\t\t</tr>\n");
	    		
	    		for(ExperimentCondition cond : manager.getConditions()){
	    			String repCodeFileName = config.getOutBase()+".all.replicationcodes.table";
	    			fout.write("\t\t\t\t\t\t<tr>\n" +
	    					"\t\t\t\t\t\t\t<td>"+cond.getName()+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td><a href='"+repCodeFileName+"'>"+repCodeFileName+"</a></td>\n");
	    			fout.write("\t\t\t\t\t\t</tr>\n");
	    		}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n\t\t\t\t<br>\n");
	    		
	    		//Replicate events
	    		fout.write("\t\t\t\t\t<div style='overflow-x:auto;'>\n" +
	    				"\t\t\t\t\t<table class='table table-bordered'>\n"+
	    				"\t\t\t\t\t\t<tr>\n" +
	    				"\t\t\t\t\t\t\t<th>Replicate</th>\n" +
	    				"\t\t\t\t\t\t\t<th>Event file</th>\n");
	    		fout.write("\t\t\t\t\t\t</tr>\n");
	    		
	    		for(ControlledExperiment rep : manager.getReplicates()){
	    			String repEventFileName = config.getOutBase()+"_"+rep.getCondName()+"_"+rep.getRepName()+".repevents.txt";
	    			fout.write("\t\t\t\t\t\t<tr>\n" +
	    					"\t\t\t\t\t\t\t<td>"+rep.getCondName()+" "+rep.getRepName()+"</td>\n" +
	    					"\t\t\t\t\t\t\t<td><a href='"+repEventFileName+"'>"+repEventFileName+"</a></td>\n");
	    			fout.write("\t\t\t\t\t\t</tr>\n");
	    		}fout.write("\t\t\t\t\t</table>\n\t\t\t\t</div>\n");
	    		fout.write("\t\t\t\t</div>\n\t\t\t</div>\n\t\t\t<br>\n");
			
			
	    		//File list of extras (histograms, etc)
	    		fout.write("\t\t\t<div class='card-body'>\n");
	    		fout.write("\t\t\t\t<h4>Miscellaneous files</h4>\n");
	    		if(config.getFindingMotifs()){
	    			for(ExperimentCondition cond : manager.getConditions())
	    				fout.write("\t\t\t\t<p><a href='intermediate-results/"+config.getOutBase()+"."+cond.getName()+".transfac'>"+cond.getName()+" subtype motifs</a></p>\n");
	    			fout.write("\t\t\t\t<p> Try inputting these motifs into <a href='http://www.benoslab.pitt.edu/stamp/'>STAMP</a> for validation</p>\n");
	    		}
	    		fout.write("\t\t\t\t<p><a href='intermediate-results/"+config.getOutBase()+".subtypeModelWidth.txt'>Subtype model width</a></p>\n");
	    		fout.write("\t\t\t\t<p><a href='intermediate-results/"+config.getOutBase()+".intraCondPeakDistances.histo.txt'>Peak-peak distance histograms (same condition)</a></p>\n");
	    		if(manager.getNumConditions()>1)
	    			fout.write("\t\t\t\t<p><a href='intermediate-results/"+config.getOutBase()+".interCondPeakDistances.histo.txt'>Peak-peak distance histograms (between conditions)</a></p>\n");
	    		if(config.getFindingMotifs())
	    			fout.write("\t\t\t\t<p><a href='intermediate-results/"+config.getOutBase()+".peaks2motifs.histo.txt'>Peak-motif distance histograms</a></p>\n");
	    		fout.write("\t\t\t</div>\n\t\t");
	    	
	    		//custom javascript functions
	    		fout.write("\t\t<script language='javascript' type='text/javascript'>\n"+	    	
	    				"\t\t\tfunction fullpopitup(url) {\n"+
	    				"\t\t\t\tnewwindow = window.open(url, 'name');\n"+
	    				"\t\t\t\tif (window.focus) { newwindow.focus() }\n"+
	    				"\t\t\t\treturn false;\n\t\t\t}\n"+
	    			
	        			"\t\t\tfunction motifpopitup(url) {\n"+
	        			"\t\t\t\tnewwindow = window.open(url, 'name', ' height=75');\n"+
	        			"\t\t\t\tif (window.focus) { newwindow.focus() }\n"+
	        			"\t\t\t\treturn false;\n"+
	    				"\t\t\t}\n\t\t</script>\n");
	    		//Optional JavaScript
	    		//jQuery first, then Popper.js, then Bootstrap JS
	    		fout.write("\t\t<script src='https://code.jquery.com/jquery-3.3.1.slim.min.js' integrity='sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo' crossorigin='anonymous'></script>\n");
	    		fout.write("\t\t<script src='https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js' integrity='sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49' crossorigin='anonymous'></script>\n");
	    		fout.write("\t\t<script src='https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/js/bootstrap.min.js' integrity='sha384-ChfqqxuZUCnJSK3+MXmPNIyE6ZbWh2IMqE241rYiqJxyMiZ6OW/JmZQ5stwEULTy' crossorigin='anonymous'></script>\n");	    	
	    		fout.write("\t\t</body>\n</html>\n");
	    		fout.close();
	    	
	    	}
	    	
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
