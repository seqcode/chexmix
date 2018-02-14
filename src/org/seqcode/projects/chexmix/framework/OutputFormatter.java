package org.seqcode.projects.chexmix.framework;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;


public class OutputFormatter {

	protected ChExMixConfig config;
	
	public OutputFormatter(ChExMixConfig c){ 
		config = c;	
	}
	
    /**
     * Plot a sets of read distributions, indexed in the Map by replicates 
     */
    public void plotAllReadDistributions(Map<ControlledExperiment, List<TagProbabilityDensity>> models){
		Color[] colors = {Color.black, Color.red, Color.blue, Color.green, Color.cyan, Color.orange};
		for(ControlledExperiment rep : models.keySet()){
			String replicateName = rep.getCondName()+"-"+rep.getRepName();
			
			double maxProb = 0;
		    List<TagProbabilityDensity> currModels = models.get(rep);
		    for (TagProbabilityDensity bm:currModels){
		    	double[] wpoints = bm.getWatsonProbabilities();
		    	double[] cpoints = bm.getCrickProbabilities();
		    	for (int i=0; i < wpoints.length ; i++)
		    		maxProb = Math.max(maxProb, Math.max(wpoints[i], cpoints[i]));	//max probability
		    }
			
			for (int i=0;i<currModels.size();i++){
				TagProbabilityDensity m = currModels.get(i);
				String filename = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+replicateName +"_"+i+ "_Read_Distributions.png";
				File f = new File(filename);
				int w = 500; //changed from 1000 to 500
				int h = 600;
				int margin= 50;
			    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
			    Graphics g = im.getGraphics();
			    Graphics2D g2 = (Graphics2D)g;
			    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
			    g2.setColor(Color.white);
			    g2.fillRect(0, 0, w, h);	
			    g2.setColor(Color.gray);
			    g2.drawLine(20, h-margin, w-20, h-margin);		// x-axis
			    g2.drawLine(w/2, margin, w/2, h-margin);	// y-axis    
			    g.setFont(new Font("Arial",Font.PLAIN,16));
			    for (int p=-2;p<=2;p++){
			    	g2.drawLine(w/2+p*200, h-margin-10, w/2+p*200, h-margin);	// tick  
			    	g2.drawString(""+p*200, w/2+p*200-5, h-margin+22);			// tick label
			    }
			    
			    double[] wpoints = m.getWatsonProbabilities();
		    	double[] cpoints = m.getCrickProbabilities();
			    g2.setColor(Color.blue);
			    g2.setStroke(new BasicStroke(4));
			    for (int p=0;p<wpoints.length-1;p++){
			    	int x1=p+(w-wpoints.length)/2;
			    	int y1=(int) (h-wpoints[p]/maxProb*(h-margin*2)*0.8)-margin;
			    	int x2=p+1+(w-wpoints.length)/2;
			    	int y2=(int) (h-wpoints[p+1]/maxProb*(h-margin*2)*0.8)-margin;
			    	g2.drawLine(x1, y1, x2, y2);	    
			    }
			    
			    g2.setColor(Color.red);
			    g2.setStroke(new BasicStroke(4));
			    for (int p=0;p<cpoints.length-1;p++){
			    	int x1=p+(w-cpoints.length)/2;
			    	int y1=(int) (h-cpoints[p]/maxProb*(h-margin*2)*0.8)-margin;
			    	int x2=p+1+(w-cpoints.length)/2;
			    	int y2=(int) (h-cpoints[p+1]/maxProb*(h-margin*2)*0.8)-margin;
			    	g2.drawLine(x1, y1, x2, y2);	    
			    }
			    g.setFont(new Font("Arial",Font.PLAIN,20));
			    g2.setColor(Color.black);
			    g2.drawString(String.format("%d", i), w-300, i*25+margin+25);
			    
			    try{
			    	ImageIO.write(im, "png", f);
			    }
			    catch(IOException e){
			    	System.err.println("Error in printing file "+filename);
			    }
			}
		}
	}

}
