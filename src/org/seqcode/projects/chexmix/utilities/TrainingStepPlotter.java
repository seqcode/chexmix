package org.seqcode.projects.chexmix.utilities;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.projects.chexmix.composite.CompositeModelComponent;
import org.seqcode.projects.chexmix.composite.ProteinDNAInteractionModel;

/**
 * TrainingStepPlotter: Used to plot a step in EM/ML over a given region
 * @author mahony
 *
 */
public class TrainingStepPlotter {
	private int trackHeight=200;
	private int trackSpacing=10;
	private int plotWidth = 1060;
	int hmargin= 50, wmargin=30;
	private int componentThickness=2;
	
	public TrainingStepPlotter(){}
	
	public BufferedImage plotCompositeEM(String outName, CompositeTagDistribution composite, ProteinDNAInteractionModel model, int [] mu, double [] pi, int trainingRound, int iteration, int trimLeft, int trimRight){
		String filename=null;
		File f = null;
		if(outName!=null){
			filename = outName + ".png";
			f = new File(filename);
		}
		int w = plotWidth;
		int h = trackHeight+trackSpacing+(hmargin*2); 
		BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);
	    Stroke defaultStroke = g2.getStroke();
	    Stroke componentStroke = new BasicStroke(componentThickness);
	    
	    //Min & max locations
	    int rstart = 0+trimLeft, rend = composite.getWinSize()-trimRight;
	    float rWidth = (float)(rend-rstart);
	    Integer startCoord = rstart - model.getCenterOffset();
	    Integer endCoord = startCoord+rend-rstart;
	    int tStart=hmargin;
	    float trackWidth = w-(wmargin*2);
	    
	    //Pi
	    double piMax = 1;
	    double maxPi=0;
	    double xlPiSum = 0;
	    for(CompositeModelComponent comp : model.getXLComponents()){
	    	if(pi[comp.getIndex()]>maxPi)
	    		maxPi = pi[comp.getIndex()];
	    	xlPiSum+=pi[comp.getIndex()];
	    }
	    if(maxPi >0.5)
	    	piMax = 1;
	    else if(maxPi >0.1)
	    	piMax = 0.5;
	    else if(maxPi >0.05)
	    	piMax = 0.1;
	    else if(maxPi >0.01)
	    	piMax = 0.05;
	    else
	    	piMax = 0.01;
	    
	    //Actual plotting
    	g2.setColor(Color.black);
    	g2.setFont(new Font("Courier", Font.PLAIN, 12));
	    FontMetrics metrics = g2.getFontMetrics();
	    g2.drawLine(wmargin, tStart+trackHeight, w-wmargin, tStart+trackHeight);		// x-axis
    	g2.drawLine(wmargin, tStart, wmargin, tStart+trackHeight);	// y-axis
    	g2.drawString(String.format("%.3f",piMax), wmargin+2, tStart);
    	int center = (model.getCenterOffset()*(int)trackWidth)/model.getWidth();
    	g2.drawLine(wmargin+center, tStart+trackHeight+8, wmargin+center, tStart+trackHeight);	// 0 tick on x-axis
    	g2.drawString(new String("0"), wmargin+center-(metrics.stringWidth(new String("0"))/2), tStart+trackHeight+metrics.getHeight()+10);
    	g2.drawString(startCoord.toString(), wmargin-(metrics.stringWidth(startCoord.toString())/2), tStart+trackHeight+metrics.getHeight()+10);
    	g2.drawString(endCoord.toString(), w-wmargin-(metrics.stringWidth(endCoord.toString())/2), tStart+trackHeight+metrics.getHeight()+10);
    	
    	g2.setColor(Color.blue);
    	g2.setStroke(componentStroke);
    	//Only plotting XL components here
    	for(CompositeModelComponent comp : model.getXLComponents()){
    		int j = comp.getIndex();
    		if(pi[j]>0){
    			float roffset = mu[j]-rstart;
    			if(roffset>=0 && roffset<rWidth){ //Some points may be out of bounds due to trimming
	    			float toffset = (roffset/rWidth)*trackWidth;
	    			double piHeight =(pi[j]/piMax)*trackHeight;  
	    			g2.drawLine((int)(wmargin+toffset), (int)(tStart+trackHeight-piHeight), (int)(wmargin+toffset), (int)(tStart+trackHeight));	// y-axis
    			}
    		}
    	}
    	g2.setStroke(defaultStroke);
    	tStart += trackHeight+trackSpacing;
	    
	    
	    //Iteration label:
	    g2.setColor(Color.black);
	    g2.setFont(new Font("Courier", Font.BOLD, 14));
	    String iString = "EM:"+trainingRound+","+iteration;
	    metrics = g2.getFontMetrics();
	    g2.drawString(iString, (w/2-(metrics.stringWidth(iString)/2)), hmargin);
	    
	    //Labels for pi of components
	    g2.setFont(new Font("Courier", Font.BOLD, 14));
	    g2.setColor(Color.GRAY);
	    metrics = g2.getFontMetrics();
	    String bString = "Back : "+String.format("%.3f", pi[model.getBackgroundComponent().getIndex()]);
	    g2.drawString(bString, (3*w/4-(metrics.stringWidth(bString)+5)), hmargin);
	    g2.setColor(Color.RED);
	    String csString = "CS   : "+String.format("%.3f", pi[model.getCSComponent().getIndex()]);
	    g2.drawString(csString, (3*w/4-(metrics.stringWidth(bString)+5)), hmargin+metrics.getHeight()+5);
	    g2.setColor(Color.blue);
	    String xlString = "XLsum: "+String.format("%.3f", xlPiSum);
	    g2.drawString(xlString, (3*w/4-(metrics.stringWidth(bString)+5)), hmargin+(2*metrics.getHeight())+10);
	    
	    
	    if(f!=null){
	    	try{
		    	ImageIO.write(im, "png", f);
		    }
		    catch(IOException e){
		    	System.err.println("Error in printing file "+filename);
		    }
	    }
    	
    	return im;
	} 
	
	public BufferedImage plotPerPointML(String outName, CompositeTagDistribution composite, ProteinDNAInteractionModel model, int [][] mu, double [][] pi, int numConditions, int trainingRound, int iteration, int trimLeft, int trimRight){
		String filename = outName + ".png";
		File f = new File(filename);
		int w = plotWidth;
		int h = ((trackHeight*(numConditions+1))+(trackSpacing*(numConditions+2))); 
		int hmargin= 50, wmargin=20;
	    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);
	    Stroke defaultStroke = g2.getStroke();
	    Stroke componentStroke = new BasicStroke(componentThickness);
	    
	    //Min & max locations
	    int rstart = 0+trimLeft, rend = composite.getWinSize()-trimRight; 
	    float rWidth = (float)(rend-rstart);	    
	    int tStart=hmargin+20;
	    float trackWidth = w-(wmargin*2);
	    
	    //Pi
	    double piMax = 1;
	    double maxPi=0;
	    for(int c=0; c<numConditions; c++)
	    	for(CompositeModelComponent comp : model.getXLComponents())
		    	if(pi[c][comp.getIndex()]>maxPi)
		    		maxPi = pi[c][comp.getIndex()];
	    if(maxPi >0.5)
	    	piMax = 1;
	    else if(maxPi >0.1)
	    	piMax = 0.5;
	    else if(maxPi >0.01)
	    	piMax = 0.1;
	    else
	    	piMax = 0.01;
	    
	    for(int c=0; c<numConditions; c++){
	    
	    	g2.setColor(Color.black);
	    	g2.drawLine(wmargin, tStart+trackHeight, w-wmargin, tStart+trackHeight);		// x-axis
	    	g2.drawLine(wmargin, tStart, wmargin, tStart+trackHeight);	// y-axis
	    	g2.drawString(String.format("%.3f",piMax), wmargin+2, tStart);
	    	
	    	g2.setColor(Color.blue);
	    	g2.setStroke(componentStroke);
	    	//Only plotting XL components here
	    	for(CompositeModelComponent comp : model.getXLComponents()){
	    		int j = comp.getIndex();
	    		if(pi[c][j]>0){
	    			float roffset = mu[c][j]-rstart;
	    			if(roffset>=0 && roffset<rWidth){ //Some points may be out of bounds due to trimming
		    			float toffset = (roffset/rWidth)*trackWidth;
		    			double piHeight =(pi[c][j]/piMax)*trackHeight;  
		    			g2.drawLine((int)(wmargin+toffset), (int)(tStart+trackHeight-piHeight), (int)(wmargin+toffset), (int)(tStart+trackHeight));	// y-axis
	    			}
	    		}
	    	}
	    	g2.setStroke(defaultStroke);
	    	tStart += trackHeight+trackSpacing;
	    }
	    
	    //Iteration label:
	    g2.setColor(Color.black);
	    g2.setFont(new Font("Courier", Font.BOLD, 14));
	    String iString = "EM:"+trainingRound+","+iteration;
	    FontMetrics metrics = g2.getFontMetrics();
	    g2.drawString(iString, (w/2-(metrics.stringWidth(iString)/2)), hmargin);
	    
	    //Labels for pi of Background and CS components
	    g2.setFont(new Font("Courier", Font.BOLD, 14));
	    g2.setColor(Color.GRAY);
	    metrics = g2.getFontMetrics();
	    String bString = "Back: "+String.format("%.3f", pi[model.getBackgroundComponent().getIndex()]);
	    g2.drawString(bString, (3*w/4-(metrics.stringWidth(bString)+5)), hmargin);
	    g2.setColor(Color.RED);
	    String csString = "CS  : "+String.format("%.3f", pi[model.getBackgroundComponent().getIndex()]);
	    g2.drawString(csString, (3*w/4-(metrics.stringWidth(bString)+5)), hmargin+metrics.getHeight()+5);
	    
    	try{
	    	ImageIO.write(im, "png", f);
	    }
	    catch(IOException e){
	    	System.err.println("Error in printing file "+filename);
	    }
    	
    	return im;
	} 
}
