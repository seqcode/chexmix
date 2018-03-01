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

import org.seqcode.genome.location.Region;


/**
 * EMStepPlotter: Used to plot a step in EM over a given region
 * @author mahony
 *
 */
public class EMStepPlotter {
	private static int trackHeight=80;
	private static int trackSpacing=10;
	private static int plotWidth = 800;
	private static int componentThickness=2;
	
	public static void execute(String outName, Region reg, int [][][][] mu, double [][] pi, double[][][][] tau, int[] numBindingType, double[][][] forCondPrior, double[][][] revCondPrior, int numConditions, int numComps, int iteration, int trimLeft, int trimRight){
		
		int [][] weightedMu = new int[numConditions][numComps];
		double newMu=0;
		for(int c=0; c<numConditions; c++){	
			for (int j=0; j <numComps;j++){
				for (int bt=0; bt <numBindingType[c];bt++){
					for (int s=0; s< 2; s++){ if(tau[c][j][bt][s]>0){
						newMu += ((double) mu[c][j][bt][s])*tau[c][j][bt][s];
					}}}
				weightedMu[c][j] = (int) Math.round(newMu);
			}}
		
		for(int c=0; c<numConditions; c++){		
			String filename = outName +"cond"+c+".png";
			File f = new File(filename);
			int w = plotWidth;
			int h = ((trackHeight*(numBindingType[c]+1))+(trackSpacing*(numBindingType[c]+2)))*2; //x2 for number of different plots
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
			int rstart = reg.getStart()+trimLeft, rend = reg.getEnd()-trimRight; 
			float rWidth = (float)(rend-rstart);
	    
			int tStart=hmargin+20;
			float trackWidth = w-(wmargin*2);
	    
			//Tau
			double tauMax = 1;
			double maxTau=0;
			for(int j=0; j<numComps; j++)
				for (int bt=0; bt <numBindingType[c];bt++)
					for (int s=0; s< 2;s++)
						if(pi[c][j]*tau[c][j][bt][s]>maxTau)
							maxTau = pi[c][j]*tau[c][j][bt][s];
			if(maxTau >0.5)
				tauMax = 1;
			else if(maxTau >0.1)
				tauMax = 0.5;
			else if(maxTau >0.01)
				tauMax = 0.1;
			else
				tauMax = 0.01;
			
			for (int bt=0; bt<numBindingType[c]; bt++){
				for (int s=0; s< 2;s++){    
					g2.setColor(Color.black);
					g2.drawLine(wmargin, tStart+trackHeight, w-wmargin, tStart+trackHeight);		// x-axis
					g2.drawLine(wmargin, tStart, wmargin, tStart+trackHeight);	// y-axis
					g2.drawString(String.format("%.3f",tauMax), wmargin+2, tStart);
	    	
					g2.setColor(Color.blue);
					g2.setStroke(componentStroke);
					for(int j=0;j<numComps;j++){ 
						if(pi[c][j]>0){
							float roffset = weightedMu[c][j]-rstart;
							if(roffset>=0 && roffset<rWidth){ //Some points may be out of bounds due to trimming
								float toffset = (roffset/rWidth)*trackWidth;
								double tauHeight =(pi[c][j]*tau[c][j][bt][s]/tauMax)*trackHeight;  
								//System.out.println("Pi: "+piHeight);
								g2.drawLine((int)(wmargin+toffset), (int)(tStart+trackHeight-tauHeight), (int)(wmargin+toffset), (int)(tStart+trackHeight));	// y-axis
							}
						}
					}
					g2.setStroke(defaultStroke);
					tStart += trackHeight+trackSpacing;
				}
			}
	    
			if(forCondPrior!=null){
				//Prior
				double priorMax = 0;
				double maxPrior=0;
				for (int bt=0; bt <numBindingType[c];bt++)
					for(int x=0; x<reg.getWidth(); x++)
						if(Math.max(forCondPrior[c][bt][x], revCondPrior[c][bt][x])>maxPrior)
							maxPrior = Math.max(forCondPrior[c][bt][x], revCondPrior[c][bt][x]);
				if(maxPrior>1)
					priorMax = maxPrior;
				else if(maxPrior >0.1)
					priorMax = 1;
				else if(maxPrior >0.01)
					priorMax = 0.1;
				else
					priorMax = 0.01;
		    
				for (int bt=0; bt<numBindingType[c]; bt++){
		    
					g2.setColor(Color.black);
					g2.drawLine(wmargin, tStart+trackHeight, w-wmargin, tStart+trackHeight);		// x-axis line
					g2.drawLine(wmargin, tStart, wmargin, tStart+trackHeight);	// y-axis line
					g2.drawString(String.format("%.3f",priorMax), wmargin+2, tStart);
		    	
					g2.setColor(Color.red);
					for(int x=0; x<rWidth; x++){
						int pos = x+trimLeft;
						if(forCondPrior[c][bt][pos]>0){
							float roffset = x;
							float toffset = (roffset/rWidth)*trackWidth;
							double priorHeight =(forCondPrior[c][bt][pos]/priorMax)*trackHeight;  
							g2.drawLine((int)(wmargin+toffset), (int)(tStart+trackHeight-priorHeight), (int)(wmargin+toffset), (int)(tStart+trackHeight));	// y-axis
						}
					}
					tStart += trackHeight+trackSpacing;
					for(int x=0; x<rWidth; x++){
						int pos = x+trimLeft;
						if(forCondPrior[c][0][pos]>0){
							float roffset = x;
							float toffset = (roffset/rWidth)*trackWidth;
							double priorHeight =(revCondPrior[c][bt][pos]/priorMax)*trackHeight;  
							g2.drawLine((int)(wmargin+toffset), (int)(tStart+trackHeight-priorHeight), (int)(wmargin+toffset), (int)(tStart+trackHeight));	// y-axis
						}
					}
					tStart += trackHeight+trackSpacing;
				}
			}
			//Iteration label:
			g2.setColor(Color.black);
			g2.setFont(new Font("Ariel", Font.BOLD, 14));
			String iString = "EM:"+iteration;
			FontMetrics metrics = g2.getFontMetrics();
			g2.drawString(iString, (w/2-(metrics.stringWidth(iString)/2)), hmargin);
	    
			try{
				ImageIO.write(im, "png", f);
			}
			catch(IOException e){
				System.err.println("Error in printing file "+filename);
			}
		} 
	}
}
