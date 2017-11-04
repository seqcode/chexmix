package org.seqcode.projects.chexmix.events;


import java.io.FileWriter;
import java.io.IOException;

import java.util.List;

import org.seqcode.gseutils.Pair;
import org.seqcode.projects.chexmix.composite.TagProbabilityDensity;

public class TagDensity extends TagProbabilityDensity{
	
	public TagDensity(int size) {
		super(size);
		// TODO Auto-generated constructor stub
	}
	public TagDensity(List<Pair<Integer, Double>> empiricalWatson, List<Pair<Integer, Double>> empiricalCrick) {
		super(empiricalWatson, empiricalCrick);
		// TODO Auto-generated constructor stub
	}
	
	//Update the influence range
	public void updateInfluenceRange(){
		Pair<Integer,Integer> intervals = probIntervalDistances(0.75);
		int longest = Math.max(Math.abs(intervals.car()), Math.abs(intervals.cdr()));
		influenceRange = longest;
	}
	
	public void printDensityToFile(String filename){
		try {
			FileWriter writer = new FileWriter(filename);
			writer.write("#watson\n");
			for (int i=0; i <watsonProbs.length ; i++)
				writer.write(watsonProbs[i]+",");
			writer.write("\n");
			writer.write("#crick\n");
			for (int i=0; i< crickProbs.length; i++)
				writer.write(crickProbs[i]+",");
			writer.write("\n");
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
	
