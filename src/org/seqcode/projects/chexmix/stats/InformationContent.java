package org.seqcode.projects.chexmix.stats;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.Arrays;
import java.util.List;

import org.seqcode.data.io.BackgroundModelIO;
import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.CountsBackgroundModel;
import org.seqcode.data.motifdb.MarkovBackgroundModel;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.data.motifdb.WeightMatrixImport;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.motifs.DrawMotifs;

/**
 * Information content of weight matrix build from aligned a set of sequence
 */
public class InformationContent {
		
	protected double[] IC;
	protected int maxPos=0;
	protected double maxScore=0.0;
	
	public InformationContent (WeightMatrix wm, MarkovBackgroundModel back){
		
		// Convert frequency matrix to log odds
		if (!wm.islogodds){
			if (back !=null){ wm.toLogOdds();}
			else{ wm.toLogOdds(back);}
		}		
		IC = new double[wm.length()];
		for (int i=0; i < wm.length(); i++){
			double v = 0.0;
			for (int c=0; c < wm.letters.length; c++){
				char letter = wm.letters[c];
				double p=wm.matrix[i][letter];
				if (p > 0) { v = v - p*Math.log(p)/Math.log(2);}
			}
			IC[i]=v;
		}
		for (int i=0; i < IC.length; i++){
			if (IC[i] > maxScore){
				maxScore=IC[i]; maxPos=i;
			}
		}		
	}
	
	public double[] getMotifIC(){return IC;}
	public int getMaxPosition(){return maxPos;}
	public double getMaxScore(){return maxScore;}
		
	public static void main(String[] args) throws ParseException{
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		List<StrandedRegion> peakReg = null;
		int win = Args.parseInteger(args,"win",50);
		// Markove background model
		String markovBackModel = Args.parseString(args, "back", null);
		MarkovBackgroundModel backMod = null;    
		try {
			if(markovBackModel == null)
		    	 backMod = new MarkovBackgroundModel(CountsBackgroundModel.modelFromWholeGenome(gconfig.getGenome()));
			else
				backMod = BackgroundModelIO.parseMarkovBackgroundModel(markovBackModel, gconfig.getGenome());		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		if(ap.hasKey("peaks")){
			peakReg = RegionFileUtilities.loadStrandedRegionsFromMotifFile(gconfig.getGenome(), Args.parseString(args, "peaks", null), win);
		}else{
			System.out.println("please provide peaks !");
			System.exit(1);
		}
		
		List<String> peakSeq = RegionFileUtilities.getSequencesForStrandedRegions(peakReg, gconfig.getSequenceGenerator());
		WeightMatrix wm = WeightMatrixImport.buildAlignedSequenceMatrix(peakSeq);

		System.out.println(wm.consensus);
		for (int c=0; c < wm.letters.length; c++){
			char letter = wm.letters[c];
			for (int i=0; i < wm.length(); i++){
				System.out.print(wm.matrix[i][letter]+",");
			}
			System.out.println();
		}
		
		InformationContent ic =  new InformationContent(wm, backMod);
		System.out.println(Arrays.toString(ic.getMotifIC()));
		System.out.println("max IC position "+ic.getMaxPosition());
		
		// Print motif 
		String motifLabel = WeightMatrix.getConsensus(wm)+", MEME";
		DrawMotifs.printMotifLogo(wm, new File("motif.png"), 75, motifLabel);
	}
}
