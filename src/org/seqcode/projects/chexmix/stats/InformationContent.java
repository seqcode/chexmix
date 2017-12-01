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
import org.seqcode.math.stats.StatUtil;
import org.seqcode.motifs.DrawMotifs;

/**
 * Information content of weight matrix build from aligned a set of sequence
 */
public class InformationContent { 
		
	protected double[] IC;
	protected int maxPos=0;
	protected double maxScore=0.0;
	
	public InformationContent (WeightMatrix fm, MarkovBackgroundModel back){
		
		// Get log odds version from frequency matrix
		WeightMatrix wm = new WeightMatrix(fm.length());
        //clone
        for (int i = 0; i < fm.length(); i++) {
        	wm.matrix[i]['A'] = fm.matrix[i]['A'];
        	wm.matrix[i]['C'] = fm.matrix[i]['C'];
        	wm.matrix[i]['G'] = fm.matrix[i]['G'];
        	wm.matrix[i]['T'] = fm.matrix[i]['T'];
        	wm.matrix[i]['a'] = fm.matrix[i]['a'];
        	wm.matrix[i]['c'] = fm.matrix[i]['c'];
        	wm.matrix[i]['g'] = fm.matrix[i]['g'];
        	wm.matrix[i]['t'] = fm.matrix[i]['t'];
        }
        if (!fm.islogodds) { 
        	wm.islogodds = true;
        	if (back !=null){ wm.toLogOdds();}
			else{ wm.toLogOdds(back);}
        }
        
		IC = new double[fm.length()];
		for (int i=0; i < fm.length(); i++){
			double v = 0.0;
			for (int c=0; c < wm.letters.length; c++){
				char letter = wm.letters[c];				
				double f=fm.matrix[i][letter];
				double p=wm.matrix[i][letter];
				v = v + f*p;
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
		
		double[] sic = StatUtil.gaussianSmoother(ic.getMotifIC(),1);
		double maxScore = 0.0;
		int maxPos = 0;
		for (int i=0; i < sic.length; i++){
			if (sic[i] > maxScore){
				maxScore=sic[i]; maxPos=i;
			}
		}
		System.out.println("smoothed ");
		System.out.println(Arrays.toString(sic));
		System.out.println("max pos after smoothing "+maxPos);
		
		
		// Print motif 
		String motifLabel = WeightMatrix.getConsensus(wm)+", MEME";
		DrawMotifs.printMotifLogo(wm, new File("motif.png"), 75, motifLabel);
	}
}
