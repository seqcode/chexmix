package org.seqcode.projects.chexmix.stats;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.data.motifdb.WeightMatrixImport;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * Information content of set of sequences
 */
public class InformationContent {
		
	public static double[] calculateMotifIC(WeightMatrix wm){ 
		double[] ic = new double[wm.length()];
		for (int j=0; j < wm.length(); j++){
			double v = 0.0;
			for (int i=0; i < wm.matrix[j].length;i++){
				double p=wm.matrix[j][i];
				if (p > 0){ v = v - p*Math.log(p)/Math.log(2);}
			}
			ic[j]=v;
		}
		return ic;
	}
	
	public static int maxPosition(double[] array){
		int pos=0; double maxScore=0.0;
		for (int i=0; i < array.length; i++){
			if (array[i] > maxScore){
				maxScore=array[i]; pos=i;
			}
		}
		return pos;
	}
	
	public static void main(String[] args) throws ParseException{
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		List<StrandedRegion> peakReg = null;
		List<String> peakSeq = new ArrayList<String>();
		int win = Args.parseInteger(args,"win",50);
		if(ap.hasKey("peaks")){
			peakReg = RegionFileUtilities.loadStrandedRegionsFromMotifFile(gconfig.getGenome(), Args.parseString(args, "peaks", null), win);
		}else{
			System.out.println("please provide peaks !");
			System.exit(1);
		}
		
		RegionFileUtilities.getSequencesForStrandedRegions(peakReg, gconfig.getSequenceGenerator());
		WeightMatrix wm = WeightMatrixImport.buildAlignedSequenceMatrix(peakSeq);
		double[] ic = calculateMotifIC(wm);
		System.out.println(ic.toString());
		System.out.println("max IC position "+maxPosition(ic));
	}
}
