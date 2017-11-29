package org.seqcode.projects.chexmix.stats;

import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.data.motifdb.WeightMatrixImport;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.gseutils.ArgParser;

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
		List<Region> peakReg = null;
		List<String> peakSeq = new ArrayList<String>();
		if(ap.hasKey("peaks")){
			peakReg = RegionFileUtilities.loadRegionsFromPeakFile(gconfig.getGenome(), "peaks", 50);
		}else{
			System.out.println("please provide peaks !");
			System.exit(1);
		}
		SequenceGenerator seqgen = gconfig.getSequenceGenerator();
		for (Region reg : peakReg){peakSeq.add(seqgen.execute(reg));}
		WeightMatrix wm = WeightMatrixImport.buildAlignedSequenceMatrix(peakSeq);
		double[] ic = calculateMotifIC(wm);
		System.out.println(ic.toString());
		System.out.println("max IC position "+maxPosition(ic));
	}
}
