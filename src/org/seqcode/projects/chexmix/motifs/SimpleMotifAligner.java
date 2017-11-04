package org.seqcode.projects.chexmix.motifs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Pair;
import org.seqcode.motifs.FreqMatrixImport;

public class SimpleMotifAligner {
	private int compareLength=-1;
	
	public SimpleMotifAligner(int compLen){
		compareLength = compLen;
	}
	
	
	/**
	 * Find best sub-alignment of compareLength using PCC distance metric and then extend alignment.
	 * This corresponds essentially to  Ungapped Smith-Waterman as implemented in STAMP.
	 * Query & target should be frequency matrices.
	 * Return the offset and orientation (true=forward) for the best alignment of target against query.
	 * @param query
	 * @param target
	 * @return
	 */
	public Pair<Integer, Double> align(WeightMatrix query, WeightMatrix target) {        
		double maxscore = -Double.MAX_VALUE;
		int bestTargetOff=0, bestQueryOff=0;
        int complen = compareLength;
        if (compareLength < 0) {
            complen = query.matrix.length;
        }
        if (target.matrix.length < complen)
            return new Pair<Integer,Double>(0,0.0);
        
        for (int queryoffset = 0; queryoffset <= query.matrix.length - complen ; queryoffset++) {
            for (int targetoffset = 0; targetoffset <= target.matrix.length - complen; targetoffset++) {
                double score = 0;
                for (int i = 0; i < complen; i++) {
                    score += PCC(query.getColumn(queryoffset + i), target.getColumn(targetoffset+i));
                }
                if (score > maxscore) {
                	maxscore = score;
                	bestTargetOff=targetoffset;
                	bestQueryOff = queryoffset;
                }
	        }
        }
        return new Pair<Integer,Double>((bestTargetOff-bestQueryOff),maxscore);
	}
	
	/**
	 * Pearson's correlation coefficient column scoring
	 * @param colA
	 * @param colB
	 * @return
	 */
	protected double PCC(double[] colA, double[] colB){
		double meanA=0, meanB=0;
		for(int i=0; i<colA.length; i++){
			meanA+=colA[i]; meanB+=colB[i];
		}
		meanA/=4; meanB/=4;
		
		//Numerator:
		double numer=0;
		for(int i=0; i<colA.length; i++){
			numer+= (colA[i]-meanA)*(colB[i]-meanB);
		}
		//Denom:
		double denomA=0, denomB=0;
		for(int i=0; i<colA.length; i++){
			denomA+=(colA[i]-meanA)*(colA[i]-meanA);
			denomB+=(colB[i]-meanB)*(colB[i]-meanB);
		}
		double denom = Math.sqrt(denomA*denomB);
		return (numer/denom);
	}
	
	// Load frequency matrix
	public static void main(String[] args) throws NumberFormatException, IOException{		
		ArgParser ap = new ArgParser(args);	
		String filename = ap.getKeyValue("motfile");
		FreqMatrixImport motifImport = new FreqMatrixImport();
		List<WeightMatrix> motifsets = new ArrayList<WeightMatrix>();
		motifsets.addAll(motifImport.readTransfacMatricesAsFreqMatrices(filename));
		
		int mlength=6;
		SimpleMotifAligner aligner = new SimpleMotifAligner(mlength);
		// Do pair wise-alignment
		for (int i=0; i < motifsets.size(); i++){
			for (int j=i+1; j< motifsets.size(); j++){
				Pair<Integer,Double> forAlignment = aligner.align(motifsets.get(i),motifsets.get(j));
				Pair<Integer,Double> revAlignment = aligner.align(motifsets.get(i),WeightMatrix.reverseComplement(motifsets.get(j)));
				System.out.println(motifsets.get(i).getName()+"\t"+motifsets.get(j).getName()+"\tforward: offset= "+forAlignment.car()+"\tmaxscore= "+forAlignment.cdr()/mlength);
				System.out.println(motifsets.get(i).getName()+"\t"+motifsets.get(j).getName()+"\treverse: offset= "+revAlignment.car()+"\tmaxscore= "+revAlignment.cdr()/mlength);
			}
		}
	}
}
