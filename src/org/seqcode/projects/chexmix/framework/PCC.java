package org.seqcode.projects.chexmix.framework;

import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;
import org.seqcode.math.stats.StatUtil;

/**
 * Calculates Pearson correlation with sliding window
 */
public class PCC {
	protected ChExMixConfig config;
	protected int numReps;
	protected int profileWidth;
	protected double[][] profile_w;	//merged watson profile
	protected double[][] profile_c;	//merged crick profile
	protected double[][] s_profile_plus_b;	//shifted watson profile
	protected double[][] s_profile_minus_b;	//shifted crick profile
	protected double maxScore;
	protected double klScore;
	protected int maxOffset;
	protected boolean isReverse;
	
	public PCC (ChExMixConfig config, int numReps){
		this.config = config;
		profileWidth=config.MAX_BINDINGMODEL_WIDTH;
		this.numReps = numReps;
		profile_w = new double[numReps][profileWidth];
		profile_c = new double[numReps][profileWidth];
		s_profile_plus_b = new double[numReps][profileWidth];
		s_profile_minus_b = new double[numReps][profileWidth];
	}
	
	public double getCorrelation(){ return maxScore;}
	public double getDivergenceKLScore(){return klScore;}
	public int getOffset(){ return maxOffset;}	
	public boolean isReverse(){return isReverse;}
	public double[][] getMergedWatsonProfile(){return profile_w;}
	public double[][] getMergedCrickProfile(){return profile_c;}
	public double[][] getShiftedWatsonProfile(){return s_profile_plus_b;}
	public double[][] getShiftedCrickProfile(){return s_profile_minus_b;}
	
	public void execute(double[][] profile_plus_a, double[][] profile_minus_a, double[][] profile_plus_b, double[][] profile_minus_b){		
		
		maxScore = -1; maxOffset = 0; isReverse = false; klScore =-1;
		// initialize shifted profile
		for (int r=0; r < numReps; r++){
			for (int w=0; w< profileWidth; w++){
				s_profile_plus_b[r][w]=0.0; s_profile_minus_b[r][w]=0.0;
			}
		}

		// Array A is fixed, concatenate to a single array
		double[][] profie_a = new double[numReps][profileWidth*2];
		double[] aveA = new double[numReps];
		for (int r=0; r< numReps ; r++){
			for (int w=0; w< profileWidth; w++){
				profie_a[r][w] = profile_plus_a[r][w];
				profie_a[r][profileWidth+w] = profile_minus_a[r][w];
			}
			//Normalize and get average
			double sumA=0, sumAn=0;
			for (int i=0; i< profie_a[r].length; i++)
				sumA += profie_a[r][i];
			for (int i=0; i< profie_a[r].length; i++)
				profie_a[r][i]/=sumA;
			for (int i=0; i< profie_a[r].length; i++)
				sumAn += profie_a[r][i];
			aveA[r] = sumAn/((double) profie_a[r].length);
		}
	
		for (int sw = 0; sw < config.PCC_SLIDING_WINDOW; sw++){	
			double[][] profile_b = new double[numReps][profileWidth*2];
			double[][] profile_br = new double[numReps][profileWidth*2];
			double[] aveB = new double[numReps];
			int offset=sw-(config.PCC_SLIDING_WINDOW/2);
			for (int r=0; r< numReps ; r++){
				for (int i=0; i < profile_b[r].length; i++){
					profile_b[r][i]=0; profile_br[r][i]=0;
				}
				for (int i=0; i< profileWidth; i++){
					if (offset+i >=0 && offset+i < profile_plus_b[r].length){
						profile_b[r][i]=profile_plus_b[r][offset+i];
						profile_b[r][profileWidth+i]=profile_minus_b[r][offset+i];
						profile_br[r][profileWidth-i-1]=profile_minus_b[r][offset+i];
						profile_br[r][profileWidth*2-i-1]=profile_plus_b[r][offset+i];
					}
				}
				// Normalize and get average				
				double sumB=0, sumBn=0;
				for (int i=0; i< profile_b[r].length;i++)
					sumB += profile_b[r][i];
				for (int i=0; i < profile_b[r].length; i++){
					profile_b[r][i]/=sumB;
					profile_br[r][i]/=sumB;
				}
				for (int i=0; i< profile_b[r].length;i++)
					sumBn += profile_b[r][i];
				aveB[r] = sumBn/((double) profile_b[r].length);
			}
		
			// calculate Pearson correlation from forward and reverse pair
			double covf = 0, covr = 0;
			double varA = 0, varBf = 0, varBr = 0;
			for (int r=0; r< numReps ; r++){
				for (int i = 0; i < profileWidth; i++){
					double ai = profie_a[r][i] - aveA[r];
					double bi_f = profile_b[r][i] - aveB[r];
					double bi_r = profile_br[r][i] - aveB[r];				
					covf += ai*bi_f;
					covr += ai*bi_r;
					varA += ai*ai;
					varBf += bi_f*bi_f;
					varBr += bi_r*bi_r;
				}	
			}
			double corrF = covf/(Math.sqrt(varA)*Math.sqrt(varBf));
			double corrR = covr/(Math.sqrt(varA)*Math.sqrt(varBr));		
			
			if (corrF > maxScore){ maxScore = corrF; maxOffset = offset; isReverse = false; }
			if (corrR > maxScore){ maxScore = corrR; maxOffset = offset-1; isReverse = true;}
		}
	
		if (!isReverse){
//			System.out.println("forward has higher correlation "+maxScore+" offset "+maxOffset);
			for (int r=0; r<numReps; r++){
				for (int i=0; i< profileWidth;i++){
					if ((i+maxOffset) >= 0 && (i+maxOffset) < profileWidth){
						s_profile_plus_b[r][i] = profile_plus_b[r][i+maxOffset];
						s_profile_minus_b[r][i] = profile_minus_b[r][i+maxOffset];
					}}}						
		}else{
//			System.out.println("Reverse has higher correlation "+maxScore+" offset "+maxOffset);
			for (int r=0; r<numReps; r++){
				for (int i=0; i< profileWidth;i++){
					if ((i-maxOffset) > 0 && (i-maxOffset) <= profileWidth){
						// need to check this !!					
						s_profile_plus_b[r][i] = profile_minus_b[r][profileWidth-i+maxOffset];	
						s_profile_minus_b[r][i] = profile_plus_b[r][profileWidth-i+maxOffset];
					}}}
		}	
		
		//create merged profile
		for (int r=0; r<numReps; r++){
			for (int i=0; i< profileWidth; i++){
				profile_w[r][i]=profile_plus_a[r][i]+s_profile_plus_b[r][i];				
				profile_c[r][i]=profile_minus_a[r][i]+s_profile_minus_b[r][i];
			}
		}
		
		double kl=0.0;
		for (int r=0; r< numReps ; r++){
			
			//trim trailing zeros
			int nz_start=0;
			int nz_end=config.MAX_BINDINGMODEL_WIDTH-1;
			for (int i=0; i < config.MAX_BINDINGMODEL_WIDTH; i++ ){
				if (profile_plus_a[r][i]+profile_minus_a[r][i]+s_profile_plus_b[r][i]+s_profile_minus_b[r][i]>0){
					nz_start=i;break;
				}
			}
			for (int i=config.MAX_BINDINGMODEL_WIDTH-1; i <= 0; i-- ){
				if (profile_plus_a[r][i]+profile_minus_a[r][i]+s_profile_plus_b[r][i]+s_profile_minus_b[r][i]>0){
					nz_end=i;break;
				}
			}	
			// copy nonzero tags
			int nonzeroWidth=nz_end-nz_start+1;
			double[] nz_aW=new double[nonzeroWidth];
			double[] nz_aC=new double[nonzeroWidth];
			double[] nz_bW=new double[nonzeroWidth];
			double[] nz_bC=new double[nonzeroWidth];
			for (int i=0; i< nonzeroWidth; i++){
				nz_aW[i]=profile_plus_a[r][i+nz_start];
				nz_aC[i]=profile_minus_a[r][i+nz_start];
				nz_bW[i]=s_profile_plus_b[r][i+nz_start];
				nz_bC[i]=s_profile_minus_b[r][i+nz_start];
			}
					
			kl += StatUtil.KL_Divergence(nz_aW, nz_bW) + StatUtil.KL_Divergence(nz_bW, nz_aW);
			kl += StatUtil.KL_Divergence(nz_aC, nz_bC) + StatUtil.KL_Divergence(nz_bC, nz_aC);
		}
		klScore = kl;		
	}
	
	/**
	 * PCC test
	 * @param args
	 */
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		econ.setPerBaseReadFiltering(false);
		econ.setLoadRead2(false);
		ChExMixConfig config = new ChExMixConfig(gcon, args);
		if(args.length==0){
			System.err.println("CompositeTagDistribution:"+
					"\t--cpoints <stranded point file>"+
					"\t--cwin <window around points>"+
					"\t--out <output file name>"+
					"Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			// parse command line arguments
			int w = Args.parseInteger(args, "cwin", 500);
			String pFile = Args.parseString(args, "cpoints", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
				
			CompositeTagDistribution maker = new CompositeTagDistribution(pts, manager, w, true);
			
			StrandedPoint firstp=pts.get(0);
			StrandedPoint currpt=pts.get(1);
			double[][] watsonA=new double [manager.getReplicates().size()][];
			double[][] crickA=new double [manager.getReplicates().size()][];
			double[][] watsonB=new double [manager.getReplicates().size()][];
			double[][] crickB=new double [manager.getReplicates().size()][];
			for (ControlledExperiment rep : manager.getReplicates()){
				int x=rep.getIndex();
				watsonA[x]=maker.getPointWatson(firstp, rep.getCondition());
				crickA[x]=maker.getPointCrick(firstp, rep.getCondition());
				watsonB[x]=maker.getPointWatson(currpt, rep.getCondition());
				crickB[x]=maker.getPointCrick(currpt, rep.getCondition());
			}
			System.out.println("profile length "+watsonA[0].length);
			
			for (ExperimentCondition cond: manager.getConditions()){
				PCC pcc = new PCC(config, cond.getReplicates().size());
				pcc.execute(watsonA, crickA, watsonB, crickB);				
			}
			manager.close();
		}
	}
}
