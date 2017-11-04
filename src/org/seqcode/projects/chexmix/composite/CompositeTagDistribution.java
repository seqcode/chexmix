package org.seqcode.projects.chexmix.composite;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;


/**
 * CompositeTagDistribution: watson/crick tag distributions for a collection of aligned points and the resulting composite
 *
 * 	Coordinates are 0-based, and the center of the distribution/alignment is defined by the centerOffset variable	
 * @author mahony
 *
 */
public class CompositeTagDistribution {
	protected ExperimentCondition condition;	
	protected List<StrandedPoint> points;
	protected int win;
	protected int centerOffset;
	protected int numPoints;
	protected double[] watson; //per-condition watson tags {condition, location}
	protected double[] crick;  //per-condition crick tags  {condition, location}
	protected double[][] perPointWatson; //per-point, per-condition watson tags  {point, condition, location}
	protected double[][] perPointCrick;  //per-point, per-condition crick tags   {point, condition, location}
	protected HashMap<StrandedPoint,Integer> pointIndex = new HashMap<StrandedPoint,Integer>();
	protected boolean isSignal;
	
	public CompositeTagDistribution(List<StrandedPoint> points, ExperimentCondition cond, int win, boolean loadSignal){
		condition = cond;
		this.win = win;
		centerOffset = win/2;
		this.points = points;
		numPoints = points.size();
		isSignal = loadSignal;
	
		watson = new double[win];
		crick = new double[win];
		perPointWatson = new double[numPoints][win];
		perPointCrick = new double[numPoints][win];
		
		for(int p=0; p<numPoints; p++)
			pointIndex.put(points.get(p), p);

		//Reset
		for(int w=0; w<win; w++){watson[w]=0; crick[w]=0;}
		for(int p=0; p<numPoints; p++)
			for(int w=0; w<win; w++){
				perPointWatson[p][w]=0; perPointCrick[p][w]=0;
			}
		
		for(ControlledExperiment rep : cond.getReplicates()){
				
			if(loadSignal || rep.hasControl()){
				//Iterate through points
				int p=0;
				for(StrandedPoint pt : points){
					//Load reads
					List<StrandedBaseCount> wReads = loadSignal ? 
							rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand()) : 
								rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand());
					List<StrandedBaseCount> cReads = loadSignal ? 
							rep.getSignal().getStrandedBases(pt.expand(win), pt.getStrand()=='+' ? '-' : '+') :
								rep.getControl().getStrandedBases(pt.expand(win), pt.getStrand()=='+' ? '-' : '+');
						
						
					if(pt.getStrand()=='+'){
						for(StrandedBaseCount sbc : wReads){
							int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
							if(sdist>=0 && sdist<win){
								watson[sdist]+=sbc.getCount();
								perPointWatson[p][sdist]+=sbc.getCount();
							}
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = sbc.getCoordinate()-pt.getLocation()+(win/2);
							if(sdist>=0 && sdist<win){
								crick[sdist]+=sbc.getCount();
								perPointCrick[p][sdist]+=sbc.getCount();
							}
						}
					}else{
						for(StrandedBaseCount sbc : wReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win){
								watson[sdist]+=sbc.getCount();
								perPointWatson[p][sdist]+=sbc.getCount();
							}
						}
						for(StrandedBaseCount sbc : cReads){
							int sdist = pt.getLocation()-sbc.getCoordinate()+(win/2);
							if(sdist>=0 && sdist<win){
								crick[sdist]+=sbc.getCount();
								perPointCrick[p][sdist]+=sbc.getCount();
							}
						}
					}						
					p++;
				}
			}
			//Normalize
			double wsum=0, csum=0;
			for(int w=0; w<win; w++){
				wsum+=watson[w]; csum+=crick[w];
			}for(int w=0; w<win; w++){
				watson[w]/=wsum; crick[w]/=csum;
			}
		}
	}
	
	//Accessors
	public int getWinSize(){return win;}
	public int getCenterOffset(){return centerOffset;}
	public double[] getCompositeWatson(){return watson;}
	public double[] getCompositeCrick(){return crick;}
	public double[] getPointWatson(StrandedPoint p){return perPointWatson[pointIndex.get(p)];}
	public double[] getPointCrick(StrandedPoint p){return perPointCrick[pointIndex.get(p)];}
	public double[] getPointWatsons(int index){return perPointWatson[index];}
	public double[] getPointCricks(int index){return perPointCrick[index];}
	public List<StrandedPoint> getPoints(){return points;}
	public StrandedPoint getPoint(int i){return points.get(i);}
	
	/**
	 * Per-condition sum of tags in composites
	 * @return
	 */
	public double getCompositeSums(){
		double sums=0;
		for(int i=0; i<win; i++)
			sums+=watson[i]+crick[i];

		return sums;
	}
	
	public String toString(){
		String out="";
		for(int w=0; w<win; w++){
			int pos = (w-centerOffset);
			out = out + pos+"\t"+watson[w]+"\t"+crick[w]+"\n";
		}
		return out;
	}
		
	//Print probs to a file
	public void printProbsToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int w=0; w<win; w++){
				int pos = (w-centerOffset);
				fout.write(pos+"\t"+watson[w]+"\t"+crick[w]+"\n");
			}
			fout.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	//Main method to make new composite distributions
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if(args.length==0){
			System.err.println("CompositeTagDistribution:"+
					"\t--cpoints <stranded point file>"+
					"\t--cwin <window around points>"+
					"Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			ExperimentManager manager = new ExperimentManager(econ);
			
			int w = Args.parseInteger(args, "cwin", 400);
			String pFile = Args.parseString(args, "cpoints", null);
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
						
			for(ExperimentCondition cond : manager.getConditions()){
				CompositeTagDistribution maker = new CompositeTagDistribution(pts, cond, w, true);
				String compositeFileName = "out_composite."+cond.getName()+".txt";
				maker.printProbsToFile(compositeFileName);
			}
			manager.close();
		}
	}
	
}

