[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexmix/badges/installer/conda.svg)](https://anaconda.org/bioconda/chexmix) 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexmix/badges/downloads.svg)](https://anaconda.org/bioconda/chexmix)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexmix/badges/platforms.svg)](https://anaconda.org/bioconda/chexmix)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/chexmix/badges/license.svg)](https://anaconda.org/bioconda/chexmix)

# ChExMix
ChExMix: the ChIP-exo mixture model

ChExMix aims to characterize protein-DNA binding subtypes in ChIP-exo experiments. ChExMix assumes that different regulatory complexes will result in different protein-DNA crosslinking signatures in ChIP-exo data, and thus analysis of ChIP-exo sequencing tag patterns should enable detection of multiple protein-DNA binding modes for a given regulatory protein. ChExMix uses a mixture modeling framework to probabilistically model the genomic locations and subtype membership of protein-DNA binding events, leveraging both ChIP-exo tag enrichment patterns and DNA sequence information. In doing so, ChExMix offers a more principled and robust approach to characterizing binding subtypes than simply clustering binding events using motif information.

Citation:
--------------
N Yamada, WKM Lai, N Farrell, BF Pugh, S Mahony .“Characterizing protein-DNA binding event subtypes in ChIP-exo data”. Bioinformatics (2019) 35(6):903-913. [doi:10.1093/bioinformatics/bty703](http://dx.doi.org/10.1093/bioinformatics/bty703).
This paper was presented at [RECOMB 2018](http://recomb2018.fr/).

Downloading Executables
--------------
Executable available from: http://mahonylab.org/software/chexmix 

Building from Source
--------------
If you want to build the code yourself, you will need to first download and build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

Dependencies:
--------------
1. ChExMix requires Java 8+. 
2. You need [MEME](http://meme-suite.org/) installed and being available in $PATH if you want to find subtypes using motif (tested with MEME version 4.11.3).
2. ChExMix loads all data to memory, so you will need a lot of available memory if you are running analysis over many conditions or large datasets.

Running ChExMix
--------------
Note that ChExMix performs a lot of EM optimization of binding events along the genome, and also integrates motif-finding via MEME. Therefore, expect it to be very time and memory intensive.

Running from a jar file:

```{r, engine='sh', count_lines}
java -Xmx20G -jar chexmix.jar <options - see below>
```

In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. If you have installed source code from github, and if all classes are in your CLASSPATH, you can run ChExMix as follows:

```{r, engine='sh', count_lines}
java -Xmx20G org.seqcode.projects.chexmix.ChExMix <options - see below>
```

Options (Required/important options are in __bold__)

__General__:

  * --__out__ \<prefix>: Output file prefix. All output will be put into a directory with the prefix name. 
  * --threads \<n\>:  Use n threads during binding event detection (default=1)
  * --verbose: Flag to print intermediate files and extra output.

__Specifying the Genome__:

  * --__geninfo__ \<genome info file\>: This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ \<file\/path\> : A fasta format file or a directory containing fasta format files corresponding to every named chromosome is required if you want to find subtypes via motifs or use a motif-prior within ChExMix.
  * --__back__ \<path\> : A file containing Markov background model for the genome is required if you want to run motif-finding or use a motif-prior within ChExMix.
  
      The Markov background model files for some species:  
      | [human](http://lugh.bmb.psu.edu/software/chexmix/backgrounds/human.back) | [mouse](http://lugh.bmb.psu.edu/software/chexmix/backgrounds/mouse.back) | [Drosophila](http://lugh.bmb.psu.edu/software/chexmix/backgrounds/fly.back) | [S. cerevisiae](http://lugh.bmb.psu.edu/software/chexmix/backgrounds/yeast.back) | [E. coli](http://lugh.bmb.psu.edu/software/chexmix/backgrounds/ecoli.back) |

__Loading Data__:

  * --__exptCONDNAME-REPNAME__ \<file\>: Defines a file containing reads from a signal experiment. Replace CONDNAME and REPNAME with appropriate condition and replicate labels.
  * --__ctrlCONDNAME-REPNAME__ \<file\>: Optional arguments. Defines a file containing reads from a control experiment. Replace CONDNAME and REPNAME with appropriate labels to match a signal experiment (i.e. to tell ChExMix which condition/replicate this is a control for). If you leave out a REPNAME, this file will be used as a control for all replicates of CONDNAME.  
  * --__format__ \<SAM/BAM/BED/IDX\>: Format of data files. All files must be the same format if specifying experiments on the command line. Supported formats are SAM/BAM, BED, and IDX index files.
 
Instead of using the above options to specify each and every ChIP-seq data file on the command-line, you can instead use a design file:
 
  * --__design__ \<file\>: A file that specifies the data files and their condition/replicate relationships. See [here](http://lugh.bmb.psu.edu/software/multigps/example.design) for an example design file. The file should be formatted to contain the following pieces of information for each data file, in this order and tab-separated:
  
    * File name
    * Label stating if this experiment is “signal” or “control”
    * File format (SAM/BAM/BED/IDX) – mixtures of formats are allowed in design files
    * Condition name
    * Replicate name (optional for control experiments – if used, the control will only be used for the corresponding named signal replicate)
    
 Limits on how many reads can have their 5′ end at the same position in each replicate. Please use the options if your experiments contain many PCR duplicates:

 * --readfilter: Flag to turn on filtering reads, recommended for highly duplicated experiments. It estimates a global per-base limit from a Poisson distribution parameterized by the number of reads divided by the number of mappable bases in the genome. The per-base limit is set as the count corresponding to the 10^-7 probability level from the Poisson. Default = no read filter
 * --fixedpb \<value\>: Fixed per-base limit.
 * --poissongausspb \<value\>: Filter per base using a Poisson threshold parameterized by a local Gaussian sliding window (i.e. look at neighboring positions to decide what the per-base limit should be). 
 * --nonunique: Flag to use non-unique reads. 
 * --mappability \<value\>: Fraction of the genome that is mappable for these experiments. Default=0.8.
 * --nocache: Flag to turn off caching of the entire set of experiments (i.e. run slower with less memory)
 
 __Scaling Data__:
 
 * --noscaling: Flag to turn off auto estimation of signal vs control scaling factor.
 * --medianscale: Flag to use scaling by median ratio of binned tag counts. Default = scaling by NCIS.
 * --regressionscale: Flag to use scaling by regression on binned tag counts. Default = scaling by NCIS.
 * --sesscale: Flag to use scaling by SES (Diaz, et al. Stat Appl Genet Mol Biol. 2012).
 * --fixedscaling \<scaling factor\>: Multiply control counts by total tag count ratio and then by this factor. Default: scaling by NCIS.
 * --scalewin \<window size\>: Window size for estimating scaling ratios. Default is 10Kbp. Use something much smaller if scaling via SES (e.g. 200bp).
 * --plotscaling: Flag to plot diagnostic information for the chosen scaling method.
 
__Running ChExMix__:

  * --round \<int\>: Max. model update rounds (default=3).
  * --nomodelupdate: Flag to turn off binding model updates.
  * --minmodelupdateevents \<int\>: Minimum number of events to support an update (default=50)
  * --prlogconf \<value\>: Poisson log threshold for potential region scanning (default=-6)
  * --fixedalpha \<int\>: Impose this alpha (default: set automatically). The alpha parameter is a sparse prior on binding events in the ChExMix model. It can be interpreted as a minimum number of reads that each binding event must be responsible for in the model. 
  * --alphascale \<value\>: Alpha scaling factor (default=1.0). Increasing this parameter results in stricter binding event calls.
  * --betascale \<value\>: Beta scaling factor (default=0.05). The beta parameter is a sparse prior on binding event subtype assignment in the ChExMix model. Increasing this parameter may result in inaccurate subtype assignment.
  * --epsilonscale \<value\>: Epsilon scaling factor (default=0.2). The epsilon parameter control a balance between motif and read distribution in subtype assignment. Increasing this parameter will increase the contribution of motif and decreases the contribution of read distribution.
  * --peakf \<file\>: File of peaks to initialize component positions. See [here](https://github.com/seqcode/chexmix/blob/master/manuscriptsppl/MCF7-Erpos.peaks) for an example peak file.
  * --galaxyhtml: Flag to produce a html output appropreate for galaxy
  * --__exclude__ \<file\>: File containing a set of regions to ignore during ChExMix training. It’s a good idea to exclude the mitochondrial genome and other ‘blacklisted’ regions that contain artifactual accumulations of reads in both ChIP-exo and control experiments. ChExMix will waste time trying to model binding events in these regions, even though they will not typically appear significantly enriched over the control (and thus will not be reported to the user). See the format of an exclude region file [here](http://lugh.bmb.psu.edu/software/multigps/support/mm9_excludes.txt) (example for mm9).
  * --excludebed \<file\>: Bed file containing a set of regions to ignore during ChExMix training.


__Binding event reporting mode__: (controls which events to put into .events file)
  * --standard [report events that pass significance threshold in condition as a whole (default mode)]
  * --lenient [report events that pass significance in >=1 replicate *or* the condition as a whole.]
  * --lenientplus [report events that pass significance in condition OR (>=1 replicate AND no signif diff between replicates)]

	
__Finding ChExMix subtypes__:

1. Using motifs:

  * --motfile \<file\>: File of motifs in transfac format to initialize subtype motifs
  * --__memepath__  \<path\>: Path to the meme bin dir (default: meme is in $PATH). MEME path is required for motif finding.
  * --nomotifs \<value\>: Flag to turn off motif-finding & motif priors
  * --nomotifprior \<value\>: Flag to turn off motif priors only
  * --mememinw \<value\>: minw arg for MEME (default=6).
  * --mememaxw \<value\>: maxw arg for MEME (default=18).
  * --memenmotifs \<int\>: Number of motifs MEME should find in each condition (default=3)
  * --memeargs \<args\>: Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)
  * --minroc \<value\>: Minimum motif ROC value (default=0.7)
  * --seqrmthres \<value\>: Filter out sequences with motifs below this threshold for recursively finding motifs (default=0.1)
  * --minmodelupdaterefs \<int\>: Minimum number of motif reference to support an subtype distribution update (default=25)
 
 2. Using ChIP-exo tag distributions:

  * --noclustering: Flag to turn off read distribution clustering
  * --pref \<value\>: Preference value for read distribution clustering (default=-0.1)
  * --numcomps \<int\>: Number of components to cluster (default=500)
  * --win \<int\>: Read profile window size (default=150)
  
__Reporting binding events__:

  * --q \<value\>: Q-value minimum (default=0.01)
  * --minfold \<value\>: Minimum event fold-change vs scaled control (default=1.5)

Example
--------------
This example runs ChExMix v0.1 on a simulated dataset. The simulated data is made by mixing mutually-exclusive binding events from yeast Reb1 and Abf1 ChIP-exo datasets. This ensures that there are two distinct binding event subtypes in the dataset. The version of ChExMix and all files required to run this analysis are in this file: [chexmix-yeast-example.tar.gz](http://lugh.bmb.psu.edu/software/chexmix/examples/chexmix-yeast-example.tar.gz)

Let’s demonstrate how ChExMix works in two different modes. Note that you will need to install [MEME](http://meme-suite.org/) to run the first example.

__Defining and assigning subtypes using both DNA motif discovery and ChIP-exo tag distributions:__

```{r, engine='sh', count_lines}
java -Xmx10G -jar chexmix.public.jar --geninfo sacCer3.info --threads 16 --expt Abf1_Reb1_merged_toppeaks.bam --ctrl masterNotag_20161101.bam --format BAM --back sacCer3.back --exclude sacCer3_exludes.txt --seq ~/group/genomes/sacCer3/ --memepath path_to_meme_bin/ --mememinw 6 --mememaxw 18 --out testchexmix_20180214_reb1_abf1_all > testchexmix_20180214_reb1_abf1_all.out 2>&1
```
Expected results can be found here: [testchexmix_20180214_reb1_abf1_all](http://lugh.bmb.psu.edu/software/chexmix/examples/testchexmix_20180214_reb1_abf1_all/ChExMix_testchexmix_20180214_reb1_abf1_all_results.html)

__Defining and assigning subtypes using only ChIP-exo tag distributions:__

```{r, engine='sh', count_lines}
java -Xmx10G -jar chexmix.public.jar --geninfo sacCer3.info --threads 16 --expt Abf1_Reb1_merged_toppeaks.bam --ctrl masterNotag_20161101.bam --format BAM --back sacCer3.back --exclude sacCer3_exludes.txt --nomotifs --out testchexmix_20180214_reb1_abf1_all_nomotifs > testchexmix_20180214_reb1_abf1_all_nomotifs.out 2>&1
```
Expected results can be found here: [testchexmix_20180214_reb1_abf1_all_nomotifs](http://lugh.bmb.psu.edu/software/chexmix/examples/testchexmix_20180214_reb1_abf1_all_nomotifs/ChExMix_testchexmix_20180214_reb1_abf1_all_nomotifs_results.html)

Note that due to some stochasticity between runs, your results may not match those above exactly, but should be broadly similar.

Output files
--------------
1. `OutName_results.html` is a html that you can open in a web browser. It summarizes the ChExMix run results including input data, binding event subtypes, and replicate information for binding events.

2. `OutName.events` is a tabular file which contains information about significant binding events. A header starts with # and contains the following information for conditions and replicates. 
  * `Name`: Condition or replicate name
  * `Index`: Index used for a condition or replicate
  * `SigCount`: Total number of tags for a condition or replicate
  * `SigCtrlScaling`: Factor used to scale signal and control experiments 
  * `SignalFraction`: Fraction of tags estimated to come from a foreground
  
Rest of the rows contains the following information:
  * `Point`: Genomic position of binding event in “chromosome:coordinates” format
  * `CondName_Sig`: Number of tags associated with binding event from signal experiments of a condition. Tags among replicates within a condition are combined. 
  * `CondName_Ctrl`: Number of tags associated with binding event from control experiments of a condition. Tags among replicates within a condition are combined.
  * `CondName_log2Fold`: Log2 fold differences of tag counts between signal and control experiments
  * `CondName_log2Q`: Log2 q-values for binding events 
  * `SubtypePoint`: Genomic position and strand of dominant subtype (subtype associated with the binding events with the highest probability) 
  * `Tau`: Probability of binding event associated with dominant subtype 
  * `SubtypeName`: Name of subtype
  * `SubtypeSequence`: Sequence associated with subtype at binding event. Sequences are only reported if subtypes are associated with motifs.
  * `SubtypeMotifSocre`: Log-likelihood motif score of a subtype at binding event. Scoring uses a 2nd-order Markov model based on whole genome nucleotide frequencies.

3. `OutName.bed` is a bed format file which contains binding event locations in a single base pair and q-value. You can load it to the UCSC genome browser. This uses [UCSC ENCODE narrowPeak format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12).
  * 1st: Chromosome
  * 2nd: Start
  * 3rd: End
  * 5th: Integer score for display. It is calculated as int(-10*log2 q-value). The value is saturated at 1000. 

4. `OutName.all.events.table` is a tabular file which contains information about potential binding events and uses a similar format to OutName.events file. This file contains the following additional information:
  * `ActiveConds`: Indices of conditions that this event is still active (having non-zero probability)

5. `OutName.replicates.counts` is a tabular file which contains information about replicate tag counts and p- and q-values associated with potential binding events
  * `Point`: Genomic position of binding event in “chromosome:coordinates” format
  * `CondName:RepName`: Number of tags associated with binding event
  * `CondName:RepName_log2P`: Log2 p-value 
  * `CondName:RepName_log2Q`: Log2 q-value 

6. `OutName_RepName.repevents.txt` is a tabular file which contains information about replicate tag counts and q-values associated with significant binding events. The format of this file is similar to Name.replicates.counts. 

7. `OutName.all.replicationcodes.table` is a tabular file which contains information about consistency of replicates. A header describes meaning of labels.
  * `BindingEvent`: Genomic position of binding event in “chromosome:coordinates” format
  * `CondName`: Label of replicate consistency information explained in header. 

Contact
--------------
For queries, please contact Naomi (nuy11@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  
Version 0.4 (2019-04-25): Updates on HTML output, heatmap and sequence plots. Create extra outputs in event file. Updates on replicate handling options. Bug fixes.

Version 0.3 (2019-01-15): Introducing --lenient and --lenientplus modes for determining alternate ways of defining final binding events when replicated experiments are present. Output files now consist of only one ".events" file (to reduce confusion), files listing binding events found in each replicate, and a file reporting on the replication status of each binding event. Fixed several bugs, including relating to the final heatmap and composite plots displayed on the output HTML page, and a bug in --noclustering. 

Version 0.2 (2018-10-03): Updates to support Galaxy integration, added command line options to set 1) user provided motifs in seeding subtypes and 2) a Markov background threshold for motif detection, updates on EM training plots and html output, and bug fixes.

Version 0.1 (2018-02-14): Initial release.
