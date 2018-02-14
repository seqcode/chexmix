# ChExMix
ChExMix: the ChIP-exo mixture model

ChExMix aims to characterize protein-DNA binding subtypes in ChIP-exo experiments. ChExMix assumes that different regulatory complexes will result in different protein-DNA crosslinking signatures in ChIP-exo data, and thus analysis of ChIP-exo sequencing tag patterns should enable detection of multiple protein-DNA binding modes for a given regulatory protein. ChExMix uses a mixture modeling framework to probabilistically model the genomic locations and subtype membership of protein-DNA binding events, leveraging both ChIP-exo tag enrichment patterns and DNA sequence information. In doing so, ChExMix offers a more principled and robust approach to characterizing binding subtypes than simply clustering binding events using motif information.

Downloading Executables
--------------
Executables will be available soon from: http://mahonylab.org/software/chexmix 
Check back soon!

Building from Source
--------------
If you want to build the code yourself, you will need to first download and build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

Dependencies:
--------------
1. ChExMix requires Java 8+. 
2. ChExMix depends on [MEME](http://meme-suite.org/) being available in $PATH (tested with MEME version 4.11.3).

Note:
--------------
ChExMix loads all data to memory, so you will need a lot of available memory if you are running analysis over many conditions or large datasets.

Citation:
--------------

Running ChExMix
--------------
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
  * --__seq__ \<path\> : A directory containing fasta format files corresponding to every named chromosome is required if you want to find subtypes run motif-finding or use a motif-prior within ChExMix.
  * --__back__ \<path\> : A file containing Markov background model for the genome is required if you want to run motif-finding or use a motif-prior within ChExMix.

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
    
 Limits on how many reads can have their 5′ end at the same position in each replicate:

 * --fixedpb \<value\>: Fixed per-base limit.
 * --poissongausspb \<value\>: Filter per base using a Poisson threshold parameterized by a local Gaussian sliding window (i.e. look at neighboring positions to decide what the per-base limit should be).
 * Default behavior is to estimate a global per-base limit from a Poisson distribution parameterized by the number of reads divided by the number of mappable bases in the genome. The per-base limit is set as the count corresponding to the 10^-7 probability level from the Poisson.
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
  * --minmodelupdateevents \<int\>: Minimum number of events to support an update (default=100)
  * --prlogconf \<value\>: Poisson log threshold for potential region scanning (default=-6)
  * --fixedalpha \<int\>: Impose this alpha (default: set automatically). The alpha parameter is a sparse prior on binding events in the ChExMix model. It can be interpreted as a minimum number of reads that each binding event must be responsible for in the model. 
  * --alphascale \<value\>: Alpha scaling factor (default=1.0). Increasing this parameter results in stricter binding event calls.
  * --betascale \<value\>: Beta scaling factor (default=0.05). The beta parameter is a sparse prior on binding event subtype assignment in the ChExMix model. Increasing this parameter may result in inaccurate subtype assignment.
  * --epsilonscale \<value\>: Epsilon scaling factor (default=0.2). The epsilon parameter control a balance between motif and read distribution in subtype assignment. Increasing this parameter will increase the contribution of motif and decreases the contribution of read distribution.
  * --mlconfignotshared: Flag to not share component configs in the ML step
  * --exclude \<file\>: File of regions to ignore
  * --peakf \<file\>: File of peaks to initialize component positions

__Finding ChExMix subtypes__:

1. Using motif:

  * --__memepath__  \<path\>: Path to the meme bin dir (default: meme is in $PATH). MEME path is required for motif finding.
  * --nomotifs \<value\>: Flag to turn off motif-finding & motif priors
  * --nomotifprior \<value\>: Flag to turn off motif priors only
  * --mememinw \<value\>: minw arg for MEME (default=6).
  * --mememaxw \<value\>: maxw arg for MEME (default=18).
  * --memenmotifs \<int\>: Number of motifs MEME should find in each condition (default=3)
  * --memeargs \<args\>: Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)
  * --minroc \<value\>: Minimum motif ROC value (default=0.7)
  * --minmodelupdaterefs \<int\>: Minimum number of motif reference to support an subtype distribution update (default=50)
 
 2. Using read distribution:

  * --noclustering: Flag to turn off read distribution clustering
  * --pref \<value\>: Preference value for read distribution clustering (default=-0.1)
  * --numcomps \<int\>: Number of components to cluster (default=500)
  * --win \<int\>: Read profile window size (default=150)
  
__Reporting binding events__:

  * --q \<value\>: Q-value minimum (default=0.01)
  * --minfold \<value\>: Minimum event fold-change vs scaled control (default=1.5)

Example
--------------
This example runs ChExMix v0.1 on simulated dataset. Simulated data to run this example can be found [here]

Command:
```{r, engine='sh', count_lines}
java -Xmx20G -jar chexmix.jar --threads 10 --expt example.bam --format BAM --memepath path-to-meme --geninfo sacCer3.info --seq path-to-genomes/sacCer3/ --verbose  --out example
```

Results can be found [here]

Contact
--------------
For queries, please contact Naomi (nuy11@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  
Version 0.1 (2018-02-14): Initial release.
