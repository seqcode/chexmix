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

Options (Required/important options are in __bold__.)

1. General:

  * --__out__ \<prefix>: Output file prefix. All output will be put into a directory with the prefix name. 
  * --threads \<n\>:  Use n threads during binding event detection. Default is 1 thread.
  * --verbose: Flag to print intermediate files and extra output
  * --memepath \<path\>: path to the meme bin dir (default: meme is in $PATH).

2. Specifying the Genome:

  * --__geninfo__ \<genome info file\>:  This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ \<path\> : A directory containing fasta format files corresponding to every named chromosome is required if you want to run motif-finding or use a motif-prior within ChExMix.

3. Loading Data:

  * --__exptCONDNAME-REPNAME__ \<file\>: Defines a file containing reads from a signal experiment. Replace CONDNAME and REPNAME with appropriate condition and replicate labels.
  * --__ctrlCONDNAME-REPNAME__ \<file\>: Optional arguments. Defines a file containing reads from a control experiment. Replace CONDNAME and REPNAME with appropriate labels to match a signal experiment (i.e. to tell ChExMix which condition/replicate this is a control for). If you leave out a REPNAME, this file will be used as a control for all replicates of CONDNAME.  
  * --__format__ \<SAM/BED/IDX\>: Format of data files. All files must be the same format if specifying experiments on the command line. Supported formats are SAM/BAM, BED, and IDX index files.
 
Instead of using the above options to specify each and every ChIP-seq data file on the command-line, you can instead use a design file:
 
  * --__design__ \<file\>: A file that specifies the data files and their condition/replicate relationships. See [here](http://lugh.bmb.psu.edu/software/multigps/example.design) for an example design file. The file should be formatted to contain the following pieces of information for each data file, in this order and tab-separated:
  
    * File name
    * Label stating if this experiment is “signal” or “control”
    * File format (SAM/BAM/BED/IDX) – mixtures of formats are allowed in design files
    * Condition name
    * Replicate name (optional for control experiments – if used, the control will only be used for the corresponding named signal replicate)
 
4. Other ChExMix options (Recommend using defaul options):

  * --mememinw \<value\>: minw arg for MEME. Default=6.
  * --mememaxw \<value\>: maxw arg for MEME. Default=13. This value should always be less than "maxScanLen".
  * --memenmotifs \<int\>: Number of motifs MEME should find in each condition (default=3)
  * --memeargs \<args\> : Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)

Example
--------------
This example runs ChExMix v0.1.2 on simulated dataset. Simulated data to run this example can be found [here]

Command:
```{r, engine='sh', count_lines}
java -Xmx20G -jar chexmix.jar --threads 10 --expt example.bam --BAM --memepath path-to-meme --geninfo sacCer3.info --seq path-to-genomes/sacCer3/ --verbose  --out example
```

Results can be found [here]

Contact
--------------
For queries, please contact Naomi (nuy11@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  
Version 0.1 (2017-11-06): Initial release.
