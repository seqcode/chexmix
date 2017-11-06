# ChExMix
ChExMix: the ChIP-exo mixture model

ChExMix aims to characterize protein-DNA binding subtypes in ChIP-exo experiments. ChExMix assumes that different regulatory complexes will result in different protein-DNA crosslinking signatures in ChIP-exo data, and thus analysis of ChIP-exo sequencing tag patterns should enable detection of multiple protein-DNA binding modes for a given regulatory protein. ChExMix uses a mixture modeling framework to probabilistically model the genomic locations and subtype membership of protein-DNA binding events, leveraging both ChIP-exo tag enrichment patterns and DNA sequence information. In doing so, ChExMix offers a more principled and robust approach to characterizing binding subtypes than simply clustering binding events using motif information.

Downloading Executables
--------------

Building from Source
--------------
If you want to build the code yourself, you will need to first download and build the seqcode-core library (https://github.com/seqcode/seqcode-core) and add its build/classes and lib directories to your CLASSPATH.

Dependencies:
--------------
1. ChExMix requires Java 8+. 
2. ChExMix depends on [MEME](http://meme-suite.org/) (tested with MEME version 4.11.3).

Citation:
--------------

Running ChExMix
--------------
Running from a jar file:

```{r, engine='sh', count_lines}
java -Xmx20G -jar ChExMix.jar <options - see below>
```

In the above, the “-Xmx20G” argument tells java to use up to 20GB of memory. If you have installed source code from github, and if all classes are in your CLASSPATH, you can run ChExMix as follows:

```{r, engine='sh', count_lines}
java -Xmx20G org.seqcode.projects.chexmix.ChExMix <options - see below>
```

Options (Required/important options are in __bold__.)

1. General:

  * --__out__ \<prefix>: Output file prefix. All output will be put into a directory with the prefix name. 
  * --threads \<n\>:  Use n threads to train SeqUnwinder model. Default is 5 threads.
  * --verbose: Flag to print intermediate files and extra output
  * --memepath \<path\>: path to the meme bin dir (default: meme is in $PATH).

2. Specifying the Genome:

  * --__geninfo__ \<genome info file\>:  This file should list the lengths of all chromosomes on separate lines using the format chrName\<tab\>chrLength. You can generate a suitable file from UCSC 2bit format genomes using the UCSC utility “twoBitInfo”. The chromosome names should be exactly the same as those used in your input list of genomic regions. 
   
      The genome info files for some UCSC genome versions:  
      | [hg18](http://lugh.bmb.psu.edu/software/multigps/support/hg18.info) | [hg19](http://lugh.bmb.psu.edu/software/multigps/support/hg19.info) | [hg38](http://lugh.bmb.psu.edu/software/multigps/support/hg38.info) | [mm8](http://lugh.bmb.psu.edu/software/multigps/support/mm8.info) | [mm9](http://lugh.bmb.psu.edu/software/multigps/support/mm9.info) | [mm10](http://lugh.bmb.psu.edu/software/multigps/support/mm10.info) | [rn4](http://lugh.bmb.psu.edu/software/multigps/support/rn4.info) | [rn5](http://lugh.bmb.psu.edu/software/multigps/support/rn5.info) | [danRer6](http://lugh.bmb.psu.edu/software/multigps/support/danRer6.info) | [ce10](http://lugh.bmb.psu.edu/software/multigps/support/ce10.info) | [dm3](http://lugh.bmb.psu.edu/software/multigps/support/dm3.info) | [sacCer2](http://lugh.bmb.psu.edu/software/multigps/support/sacCer2.info) | [sacCer3](http://lugh.bmb.psu.edu/software/multigps/support/sacCer3.info) |
  * --__seq__ \<path\> : A directory containing fasta format files corresponding to every named chromosome is required.

3. Input Genomic Regions:

  * --__genregs__ \<file\>: Genomic regions with annotations filename OR --__genseqs__\<file\>: Sequences with annotations filename. A tab delimited file of a list of genomic points/sequences and corresponding annotations/labels. A simple example :
      ```{r, engine='sh', count_lines}
	GenRegs file:
	chr10:100076604	enhancer;shared
	chr6:100316177	promoter;celltypeA

	GenSeqs file:
	ATTGC....TTA	enhancer;shared
	CGTAA....GGT	promoter;celltypeA
      ```
  * --win \<int\>:  Size of the genomic regions in bp. Default = 150.
  * --makerandregs: Flag to make random genomic regions as an extra outgroup class in classification (only applicable when genome is provided).
  
4. Other ChExMix options (Recommend using defaul options):

  * --minscanlen \<value\>: Minimum length of the window to scan *K*-mer models. Default=8.
  * --maxscanlen \<value\>: Maximum length of the window to scan *K*-mer models. Default=14.
  * --hillsthresh \<value\>: Scoring threshold to identify hills. Default=0.1.
  * --mememinw \<value\>: minw arg for MEME. Default=6.
  * --mememaxw \<value\>: maxw arg for MEME. Default=13. This value should always be less than "maxScanLen".
  * --memenmotifs \<int\>: Number of motifs MEME should find in each condition (default=3)
  * --memeargs \<args\> : Additional args for MEME (default:  -dna -mod zoops -revcomp -nostatus)
  * --memesearchwin \<value\>: Window around hills to search for discriminative motifs. Default=16
  * --a \<int\>: Maximum number of allowed ADMM iterations. Default=400.


Example
--------------
This example runs ChExMix v0.1.2 on simulated dataset. Simulated sequences to run this example can be found [here]

Command:
```{r, engine='sh', count_lines}
java -Xmx20G -jar chexmix.jar --out example --threads 10 --debug --memepath path-to-meme --geninfo mm10.info --seq path-to-genomes/mm10/ --genseqs simulateOverlap.seqs --win 150 --mink 4 --maxk 5 --r 10 --x 3 --maxscanlen 15
```

Results can be found [here](http://lugh.bmb.psu.edu/software/sequnwinder/example/SeqUnwinder_results.html)

Contact
--------------
For queries, please contact Naomi (nuy11@psu.edu) or Shaun Mahony (mahony@psu.edu).

Major History:
--------------  
Version 0.1 (2017-11-06): Initial release.
