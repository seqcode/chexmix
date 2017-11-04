# ChExMix
ChExMix: the ChIP-exo mixture model

ChExMix aims to characterize protein-DNA binding subtypes in ChIP-exo experiments. ChExMix assumes that different regulatory complexes will result in different protein-DNA crosslinking signatures in ChIP-exo data, and thus analysis of ChIP-exo sequencing tag patterns should enable detection of multiple protein-DNA binding modes for a given regulatory protein. ChExMix uses a mixture modeling framework to probabilistically model the genomic locations and subtype membership of protein-DNA binding events, leveraging both ChIP-exo tag enrichment patterns and DNA sequence information. In doing so, ChExMix offers a more principled and robust approach to characterizing binding subtypes than simply clustering binding events using motif information.
