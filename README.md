
 ### Antimicrobial resistance genes pipeline
The pipeline indicated below is to be used in identifying antimicrobial resistance genes using long reads sequenced from Oxford nanopore.
The sequenced samples were previously identified as being resistant to commonly used antibiotics using antibiotics susceptibility tests (ASTs). 

The tools used in the pipeline are listed below and used agter guppy was used to basecall and barcode the reads.

Fastqc    - Check the quality of the reads.

Porechop  -  removal of adapters.

Filtlong - Filtering of reads. 

Flye     - de novo assembly of the reads.

racon + medaka - polishing of the reads. 

