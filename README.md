
 ### Antimicrobial resistance genes pipeline
The pipeline indicated below is to be used in identifying antimicrobial resistance genes using long reads sequenced from Oxford nanopore.
The sequenced samples were previously identified as being resistant to commonly used antibiotics using antibiotics susceptibility tests (ASTs). 

The tools used in the pipeline are listed below and used after guppy was used to basecall and barcode the reads.

### Quality control
 
The quality of the nanopore reads were first assessed using FastQC. 



Fastqc    - Check the quality of the reads.

Porechop  -  removal of adapters.

Filtlong - Filtering of reads. 

Flye     - de novo assembly of the reads.

Quast - Checking the quality of the assembled genome.

Bwa  - indexing the reads and mapping reads to assembled genome.

racon + medaka - polishing of the reads. 

Abricate - identify AMR genes and plasmid



