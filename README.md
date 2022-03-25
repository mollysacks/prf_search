#PRF Search

This repository contains scripts that can be used to search through DNA sequences for the features that stimulate -1 programmed ribosomal frameshifting (PRF). 

This pipeline was run on every annotated ORF in E. coli. If there was an intergenic region between one ORF and the next ORF, we added that region onto the 3' end of the first ORF. 

#Dependencies

Rscape
HMMER
ViennaRNA

Python3:

argparse
pandas
json
Bio.Seq
re
numpy
math
RNA (ViennaRNA Python wrapper)
matplotlib
shutil
mlines
datetime
random

#Run commands

To run the whole pipeline:

`bash search.sh -o E_coli_3p_utr -q mrna_plus_3utr.fa -d bacteria.1236.1.genomic.fna -r /path/to/your/installation/of/rscape_v1.6.1/bin/R-scape -p /path/to/this/repo/prf_search`

To score sequences without running HMMER or CaCoFold:

`python3 random_probability_eval.py -I All-genes-of-E.-coli-K-12-substr.-MG1655.fasta`

To run pipeline on 4000 randomly generated sequences following E. coli gene length distribution:

`python3 random_probability_eval.py -N 4000 -F All-genes-of-E.-coli-K-12-substr.-MG1655.fasta`


