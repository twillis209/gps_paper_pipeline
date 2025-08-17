#!/bin/bash
wget -O md_ex.tsv.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003219/harmonised/26691988-GCST003219-EFO_0001365.h.tsv.gz
# Smaller study is Europeans only
#wget -O glaucoma_ex.tsv.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009722/harmonised/31959993-GCST009722-EFO_0000516.h.tsv.gz
wget -O glaucoma_ex.tsv.gz http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011770/GCST90011770_buildGRCh37.tsv.gz
