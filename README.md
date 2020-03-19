# Spacer-Seq
MiSeq data ingestion and processing for targeted CRISPR spacer acquisition experiments.

Necessary Python packages:
1. fuzzysearch
2. Biopython
3. numpy
4. xlsxwriter
5. xlrd

## unzip_illumina.py
Takes in raw files from Illumina run and converts into FASTQ files.

## trimmomatic-0.39.jar
Java executable that takes in FASTQ files and trims length based on quality score of bases.

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

http://www.usadellab.org/cms/?page=trimmomatic

## trimFASTQs.py
Runs trimmomatic on all unzipped, indexed FASTQs

## Fed_SPCRs_v2.py
Takes in trimmed FASTQ files and identifies newly acquired spacers. Generates several FASTQ files, filled with sequences of newly acquired spacers, grouped by shared attributes.

## SPCR_blast_v2.py
Takes in trimmed FASTQ files and performs local BLAST against sequence data of interest (eg. plasmid map, host cell genome, etc.). Returns FASTQ files filled with new spacer sequences, based on alignment.

## oComp_Ordering_SBK.py
Takes in trimmed FASTQ files and three different protospacer sequences, and returns an Excel file with data regarding position of specified protospacers in sequenced CRISPR arrays.

## SPCR_hunt.py
