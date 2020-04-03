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

Instructions: Run with python3. Requires one line argument. Provide path of "FASTQ_Generation..." folder downloaded from MiSeq run as line argument. Output files are generated in a folder titled "rawIndexedFASTQs".

Example:
```
python3 unzip_illumina.py ~/MiSeq_Data_Dev/msSBK_4_reIndex/FASTQ_Generation_2020-02-14_17_03_25Z-207492285
```

## trimmomatic-0.39.jar
Java executable that takes in FASTQ files and trims length based on quality score of bases.

Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

http://www.usadellab.org/cms/?page=trimmomatic

## trimFASTQs.py
Runs trimmomatic on all unzipped, indexed FASTQs.

Instructions: Run with python3. Requires one line argument. Provide path of folder containing unzipped FASTQs aka. "rawIndexedFASTQs". Output files are generated in a folder titled "trimmedFASTQs".

Example:
```
python3 trimFASTQs.py ~/MiSeq_Data_Dev/msSBK_4_reIndex/rawIndexedFASTQs
```

## Fed_SPCRs_v2.py
Takes in trimmed FASTQ files corresponding to a single index and identifies newly acquired spacers. Generates several FASTQ files, filled with sequences of newly acquired spacers, grouped by shared attributes.

## SPCR_blast_v2.py
Takes in trimmed FASTQ files and performs local BLAST against sequence data of interest (eg. plasmid map, host cell genome, etc.). Returns FASTQ files filled with new spacer sequences, based on alignment. This step is not necessary if you are searching for a known spacer. It can be helpful when trying to characterize  new spacers of unknown origin.

## oComp_Ordering_SBK.py
Takes in trimmed FASTQ files corresponding to a single index and three different protospacer sequences, and returns an Excel file with data regarding position of specified protospacers in sequenced CRISPR arrays.

## analyzeTrimmedIndices.py
Runs Fed_SPCRs_v2.py and/or oComp_Ordering_SBK.py on all trimmed multi-FASTQs, in a given folder, generated by trimFASTQs.py.

Instructions: Run with python3. Requires three line arguments.
* Line argument one: path of folder containing trimmed FASTQs aka. "trimmedFASTQs".
* Line argument two: Instruction to extract and/or tabulate order of spacers. `e` = extract; `o` = order; `v` = analyze variable regions. You can provide combinations of these arguments to do multiples things. Ex: `eov` = extract + order + analyze variable regions.
* Line argument three: path and filename of xlsx (Excel) ordering config file containing spacers of interest for ordering function. If not doing this analysis, type `None`.
* Line argument four: path and filename of xlsx (Excel) ordering config file containing variable region options for analysis. If not doing this analysis, type `None`.

Example:
```
python3 analyzeTrimmedIndices.py ~/MiSeq_Data_Dev/msSBK_4_reIndex/trimmedFASTQs ~/MiSeq_Data_Dev/msSBK_4_reIndex/msSBK_4_ordering_config.xlsx ord
```
