#!/bin/bash

export OMP_NUM_THREADS=1

for /L %%N in (1,1,9) do (
	echo "msSBK_2_0%%N"
	java -jar ./Trimmomatic-0.39.jar SE msSBK_2/msSBK_2_0%%N.fastq msSBK_2/msSBK_2_0%%N_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:6:25 MINLEN:60 )

for /L %%N in (10,1,42) do (
	echo "msSBK_2_%%N"
	java -jar ./Trimmomatic-0.39.jar SE msSBK_2/msSBK_2_%%N.fastq msSBK_2/msSBK_2_%%N_trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:6:25 MINLEN:60 )

cmd /k
#After creating bsub.sh, need to fix windows characters w/ "dos2unix bsub.sh", then make it executable w/ "chmod +x bsub.sh"
#run with ./bsub.sh
