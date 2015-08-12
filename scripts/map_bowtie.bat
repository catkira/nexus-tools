REM 1st parameter: Preprocessed FASTQ
REM 2nd parameter: Path to bowtie reference genome
cls
python P:\bowtie-1.1.2\bowtie -S -p 4 --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata %2 %1 %1.sam
