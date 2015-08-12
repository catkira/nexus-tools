echo off
cls
echo source files:  %1 &
REM del %1_out.fa
seqan_flexbar.exe %1 -er 0.2 -ol 4 -fm 22 -ml 18 -t -tt -tl 5 -tnum 4 -ss -b P:\git\chip-nexus\data\barcodes.fa -a P:\git\chip-nexus\data\adapter.fa -o %1_out.fastq
REM P:\git\chip-nexus\data\sra\SRR1175698.fastq.gz %2 -fm 22 -ml 18 -t -tn 5 -app -fr 10000 -tnum 1 -b P:\git\chip-nexus\data\barcodes.fa -a P:\git\chip-nexus\data\adapter.fa -ol 4 -er 0.2 -o P:\result.fastq.gz
pause