#!/bin/bash

width=$1
reffile=$2
peakfile=$3
outfile=$4

#extend width of peaks
less $peakfile | awk -v width=$width 'BEGIN {OFS="\t"} $1!="track" {print $1,$2-width,$3+width,$4}' > temp.bed

#get sequence for peaks
bedtools getfasta -fi $reffile -bed temp.bed -fo temp.fa

#less temp.fa | awk 'BEGIN {OFS="\t"} {if(substr($1,1,1)==">"){print $0} else{gsub(/[a-z]/,"N",$0); print $0}}'

#convert repeats to NNNNNN
less temp.fa | awk 'BEGIN {OFS="\t"} {if(substr($1,1,1)==">"){print $0} else{gsub(/[a-z]/,"N",$0);print $0}}' > $outfile

#find top10 motifs
dreme -m 10 -maxk 10 -p $outfile

#find top motif in peak sequences
fimo --motif 1 ./dreme_out/dreme.txt $reffile

#intersect with peaks
#bedtools intersect -u $peakfile -b
