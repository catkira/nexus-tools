#!/bin/bash
#Author Benjamin Menkuec
#License: LGPL

width=$1
reffile=$2
peakfile=$3
selected_motif=$4

#extend width of peaks
less $peakfile | awk -v width=$width 'BEGIN {OFS="\t"} $1!="track" \
	{print $1,$2-width,$3+width,$4}' > temp.bed

#get sequence for peaks
bedtools getfasta -fi $reffile -bed temp.bed -fo temp.fa

#convert repeats to NNNNNN
less temp.fa | awk 'BEGIN {OFS="\t"} \
	{if(substr($1,1,1)==">"){print $0} \
	else{gsub(/[a-z]/,"N",$0);print $0}}'\
	> temp_no_repeats.fa

#find top N motifs until the selected one, maximum motif size = 10
dreme -m $selected_motif -maxk 10 -p temp_no_repeats.fa

#select m01 for first motif, m02 for second ....
motif_idline=$(grep "m0$selected_motif" ./dreme_out/dreme.xml)
motif=$(awk 'BEGIN{FS="\""}{print $4}' <<< $motif_idline)
echo "selected motif: " $motif

#find top motif in peak sequences
fimo --motif $motif ./dreme_out/dreme.txt temp_no_repeats.fa

#generate centered motif locations
less ./fimo_out/fimo.gff | grep "^chr"| awk 'BEGIN{OFS="\t"} \
	{split($1,chr,":"); split(chr[2],loc,"-"); \
	offset = int(($4+$5)/2); \
	print chr[1],loc[1]+offset,loc[1]+offset+1,$9}' \
	> centered_motif_locations.bed

#delete temporary files
rm temp.bed
rm temp.fa
rm temp_no_repeats.fa
