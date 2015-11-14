#!/bin/bash
#Author Benjamin Menkuec
#License: LGPL

width=$1
reffile=$2
peakfile=$3
selected_motif=$4
motif_file=$5
motif_threshold=$6

#extend width of peaks
less $peakfile | awk -v width=$width 'BEGIN {OFS="\t"} $1!="track" \
	{\
    start = $2-width;\
    if (start < 0) \
        start = 0; \    
    print $1,start,$3+width,$4}' > temp.bed

#get sequence for peaks
bedtools getfasta -fi $reffile -bed temp.bed -fo temp.fa

#convert repeats to NNNNNN
less temp.fa | awk 'BEGIN {OFS="\t"} \
	{if(substr($1,1,1)==">"){print $0} \
	else{gsub(/[a-z]/,"N",$0);print $0}}'\
	> temp_no_repeats.fa

#find top N motifs until the selected one, maximum motif size = 10
if [ "$motif_file" == "" ] ; then
    dreme -m $selected_motif -maxk 10 -p temp_no_repeats.fa
    motif_file=./dreme_out/dreme.xml
fi
motif=$selected_motif

echo $motif > current_motif.txt
echo -------------------------------------------------
echo using motif file $motif_file with motif $motif for $peakfile
echo -------------------------------------------------

#find top motif in peak sequences
#set threshold so that only exact matches will be reported
fimo --qv-thresh --thresh $motif_threshold --max-stored-scores 20000000 --motif $motif $motif_file temp_no_repeats.fa

#generate centered motif locations
#	offset = int(($4+$5)/2); \
less ./fimo_out/fimo.gff | grep "^chr"| awk 'BEGIN{OFS="\t"} \
	{split($1,chr,":"); split(chr[2],loc,"-"); \
    offset = int(($4+$5)/2); \
    if ($7=="+") \
	print chr[1],loc[1]+offset,loc[1]+offset+1,$9}' \
	> centered_motif_locations.bed

#delete temporary files
#rm temp.bed
#rm temp.fa
#rm temp_no_repeats.fa
