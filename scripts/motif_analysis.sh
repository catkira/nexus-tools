#!/bin/bash
#Author Benjamin Menkuec
#License: LGPL

width=$1
reffile=$2
peakfile=$3
selected_motif=$4
motif_file=$5

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

re='^[0-9]+$'
if [[ $selected_motif =~ $re ]] ; then
    #find top N motifs until the selected one, maximum motif size = 10
    if [ "$motif_file" == "" ] ; then
        dreme -m $selected_motif -maxk 10 -p temp_no_repeats.fa
        motif_file=./dreme_out/dreme.xml
    fi

    #select m01 for first motif, m02 for second ....
    motif_idline=$(grep "m0$selected_motif" $motif_file)
    motif=$(awk 'BEGIN{FS="\""}{print $4}' <<< $motif_idline)
else
    motif=$selected_motif
    if [ "$motif_file" == "" ] ; then
        motif_file=./dreme_out/dreme.xml
    fi
fi

echo $motif > current_motif.txt
echo -------------------------------------------------
echo using motif file $motif_file with motif $motif
echo -------------------------------------------------

#find top motif in peak sequences
fimo --motif $motif $motif_file temp_no_repeats.fa

#generate centered motif locations
less ./fimo_out/fimo.gff | grep "^chr"| awk 'BEGIN{OFS="\t"} \
	{split($1,chr,":"); split(chr[2],loc,"-"); \
	offset = int(($4+$5)/2); \
	print chr[1],loc[1]+offset,loc[1]+offset+1,$9}' \
	> centered_motif_locations.bed

#delete temporary files
#rm temp.bed
#rm temp.fa
#rm temp_no_repeats.fa
