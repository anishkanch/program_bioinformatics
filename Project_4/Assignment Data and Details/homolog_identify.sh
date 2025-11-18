#!/bin/bash

#usage of your script should be: 
#homolog_identify.sh <query.faa> <subject.faa> <bedfile> <outfile>

#load in file from usage of script
query_faa="$1"
subject_faa="$2"
bedfile="$3"
outfile="$4"

#blastn, $1 will be pident, $2 will be length, $3 will be qlen
#we want $1 (pident) to be 100 and for $2(length) to equal $3(qlen), did some research on how to use two conditions with awk via stackoverflow: Awk AND operator
tblastn -query "$query_faa" -subject "$subject_faa" -outfmt "6 sseqid pident length qlen sstart send" | awk '$2 > 30 && $3>($4*0.9)' > "tmp.txt"

#conditionals to determine which genes in a BED file contain the matches listed tmp.txt
while read sseqid pident alength qlen sstart send; do

    if ((sstart < send)); then
        start=$sstart
        end=$send
    elif ((sstart > send)); then
        start=$send
        end=$sstart
    fi
    
    while read chrom gene_start gene_end gene_name score strand; do
        if [[ "$sseqid" == "$chrom" && $start -ge $gene_start && $end -le $gene_end ]]; then
            echo "$gene_name"
        fi
    done < "$bedfile"

#outfile will contain the genes in BED that are homologs
done < "tmp.txt" | sort -u > "$outfile"


