#!/usr/bin/env python3

import subprocess
import tempfile

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    """
    This function will run blastn-short on the files and isolate pident, length, and qlen (defined below.
    Then, it will return "hits" which correspond to 80% accuracy and same length.
    """
    
    #https://stackoverflow.com/questions/26445313/blast-using-python-subprocess-call-about-alignment-def-and-no-idea-what-is-wrong
    blast = subprocess.run(['blastn', '-query', primer_file, '-subject', assembly_file, '-task', 'blastn-short', '-outfmt', '6 std qlen'], capture_output=True, text=True)
    lines = blast.stdout.strip().split("\n")

    hits = []
    for line in lines:
        parts = line.split("\t")
        pident = float(parts[2])
        length = int(parts[3])
        qlen = int(parts[12])
        if pident >= 80 and length == qlen:
            hits.append(parts)
    
    hits.sort(key=lambda x: int(x[8]))
    return(hits)

def step_two(sorted_hits: list[str], max_amplicon_size: int) -> list[tuple[list[str]]]:   
    """
    This function will first compares hits that are next to each other within the sorted list. 
    It looks for hits that are on the same sseqid, face towards each other, and have the 3' end within the given range. 
    """
    
    pairs = []
    num_hit = len(sorted_hits) - 1
    for hit in range(num_hit):
       
        #comparing the ones directly next to each other
        first = sorted_hits[hit]
        second = sorted_hits[hit + 1]

        #start and end for first and second (4 lines)
        first_start = int(first[8])
        first_end = int(first[9])
        second_start = int(second[8])
        second_end = int(second[9])

        #is true if either are first or second are the forward strand
        first_forward = first_start < first_end
        second_forward = second_start < second_end
        
        #skip if facing same way (as per document)
        if first_forward == second_forward:
            continue

        #obtain 3 prime ends to compare for distance
        if first_forward:
            first_3 = first_end
        else:
            first_3 = first_start
        
        if second_forward:
            second_3 = second_end
        else:
            second_3 = second_start
        
        #compare to distance
        if abs(first_3 - second_3) <= max_amplicon_size:
            pairs.append((first,second))
    
    return(pairs)

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    """
    This function will use seqtk to extract the DNA sequences of each primer hit.
    It will return all sequencces in FASTA format. 
    """
    
    amplicon = []
    for first, second in hit_pairs:
       
        #same as before, get start and ends
        first_start = int(first[8])
        first_end = int(first[9])
        second_start = int(second[8])
        second_end = int(second[9])   

        #same as before, T based on whats first
        first_forward = first_start < first_end
        second_forward = second_start < second_end

        #find the inner boundaries (remove primer)
        if first_forward:
            first_3 = first_end + 1  
        else:
            first_3 = first_end - 1  

        if second_forward:
            second_3 = second_end + 1   
        else:
            second_3 = second_end - 1  

        #obtain start and end (start is the smallest out of the ends, end is the biggest)
        start = min(first_3, second_3)
        end = max(first_3, second_3)
        chrom = first[1] #same for first and second as they are the same contig

        #make list per format in https://github.com/lh3/seqtk 
        amplicon.append(f"{chrom}\t{start-1}\t{end}\n")

    #per document, seqtk needs a file
    #https://stackoverflow.com/questions/39983886/python-writing-and-reading-from-a-temporary-file 
    with open("tmp.text", "w") as f:
        f.writelines(amplicon)
        
    #run seqtk
    seqtk = subprocess.run(["seqtk", "subseq", assembly_file, "tmp.text"], capture_output=True, text=True)
    return(seqtk.stdout)


