#!/usr/bin/env python3

#usage should be get_homolog_seqs.py <blast file> <bed file> <assembly file> <output file>

import sys
def homologs(blastfile: str, bedfile: str, assemblyfile: str, outputfile: str) -> None:
    """
    This script will identify homologous genes from a blast result.
    It will then match then to genes in a BED file.
    For minus strands, it will reverse and obtain the compliment.
    Finally, it will write the homologous genes to an output FASTA file. 
    """
    # read blast file
    hits = []
    with open(blastfile) as fin:
        for line in fin: # .readlines() is the default iter method for the open file class
            
            # unpack and convert types of desired columns. This is ugly. We'll revisit later...
            _, sid, pcnt, matchlen, _, _, _, _, sstart, send, _, _, qlen = line.split()
            pcnt = float(pcnt)
            matchlen = int(matchlen)
            sstart = int(sstart)
            send = int(send)
            qlen = int(qlen)
        
            # Keep hits that could be homologs
            if pcnt > 30 and matchlen > 0.9*qlen:
                # We could store matches as a list or tuple.
                # We won't want to modify the elements so a tuple is "safer" in that we then can't modify it by mistake
                hits.append((sid, sstart, send))

    # Now read the bed file
    feats = []
    with open(bedfile) as fin:
        for line in fin:
            bed_sid, bed_start, bed_end, gene, *_ = line.split() # an asterisk before a variable name when unpacking makes that variable store remaining elements
            bed_start = int(bed_start)
            bed_end = int(bed_end)
            
            feats.append((bed_sid, bed_start, bed_end, gene))

    # Now we have our two datasets read in, we can loop over them to find matches
    homologs = []
    for blast_sid, blast_sstart, blast_send in hits: # unpack our blast data
        for bed_sid, bed_start, bed_end, gene in feats:
            # Don't bother checking the rest if the sid doens't match
            if blast_sid != bed_sid:
                continue
            
            # Once we are dealing with features at higher index locations than our hit, go to the next hit (break loop over feats)
            if blast_sstart <= bed_start or blast_send <= bed_start:
                break
            
            # Otherwise, check if the hit is inside the feature
            if (blast_sstart > bed_start
                and blast_sstart <= bed_end
                and blast_send > bed_start
                and blast_send <= bed_end
            ):
                homologs.append(gene)
                break # Each BLAST hit will only be in one feature so move to next hit once you've found it

    # Get the unique homologs using a set()
    unique_homologs = set(homologs)

    #3. Extract the DNA sequence of identified homologous genes from the assembly sequence
    assembly_name = []
    assembly_seq = []

    with open(assemblyfile) as fin:
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                assembly_name.append(line[1:].split()[0])
                assembly_seq.append("")
            else:
                assembly_seq[-1] = assembly_seq[-1] + line
    #both should be the same length. first string will have names. the second will have dna sequences. indexes should match up.

    #4. If the gene is encoded on the - strand (indicated in the BED file), reverse complement the sequence
    reverse_gene = []
    reverse_seq = []
    bed_total = []

    with open(bedfile) as fin:
        for line in fin:
            line = line.strip()
            columns = line.split()

            bed_name = columns[0]
            bed_start = int(columns[1])
            bed_end = int(columns[2])
            bed_gene = columns[3]
            bed_strand = columns[5]

            bed_total.append((bed_gene, bed_name, bed_start, bed_end, bed_strand))

    #look through genes if its only in the unique homologs from sample script
    for gene in unique_homologs:
        match = None
        for g in bed_total:
            if g[0] == gene:
                match = g
                break
        if match is None:
            continue          
        bed_gene, bed_name, bed_start, bed_end, bed_strand = match 
            
        #check for assembly to bed matches
        if bed_name not in assembly_name:
            continue
            
        index = assembly_name.index(bed_name)
        seq = assembly_seq[index][bed_start:bed_end]
        
        #reverse and compliment if negative, leave alone if positive
        if bed_strand == "-":
            compliment = ""
            for b in seq:
                if b == "G":
                    compliment = compliment + "C"
                elif b == "C":
                    compliment = compliment + "G"
                elif b == "A":
                    compliment = compliment + "T"
                elif b == "T":
                    compliment = compliment + "A"
            reverse = compliment[::-1] #https://stackoverflow.com/questions/931092/how-do-i-reverse-a-string-in-python
            reverse_gene.append(bed_gene)
            reverse_seq.append(reverse)
        else:
            reverse_gene.append(bed_gene)
            reverse_seq.append(seq)


    #5. Write the sequences of the homologous genes to the specified output file (specified using command line input) in FASTA format (gene name as header).
    with open(outputfile, "w") as fout:
        for i in range(len(reverse_gene)):
            fout.write(">" + reverse_gene[i]+ "\n") 
            fout.write(reverse_seq[i] + "\n")

homologs(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])