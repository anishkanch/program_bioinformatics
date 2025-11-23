#!/usr/bin/env python3

import subprocess
import argparse
from magnumopus.sam import SAM

#command line arguments - usage of your script should be map_consensus.py -1 <READ1> -2 <READ2> -r <REF_SEQS> [-s <SEQ_NAME>] 
cl = argparse.ArgumentParser("Map the provided reads to the provided reference sequence")
cl.add_argument("-1", help = "Read 1", required=True, dest="read_file")
cl.add_argument("-2", help = "Read 2", required=True, dest="second_read_file")
cl.add_argument("-r", help = "Reference File", required=True, dest="reference_file")
#If the user specified a sequence name in the reference, print the consensus of reads that mapped to the named sequence
cl.add_argument("-s", help = "Sequence Name", required=False, dest="seq_name")

arguments = cl.parse_args()

#mapping with minimap2
minimap = subprocess.run(["minimap2", "-ax", "sr","-B", "1", "--score-N", "0", "-k", "10", arguments.reference_file, arguments.read_file, arguments.second_read_file], capture_output=True, text=True)

#write sam output
with open("sam_file", "w") as file:
      file.write(minimap.stdout)

#load obj
sam_object = SAM.from_sam("sam_file")

#If the user did not specify a sequence name, print the consensus of reads that mapped to the reference with the best mapping
if arguments.seq_name:
      consensus = sam_object.consensus(arguments.seq_name)
      header = ">"+arguments.seq_name+"_consensus"

else:
      consensus = sam_object.best_consensus()
      header = ">best_consensus"

#return
print(header)
print(consensus)