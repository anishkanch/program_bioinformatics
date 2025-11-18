#!/usr/bin/env python3

from ispcr import step_one

sorted_good_hits = step_one(
	primer_file="data/general_16S_515f_806r.fna",
	assembly_file="data/Vibrio_cholerae_N16961.fna"
)

for hit in sorted_good_hits:
	print(hit)
