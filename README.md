BIOL 7200 — Programming for Bioinformatics (Fall 2025)
------------------------------------------------------
This repository contains all programming assignments and projects completed for BIOL 7200, a course focused on computational problem-solving, sequence analysis, and biological data processing.
Assignments progress from foundational Bash scripting to full Python software modules capable of performing in-silico PCR, global alignment, SAM parsing, and consensus sequence construction.
The repository emphasizes clean organization, reproducibility, and portfolio-quality code.

Each Project_X folder contains the full code for that assignment (Bash or Python), plus any helper scripts or modules.
The Assignment Data and Details directory contains original instructions and raw data needed to test the code.
- Data folder inside contains data that can be used to rest scripts
- No data is modified; assignments reference this folder to ensure reproducible execution.

Skills Developed Across All Projects
------------------------------------
1. Bash + Unix-Based Data Processing
- Navigating and manipulating large biological files
- Using grep, sed, awk, pipes, and regex
- Text processing at scale (FASTA, BLAST, multi-line formats)
- Automating workflows with Bash scripts

3. Python Programming for Bioinformatics
- Parsing FASTA, FASTQ, BLAST, and SAM formats
- Building reusable packages
- Handling biological strings, indices, quality scores, and alignments
- Designing algorithms from scratch (e.g., Needleman–Wunsch)
- Creating command-line interfaces with argparse
- Writing maintainable, well-structured code
  
4. Sequence Analysis & Algorithms
- Core bioinformatics algorithms implemented include:
  - isPCR — identify primer matches and simulate amplification
  - Needleman–Wunsch global alignment — dynamic programming, scoring, traceback
  - SAM parsing — interpreting CIGAR strings, flags, mapping direction, and base calls
  - Consensus sequence reconstruction — integrating multiple reads
 
5. Data Analytics & Visualization
- Assignments incorporate:
  - Data manipulation
  - Plotting with Python (Matplotlib)
  - Quantitative reasoning
  - Interpretation of biological datasets
