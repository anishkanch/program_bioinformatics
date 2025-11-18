

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    """
    Needleman-Wunsch algorithm - maximizing the number of identical bases that line up in the sequences while minimizing gaps and mismatches
    Will return aligned sequence 1, 2, and a score
    """

#Obtain the lengths of the two sequences
    a = len(seq_a)
    b = len(seq_b)

#Initialization: construct a matrix
    rows = b+1
    cols = a+1
    matrix = []
    for r in range(rows):
        row = []
        for c in range(cols):
            row.append(0)

        matrix.append(row)
    
#Initialization: fill top row and left column as gaps
    gap = -1
    for i in range(1, cols):
        matrix[0][i] = i * gap
    for i in range(1, rows):
        matrix[i][0] = i * gap
    
#Matrix-filling: fill the matrix
    match = 1
    mismatch = -1
    gap = -1

    for r in range(1, rows):
        for c in range(1, cols):
            #diagonal, match bonus or mismatch penalty
            if seq_a[c-1] == seq_b[r-1]:
                diagonal = matrix[r-1][c-1] + match
            else:
                diagonal = matrix[r-1][c-1] + mismatch
            #up with gap penalty
            up = matrix[r-1][c] + gap
            #left with gap penalty 
            left = matrix[r][c-1] + gap
            
            #greatest our of diagnonal, up, left will fill
            matrix[r][c] = max(diagonal, up, left)
    
#Back-tracing: Follow the arrows back
    align_a = ""
    align_b = ""

    #start at max value (bottom right) and work backwards
    r = rows -1
    c = cols -1

    while r>0 or c>0: #top left corner 
        if r>0 and c>0: #both seq should have letters to compare
            score_current = matrix[r][c]
            score_diag = matrix[r-1][c-1]

            if seq_a[c-1] == seq_b[r-1]: #match bonus added
                score_diag += match
            else:
                score_diag += mismatch #mismatch penalty
            
            #going diagonal as tie breaker
            if score_current == score_diag:
                #reverses, adding as we go backwards
                align_a = seq_a[c-1] + align_a
                align_b = seq_b[r-1] + align_b
                #move up and do the left
                r = r-1
                c= c-1
                continue
            
            #if no diagonal, can go in two directions

            #left means seq_b has a gap, add the dash as a gap
            if c>0 and (r==0 or matrix[r][c] == matrix[r][c-1] + gap):
                align_a = seq_a[c-1] + align_a
                align_b = '-' + align_b
                c= c-1
            #up means seq_a has a gap, add the dash as a gap
            else:
                align_a = '-' + align_a
                align_b = seq_b[r-1] + align_b
                r = r-1
    
#Finishing: alignment score
    score = matrix[rows - 1][cols - 1] #final alignment score is in the bottom right 
    return (align_a, align_b), score    
    
