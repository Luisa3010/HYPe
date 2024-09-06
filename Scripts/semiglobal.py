#### create semiglobal alignment of two sequences - made by Vincent Hanlon



import numpy as np

def semiglobal_matrix(query, ref, match=1, mismatch=-1, gap=-2):
# Fill in  a matrix for semi-global alignment of two sequences
# (A modification of needleman-wunsch)

    nrows=len(query)+1
    ncols=len(ref)+1
    matrix=np.zeros(shape=[nrows, ncols])

    # Whereas the first row needs to be zeros (to allow gaps in front of the query,
    # ...the first column needs successive gap penalties, representing successive
    # deletions of the beginning of the query
    for row in range(1, nrows):
        matrix[row, 0]=matrix[row-1, 0] + gap

    # Filling in the rest of the matrix
    for row in range(1, nrows):
        for col in range(1, ncols):

            # Standard NW procedure: moving from top or left incurs a gap
            topscore=matrix[row-1,col] + gap
            leftscore=matrix[row,col-1] + gap

            # A diagonal move incurs a match or mismatch, depending on the sequence
            if query[row-1] == ref[col-1]:
                diagscore=matrix[row-1,col-1] + match
            else:
                diagscore=matrix[row-1,col-1] + mismatch

            # Each cell is filled with the best-scoring route
            matrix[row,col]=max(topscore, leftscore, diagscore)

    # In the resulting matrix, the last row (minus the first column) gives the best
    # alignment ending at each position in the ref sequence
    return(matrix)
