def unrestricted_alignment(seq1, seq2):
    # Initialize the matrix
    matrix = [[0 for i in range(len(seq2) + 1)] for j in range(len(seq1) + 1)]

    # Initialize the first row and column
    for i in range(len(seq1) + 1):
        matrix[i][0] = i * INDEL
    for j in range(len(seq2) + 1):
        matrix[0][j] = j * INDEL

    # Fill in the rest of the matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            # Compute the score for a match, mismatch, and indel
            match = matrix[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB)
            indel1 = matrix[i - 1][j] + INDEL
            indel2 = matrix[i][j - 1] + INDEL

            # Take the minimum of the three scores
            matrix[i][j] = min(match, indel1, indel2)

    # Compute the alignment
    i = len(seq1)
    j = len(seq2)
    alignment1 = ''
    alignment2 = ''
    while i > 0 and j > 0:
        # If the score is a match, add the characters to the alignment
        if matrix[i][j] == matrix[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else SUB):
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = seq2[j - 1] + alignment2
            i -= 1
            j -= 1
        # If the score is an indel, add a gap to the alignment
        elif matrix[i][j] == matrix[i - 1][j] + INDEL:
            alignment1 = seq1[i - 1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j - 1] + alignment2
            j -= 1

    # Add the remaining characters to the alignment
    while i > 0:
        alignment1 = seq1[i - 1] + alignment1
        alignment2 = '-' + alignment2
        i -= 1
    while j > 0:
        alignment1 = '-' + alignment1
        alignment2 = seq2[j - 1] + alignment2
        j -= 1

    return matrix[len(seq1)][len(seq2)], alignment1, alignment2




def banded_alignment(seq1, seq2):
    # Initialize the matrix
    k = 2 * MAXINDELS + 1
    n = len(seq1) + 1
    matrix = [[0 for i in range(n)] for j in range(k)]

    # Initialize the first row
    for j in range(1, n):
        matrix[0][j] = j * INDEL

    # Fill in the rest of the matrix
    for i in range(1, k):
        for j in range(i, n):
            # Compute the score for a match, mismatch, and indel
            if seq1[j - 1] == seq2[j - i]:
                score = MATCH
            elif seq1[j - 1].lower() == seq2[j - i].lower():
                score = SUB / 2
            else:
                score = SUB
            match = matrix[i - 1][j - 1] + score
            indel1 = matrix[i - 1][j] + INDEL
            indel2 = matrix[i][j - 1] + INDEL

            # Take the minimum of the three scores
            matrix[i][j] = min(match, indel1, indel2)

    # Compute the alignment
    i = k - 1
    j = n - 1
    alignment1 = ''
    alignment2 = ''
    while i > 0 and j > 0:
        # If the score is a match, add the characters to the alignment
        if matrix[i][j] == matrix[i - 1][j - 1] + (MATCH if seq1[j - 1] == seq2[j - i] else SUB):
            alignment1 = seq1[j - 1] + alignment1
            alignment2 = seq2[j - i] + alignment2
            i -= 1
            j -= 1
        # If the score is an indel, add a gap to the alignment
        elif matrix[i][j] == matrix[i - 1][j] + INDEL:
            alignment1 = seq1[j - 1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j - i] + alignment2
            j -= 1

    # Add the remaining characters to the alignment
    while i > 0:
        alignment1 = seq1[j - 1] + alignment1
        alignment2 = '-' + alignment2
        i -= 1
    while j > 0:
        alignment1 = '-' + alignment1
        alignment2 = seq2[j - i] + alignment2
        j -= 1

    return matrix[k - 1][n - 1], alignment1, alignment2

