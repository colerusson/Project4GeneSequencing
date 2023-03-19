#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


# implement the unrestricted alignment using the Needleman-Wunsch algorithm
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

    # Add the remaining characters to the   alignment
    while i > 0:
        alignment1 = seq1[i - 1] + alignment1
        alignment2 = '-' + alignment2
        i -= 1
    while j > 0:
        alignment1 = '-' + alignment1
        alignment2 = seq2[j - 1] + alignment2
        j -= 1

    return matrix[len(seq1)][len(seq2)], alignment1, alignment2


# implement the banded alignment using the Needleman-Wunsch algorithm
def banded_alignment(seq1, seq2):
    # TODO: Implement the banded alignment
    # Compute the bandwidth
    bandwidth = 2 * MAXINDELS + 1

    # Initialize the matrix
    matrix = [[0 for i in range(bandwidth)] for j in range(len(seq1) + 1)]

    # Initialize the first row
    for i in range(len(seq1) + 1):
        matrix[i][0] = i * INDEL

    # Fill in the rest of the matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, bandwidth):
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


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean
    # that tells you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        ###################################################################################################
        if not self.banded:
            score, alignment1, alignment2 = unrestricted_alignment(seq1[:align_length], seq2[:align_length])
        # Perform banded alignment
        else:
            score, alignment1, alignment2 = banded_alignment(seq1[:align_length], seq2[:align_length])
        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
