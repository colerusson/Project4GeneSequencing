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


def unrestricted_alignment(s1, s2):
    # Initialize the matrix with zeros
    matrix = [[0 for x in range(len(s2) + 1)] for y in range(len(s1) + 1)]

    # Fill in the first row and column with gap penalties
    for i in range(1, len(s1) + 1):
        matrix[i][0] = i * -5
    for j in range(1, len(s2) + 1):
        matrix[0][j] = j * -5

    # Fill in the rest of the matrix with edit distances
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            substitution_cost = 1 if s1[i - 1] != s2[j - 1] else -3
            matrix[i][j] = max(matrix[i - 1][j - 1] + substitution_cost,
                               matrix[i - 1][j] - 5,
                               matrix[i][j - 1] - 5)

    # Traceback to get the aligned sequences
    aligned_s1 = ""
    aligned_s2 = ""
    i = len(s1)
    j = len(s2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (1 if s1[i - 1] != s2[j - 1] else -3):
            aligned_s1 = s1[i - 1] + aligned_s1
            aligned_s2 = s2[j - 1] + aligned_s2
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] - 5:
            aligned_s1 = s1[i - 1] + aligned_s1
            aligned_s2 = "-" + aligned_s2
            i -= 1
        else:
            aligned_s1 = "-" + aligned_s1
            aligned_s2 = s2[j - 1] + aligned_s2
            j -= 1

    return matrix[-1][-1], aligned_s1, aligned_s2


def banded_alignment(seq1, seq2, MaxCharactersToAlign):
    m, n = len(seq1), len(seq2)
    k = 7  # bandwidth is 2k+1
    F = [[0] * n for _ in range(m)]
    for i in range(m):
        for j in range(max(0, i - k), min(n, i + k + 1)):
            if i == 0 and j == 0:
                F[i][j] = 0
            elif i == 0:
                F[i][j] = F[i][j - 1] - 5
            elif j == 0:
                F[i][j] = F[i - 1][j] - 5
            else:
                match = -3 if seq1[i] == seq2[j] else 1
                F[i][j] = max(F[i - 1][j - 1] + match, F[i][j - 1] - 5, F[i - 1][j] - 5)
    score = F[m - 1][n - 1]
    alignment1, alignment2 = '', ''
    i, j = m - 1, n - 1
    while i >= 0 and j >= 0:
        if i > 0 and j > 0 and F[i][j] == F[i - 1][j - 1] - 3 and seq1[i] == seq2[j]:
            alignment1 = seq1[i] + alignment1
            alignment2 = seq2[j] + alignment2
            i -= 1
            j -= 1
        elif j > 0 and F[i][j] == F[i][j - 1] - 5:
            alignment1 = '-' + alignment1
            alignment2 = seq2[j] + alignment2
            j -= 1
        elif i > 0 and F[i][j] == F[i - 1][j] - 5:
            alignment1 = seq1[i] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        else:
            break
    alignment1 = alignment1[MaxCharactersToAlign]
    alignment2 = alignment2[MaxCharactersToAlign]
    return score, alignment1, alignment2


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, seq1, seq2, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length

        ###################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
        # score = random.random()*100
        # alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
        #	len(seq1), align_length, ',BANDED' if banded else '')
        # alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
        #	len(seq2), align_length, ',BANDED' if banded else '')

        if not self.banded:
            score, alignment1, alignment2 = unrestricted_alignment(seq1[:align_length], seq2[:align_length])
        # Perform banded alignment
        else:
            score, alignment1, alignment2 = banded_alignment(seq1[:align_length], seq2[:align_length], self.MaxCharactersToAlign)

        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
