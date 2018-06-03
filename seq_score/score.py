#!/usr/bin/python3
import sys
import math
from Bio.SubsMat import MatrixInfo

# check if the correct number of arguments have been passed
if len(sys.argv) < 5:
    print("Usage:")
    print("python3 score.py <seq 1 with - for gap> <seq 2 with - for gap> <gap open penalty(positive)> <gap extension penalty(positive)>")
    print("Note: this algorithm uses affine gap scoring")
    print("Note: if one sequence is shorter than the other, only the paired base pairs and gaps will be scored.")
    exit(0)

# we continue!
seq1 = sys.argv[1] # sequence 1
seq2 = sys.argv[2] # sequence 2
d = float(sys.argv[3]) # open gap penalty
e = float(sys.argv[4]) # extension gap penalty
score = 0
curr_gap = False
for i in range(min(len(seq1, seq2))): # loop over every base pair, and add to the accumulative score
    if seq1[i:i+1] == '-' or seq2[i:i+1] == '-': # test for a gap
        if curr_gap: # if we are in the middle of a gap, subtract the extension penalty
            score -= e
            continue
        else: # if we are starting a gap, subtract the open penalty
            score -= d
            curr_gap = True # make sure to remember that we are in the middle of a gap
            continue

    # if there is no gap, look up the base pair in the BLOSUM50 matrix
    score += MatrixInfo.blosum50[(seq1[i:i+1], seq2[i:i+1])]

print("SCORE: ", str(int(score)))