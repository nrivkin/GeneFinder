# -*- coding: utf-8 -*-
"""
metagene_finder.py
includes functions for finding the longest protien sequence shared by the
metagenome and the nitrogenase sequence
needs amino_acids.py, load.py, and gene_finder

All code is either from SoftDes class or personally written

some of the code is based on information found on the wikipedia page on the
longest substring problem,
https://en.wikipedia.org/wiki/Longest_common_substring_problem

@nrivkin: Noah Rivkin

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
import load
# from gene_finder import get_complement, get_reverse_complement, find_all_ORFs
# from gene_finder import find_all_ORFs_oneframe, find_all_ORFs_both_strands, coding_strand_to_AA # these are useful, though possibly not overly efficient


# Use find_all_ORFs and coding_strand_to_AA to find possible, then compare?
def gen_array(n,m):
    """
    creates a list of n lists, each of which has m zeros
    this can be used for the dynamic method of finding substrings
    """
    array = []
    for i in range(n):
        array.append([])
        for j in range(m):
            array[i].append(0)
    return array


def find_common_sub(str1,str2):
    """
    finds the longest common substring from str1 and str2
    with dynamic method
    """
    l1 = len(str1)
    l2 = len(str2)
    longest = 0
    longest_loc = []
    matrix = gen_array(l1,l2)
    # marks location in matrix where there is a match with length of common suffix
    for i in range(l1):
        for j in range(l2):
            if str1[i] == str2[j]:
                if i == 0 or j == 0:
                    matrix[i][j] = 1
                    if matrix[i][j] > longest:
                        longest = matrix[i][j]
                else:
                    matrix[i][j] = matrix[i - 1][j -1] + 1
                    if matrix[i][j] > longest:
                        longest = matrix[i][j]
                        longest_loc = i
    common = ''
    for k in range(longest):
        pass
        common = common + str1[longest_loc - k]
    return common

nitrogenase = load_nitrogenase_seq()
metagenome = load_metagenome()
find_common_sub(nitrogenase, metagenome[90][1])
