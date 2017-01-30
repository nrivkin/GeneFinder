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
from gene_finder import get_complement, get_reverse_complement, find_all_ORFs
from gene_finder import find_all_ORFs_oneframe, find_all_ORFs_both_strands, coding_strand_to_AA # these are useful, though possibly not overly efficient


# Use find_all_ORFs and coding_strand_to_AA to find possible, then compare?