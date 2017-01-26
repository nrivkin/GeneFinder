# -*- coding: utf-8 -*-
"""
gene_finder.py
includes functions for finding ORFs in genes
needs amino_acids.py and load.py

All code is either from SoftDes class or personally written

@nrivkin: Noah Rivkin

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == 'A':
        return('T')
    elif nucleotide == 'T':
        return('A')
    elif nucleotide == 'C':
        return('G')
    elif nucleotide == 'G':
        return('C')
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    length = len(dna)
    rev_comp = ''
    for i in range(length):
        rev_comp = get_complement(dna[i]) + rev_comp
    return rev_comp
    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    sequence = ''
    for index in range(len(dna)):
        if index % 3 == 0 and index < len(dna) - 2:
            codon = dna[index] + dna[index + 1] + dna[index + 2]
            if codon in codons[10]:
                return sequence
            else:
                sequence += codon
    return dna
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    ORFs = []
    index = 0
    while index < len(dna) - 2:
        codon = dna[index] + dna[index + 1] + dna[index + 2]
        if codon in codons[3]:
            ORF = rest_of_ORF(dna[index:len(dna)])
            ORFs.append(ORF)
            index += len(ORF)
        else:
            index += 3
    return(ORFs)
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for index in range(3):
        test = dna[index:len(dna)]
        all_ORFs = all_ORFs + find_all_ORFs_oneframe(test)
    return all_ORFs
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    all_ORFs_both_strands = []
    all_ORFs_both_strands = all_ORFs_both_strands + find_all_ORFs(dna)
    second_strand = get_reverse_complement(dna)
    all_ORFs_both_strands = all_ORFs_both_strands+find_all_ORFs(second_strand)
    return all_ORFs_both_strands
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    longest = ''
    ORFs = find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF) > len(longest):
            longest = ORF
    return longest
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = 0
    for i in range(num_trials):
        shuffled = shuffle_string(dna)
        length = len(longest_ORF(shuffled))
        if length > longest:
            longest = length
    return longest
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    index = 0
    a_strand = ''
    while index < len(dna) - 2:
        codon = dna[index] + dna[index + 1] + dna[index + 2]
        a_strand = a_strand + aa_table[codon]
        index += 3
    return a_strand
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    sequences = []
    for ORF in ORFs:
        if len(ORF) > threshold:
            sequences.append(coding_strand_to_AA(ORF))
    return sequences
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()


dna = load_seq("./data/X73525.fa")
gene_finder(dna)
