aa = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
      '|', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R',
      'G']

codons = [['TTT', 'TTC'],                              # 0
          ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],  # 1
          ['ATT', 'ATC', 'ATA'],                       # 2
          ['ATG'],                                     # 3     start codon
          ['GTT', 'GTC', 'GTA', 'GTG'],                # 4
          ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],  # 5
          ['CCT', 'CCC', 'CCA', 'CCG'],                # 6
          ['ACT', 'ACC', 'ACA', 'ACG'],                # 7
          ['GCT', 'GCC', 'GCA', 'GCG'],                # 8
          ['TAT', 'TAC'],                              # 9
          ['TAA', 'TAG', 'TGA'],                       # 10    stop codon
          ['CAT', 'CAC'],                              # 11
          ['CAA', 'CAG'],                              # 12
          ['AAT', 'AAC'],                              # 13
          ['AAA', 'AAG'],                              # 14
          ['GAT', 'GAC'],                              # 15
          ['GAA', 'GAG'],                              # 16
          ['TGT', 'TGC'],                              # 17
          ['TGG'],                                     # 18
          ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  # 19
          ['GGT', 'GGC', 'GGA', 'GGG']]                # 20

# create a dictionary lookup table for mapping codons into amino acids
aa_table = {}
for i in range(len(aa)):
    for codon in codons[i]:
        aa_table[codon] = aa[i]
