import numpy as np
from Bio.Seq import Seq


def find_largest_polypeptide_in_DNA(seq, translationTable=1):
    longest_DNA = ''
    longest_amino_acid_sequence = 0
    Data = []
    for strand, nuc in [(1, Seq(seq)), (-1, Seq(seq).reverse_complement())]:
        # Check all three reading frames in this direction.
        for frame in range(3):
            trans = str(nuc[frame:].translate(translationTable))
            cut_codons = 0
            while 'M' in trans:
                codons_before_Met = trans.find('M')
                cut_codons += codons_before_Met
                trans = trans[codons_before_Met:]
                if '*' in trans:
                    length = trans.find('*') + 1
                    if length > longest_amino_acid_sequence:
                        longest_amino_acid_sequence = length
                        first_bp = frame + 3 * cut_codons
                        last_bp = frame + 3 * cut_codons + 3 * (length)
                        longest_DNA = str(nuc[first_bp:last_bp + 1])
                        Data = [trans, longest_DNA, length, first_bp, last_bp, strand, frame,
                                longest_amino_acid_sequence]
                    trans = trans[length:]
                else:
                    # Ignore sequence M... if ORF extends beyond FASTA
                    trans = ''
    return Data


def complement(s):
    """relationship between two structures"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s)
    letters = [complement[base] for base in letters]
    return ''.join(letters)


def revcomplement(s):
    return complement(s[::-1])


def random_dna_sequence(length):
    """Function generate random DNA with a certain probability"""
    return ''.join(np.random.choice(BASES, p=P) for _ in range(int(length)))


# constants
BASES = ('A', 'C', 'T', 'G')
TABLE = 1
MIN_PRO_LEN = 0
RESULTS = []
while True:
    length, frequency = input('Input length of DNA (numeric between 100 and 1000)'
                              '\nand percentage of GC composition (numeric between 20 and 80)'
                              '\nExample: 300 80').split(" ")
    if not length.isnumeric() or not frequency.isnumeric:
        print("You didn't enter a number. Try again: ")
    elif not 100 <= int(length) <= 1000 or not 20 <= int(frequency) <= 80:
        print("You entered the wrong range"
              "\nFor length range (100-1000), for frequency (20-80),try again: ")
    else:
        break

length = int(length)
frequency = int(frequency)

CG = frequency / 100 / 2
AT = (1 - frequency / 100) / 2
"""Probability. Rule of molecular biology"""
P = [AT, CG, AT, CG]

print(find_largest_polypeptide_in_DNA(random_dna_sequence(length)))
