import numpy as np
from Bio.Seq import Seq
import regex as re


def find_largest_polypeptide_in_DNA(seq, translationTable=1):
    allPossibilities = []
    for frame in range(3):
        print("CHECK OUT")
        print(str(Seq(seq[frame:]).translate(translationTable)))
        print("STOP")
    """
      framePossibilitiesF = [i[i.find("M"):] for i in trans.split("*") if "M" in i]
      allPossibilities += framePossibilitiesF
  allPossibilitiesLengths = [len(i) for i in allPossibilities]

  if len(allPossibilitiesLengths) == 0:
      raise Exception("no candidate ORFs")

  proteinAsString = allPossibilities[allPossibilitiesLengths.index(max(allPossibilitiesLengths))]

  return Seq(proteinAsString, alphabet=ProteinAlphabet)
"""


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



find_largest_polypeptide_in_DNA(random_dna_sequence(length))
seq = Seq(random_dna_sequence(length))
print("ASDADSADSADSADSADASDSADSADASDSAD")
startP = re.compile('ATG')
max_len = int(0);
for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
    for frame in range(3):
        for pro in nuc[frame:].translate(TABLE).split("*"):
            if len(pro) >= max_len:
                max_len = len(pro)
                RESULTS = [pro[:30], pro[-3:], len(pro), strand, frame]
                # print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))
print(seq)
print(*RESULTS)
