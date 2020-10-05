import numpy as np


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
    return ''.join(np.random.choice(BASES, p=P) for _ in range(length))


BASES = ('A', 'C', 'T', 'G')  # constants

while True:
    length, frequency = input('Input length of DNA (numeric between 100 and 1000)'
                              '\nand percentage of GC composition (numeric between 20 and 80) '
                              '\nExample: 300 80').split(" ")
    if not length.isnumeric() or not frequency.isnumeric:
        print("You didn't enter a number. Try again: ")
    elif not 100 <= int(length) <= 1000 or not 20 <= int(frequency) <= 80:
        print("You entered the wrong range"
              "\nFor length range (100-1000), for frequency (20-80),try again: ")
    else:
        break

CG = int(frequency) / 100 / 2
AT = (1 - int(frequency) / 100) / 2
"""Probability. Rule of molecular biology"""
P = [AT, CG, AT, CG]

print(random_dna_sequence(int(length)))


