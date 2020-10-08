import numpy as np
from Bio.Seq import Seq


def random_dna_sequence(length, frequency):
    """Function generate random DNA with a certain probability"""
    CG = frequency / 100 / 2
    AT = (1 - frequency / 100) / 2
    """Probability. Rule of molecular biology"""
    P = [AT, CG, AT, CG]
    return ''.join(np.random.choice(BASES, p=P) for _ in range(int(length)))


def find_orf(seq):
    for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            if len(max(nuc[frame:].translate(TABLE).split("*"))) > MIN_PRO_LEN:
                return True


# constants
BASES = ('A', 'C', 'T', 'G')
LENGHT = int(1000)
RESULTS = []
MIN_PRO_LEN = int(10)
TABLE = 1
freq = int(30)


print(random_dna_sequence(LENGHT, freq))
print("SPACE")
find_orf(Seq(random_dna_sequence(LENGHT, freq)))
