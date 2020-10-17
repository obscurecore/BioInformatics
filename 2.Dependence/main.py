import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
import re

Dict, BASES, LENGHT, COUNT, StartGC, EndGC, TABLE, MIN_PRO_LEN = {}, ('A', 'C', 'T', 'G'), int(100), int(10), int(
    20), int(80), int(1), int(10)

PATTERN = re.compile("ATG(?:...).*?(?:TAG|TGA|TAA)")


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
            ORF = PATTERN.findall(str(nuc[frame:]))
            if ORF and len(max(ORF)) > MIN_PRO_LEN:
                return True
            return False


# Main
for freq in range(StartGC, EndGC):
    x = 0
    for i in range(COUNT):
        if find_orf(Seq(random_dna_sequence(LENGHT, freq))):
            x += 1
    Dict[freq] = x / COUNT * 100

# Histogram
plt.plot(list(Dict.keys()), Dict.values(), color='g')
plt.show()
