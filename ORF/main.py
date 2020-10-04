import numpy as np

BASES = ('A', 'C', 'T', 'G')
P = (0.2, 0.2, 0.3, 0.3)

def random_dna_sequence(length):
    return ''.join(np.random.choice(BASES, p=P) for _ in range(length))
print(random_dna_sequence(10))
