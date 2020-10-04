import numpy as np

BASES = ('A', 'C', 'T', 'G')  # constants
print(20 / 100 / 2, 20 / 100 / 2, (1 - 20 / 100) / 2, (1 - 20 / 100) / 2)

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

CG = frequency / 100 / 2
AT = (1 - frequency / 100) / 2
P = [AT, CG, AT, CG]  # probability

def random_dna_sequence(length):
    """Function generate random DNA with a certain probability"""
    return ''.join(np.random.choice(BASES, p=P) for _ in range(length))


print(random_dna_sequence(10))
