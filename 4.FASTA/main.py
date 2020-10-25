from Bio import SeqIO

LIST = [len(rec) for rec in SeqIO.parse(input('Enter the file name = '), "fasta")]
LIST.sort()

N50 = int(0)
FULL_LENGTH = sum(LIST)
COUNT = len(LIST)
AVE_LENGTH = FULL_LENGTH / COUNT
temp = int(0)

for i in LIST:
    temp += i
    if temp > FULL_LENGTH / 2:
        N50 = i
        break

print("COUNT = ", COUNT, "\nMIN L = ", min(LIST), "\nMAX L =", max(LIST), "\nAVERAGE L = ", AVE_LENGTH, "\nN50 = ", N50)
