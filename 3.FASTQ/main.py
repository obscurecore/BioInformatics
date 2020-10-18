from Bio import SeqIO

Q_SCORE, NUMBER_OF_NUC, REST, READ = int(0), int(0), int(0), int(0)

for record in SeqIO.parse("test_example3.fastq", "fastq"):
    READ += 1
    for x in record.letter_annotations['phred_quality']:
        if x >= 30:
            Q_SCORE += 1
        else:
            REST += 1

NUMBER_OF_NUC = Q_SCORE + REST

print("1 = ", READ, "\n2 = ", NUMBER_OF_NUC, "\n3 = ", NUMBER_OF_NUC // READ, "\n4 = ",
      "{:%}".format(Q_SCORE / NUMBER_OF_NUC))

