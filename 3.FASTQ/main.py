from Bio import SeqIO

Q_SCORE, NUMBER_OF_NUC, REST, READ = int(0), int(0), int(0), int(0)

for record in SeqIO.parse(input('Enter the file name = '), "fastq"):
    READ += 1
    for x in record.letter_annotations['phred_quality']:
        if x >= 30:
            Q_SCORE += 1
        else:
            REST += 1

NUMBER_OF_NUC = Q_SCORE + REST

print("Number of reads = ", READ, "\nTotal number NUC = ", NUMBER_OF_NUC, "\nThe average length of the reads = ", NUMBER_OF_NUC // READ, "\nQSCORE = ",
      "{:%}".format(Q_SCORE / NUMBER_OF_NUC))

