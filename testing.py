from Bio import SeqIO
from primer_algorithms import temp_calc
import re

for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
    sequence2 = sequence1.seq
    sequence3 = sequence1.reverse_complement().seq

primer = "ATGGGC"  # temp should be 20
a = temp_calc(primer)

pattern = re.compile("(C|G)(C|G)(C|G)")

s = sequence2[:1]
print(s)

print(s == "G")

