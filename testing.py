from Bio import SeqIO
from Bio.Restriction import *
from Bio.Seq import Seq
from Sequence import Sequence

for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
    sequence2 = sequence1.seq
    sequence3 = sequence1.reverse_complement().seq


def EcoRI(sequence):
    sekvens = sequence
    enzymes = RestrictionBatch(["EcoRI"])
    result = enzymes.search(sekvens)
    locations = list(result.values())
    sekvens = str(sequence)
    location_single = [i-1 for sublist in locations for i in sublist]
    fragments = [sekvens[v1:v2] for v1, v2 in zip([0]+location_single, location_single+[None])]

    return fragments

print("\n")

ecori_test = Sequence("test.fasta")
print(ecori_test.get_frw_sequence()[1])


