from Bio import SeqIO
from Bio.Restriction import *
from Bio.Seq import Seq


for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
    sequence2 = sequence1.seq
    sequence3 = sequence1.reverse_complement().seq

sekvens = Seq("ATTCGGAATTCACGA")
enzymes = RestrictionBatch(["EcoRI"])
result = enzymes.search(sekvens)
print(result)
string = str(sekvens)
print(string[0:6], string[6:])
print(type(result))

locations = list(result.values())
print(locations)
location_single = [i-1 for sublist in locations for i in sublist]
print(location_single)

fragments = [sekvens[v1:v2] for v1, v2 in zip([0]+location_single, location_single+[None])]




def EcoRI(sequence: str):
    sekvens = Seq(sequence)
    enzymes = RestrictionBatch(["EcoRI"])
    result = enzymes.search(sekvens)
    locations = list(result.values())
    location_single = [i-1 for sublist in locations for i in sublist]
    fragments = [sequence[v1:v2] for v1, v2 in zip([0]+location_single, location_single+[None])]

    return fragments

print("\n")

abc = EcoRI("ATTCGGAATTCACGA")
print(abc[1])
print(abc[0])
print(abc[-1])