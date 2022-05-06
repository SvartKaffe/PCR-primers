from Bio.Seq import Seq
from Bio import SeqIO
import re
from primer_algorithm import find_primers


class Sequence:
    def __init__(self):
        self.frw_sequence = None
        self.frw_primers = []
        self.rvs_sequence = None
        self.rvs_primers = []
        self.sequence_length = 0

    def read_sequence(self, file):
        for sequence in SeqIO.parse(file, "fasta"):
            self.frw_sequence = sequence.seq
            self.rvs_sequence = sequence.reverse_complement()
            self.sequence_length = len(sequence)

    def primers(self, temp):
        self.frw_primers = find_primers(self.frw_sequence, temp)
        self.rvs_primers = find_primers(self.rvs_sequence, temp, reverse=True)


if __name__ == "__main__":
    object1 = Sequence()
    object1.read_sequence("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta")
    object1.primers(60)
    print(object1.frw_primers, "\n")
    print(object1.rvs_primers)
