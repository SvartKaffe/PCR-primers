from Bio.Seq import Seq
from Bio import SeqIO
import re
from primer_algorithms import find_primers, unique_primers


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

    def build_trie(self):



if __name__ == "__main__":
    object1 = Sequence()
    object1.read_sequence("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta")
    object1.primers(60)
    print(len(object1.frw_primers), len(object1.rvs_primers))
    print(list(object1.frw_primers.keys())[1], list(object1.rvs_primers.keys())[1])
    #object1.align_primer(10)
    #print(len(object1.frw_primers), len(object1.rvs_primers))
