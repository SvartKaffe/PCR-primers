from Bio.Seq import Seq
from Bio import SeqIO


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

    def primers(self):
        print()


if __name__ == "__main__":
    object1 = Sequence()
    object1.identify_primers()
    object1.read_sequence("test.fasta")



