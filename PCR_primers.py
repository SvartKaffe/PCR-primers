from Bio.Seq import Seq
from Bio import SeqIO
from primer_algorithms import trie_primers
from Trie import Trie, TrieNode


class Sequence:
    def __init__(self):
        self.frw_sequence = None
        self.frw_primers = []
        self.rvs_sequence = None
        self.rvs_primers = []

    def read_sequence(self, file):
        for sequence in SeqIO.parse(file, "fasta"):
            self.frw_sequence = sequence.seq
            self.rvs_sequence = sequence.reverse_complement()

    def build_trie(self, length):
        forward_primers = trie_primers(self.frw_sequence, length)
        reverse_primers = trie_primers(self.rvs_sequence, length)
        trie = Trie()
        # add sequences to trie
        for i in forward_primers:
            trie.insert(i)
        for j in reverse_primers:
            trie.insert(j)

        return trie


if __name__ == "__main__":
    object1 = Sequence()
    object1.read_sequence("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta")
    object1.build_trie(20)

