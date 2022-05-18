from Bio.Seq import Seq
from Bio import SeqIO
from primer_algorithms import trie_primers, find_primers
from Trie import Trie, TrieNode
import time


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
    start = time.time()
    object1 = Sequence()
    object1.read_sequence("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta")
    trie = object1.build_trie(20)
    forward = find_primers(object1.frw_sequence, 60)
    reverse = find_primers(object1.rvs_sequence, 60, reverse=True)
    print(len(forward))
    print(len(reverse))

    good_frw_primers = {}
    good_rvs_primers = {}
    for primer, values in forward.items():
        result = trie.hamming_distance(primer, 18)
        if len(result) == 0:
            good_frw_primers[primer] = values
    for primer, values in reverse.items():
        result = trie.hamming_distance(primer, 18)
        if len(result) == 0:
            good_rvs_primers[primer] = values
    end = time.time()
    print(end-start)
    print(f"good forward primers:{len(good_frw_primers)}, good reverse primers: {len(good_rvs_primers)}")



