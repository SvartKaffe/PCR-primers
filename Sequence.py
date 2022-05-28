from Bio import SeqIO
from primer_algorithms import trie_primers, find_primers
from Trie import Trie, TrieNode
import time


class Sequence:
    """
    The Sequence class is used to read the given genome and to build the trie. It can also return the forward and
    reverse complement of the given genome.
    """
    def __init__(self, file):
        """
        When called it reads the sequence from the given fasta file.
        :param file: file or the file location of the fasta file
        """
        for sequence in SeqIO.parse(file, "fasta"):
            self.frw_sequence = sequence.seq
            self.rvs_sequence = sequence.reverse_complement().seq

    def get_frw_sequence(self) -> "biopython Seq object":
        """
        Returns the forward sequence.
        :return: the forward sequence
        """
        return self.frw_sequence

    def get_rvs_sequence(self) -> "biopython Seq object":
        """
        Returns the reverse complement of the sequence
        :return: the reverse complement
        """
        return self.rvs_sequence

    def build_trie(self, length: int) -> "Trie object":
        """
        Calls the trie_primers() function to generate primers of size length and usees them to build the trie.
        :param length: length of the primers
        :return: A Trie object
        """
        forward_primers = trie_primers(self.frw_sequence, length)
        reverse_primers = trie_primers(self.rvs_sequence, length)
        trie = Trie()
        # add primers to trie
        for i in forward_primers.keys():
            trie.insert(i)
        for j in reverse_primers.keys():
            trie.insert(j)

        return trie


if __name__ == "__main__":
    test = "TGACTGACTCACGGTCGTTTGTGCACGGCTTATCGCTAACCGGTGTCTGCGCACCCGGTCAATCTTTAGCGACAATACACAACCTGGTTGACAATCGCTATGCT" \
           "GGCGTTTTCACCCTATTTGTACCGACCATAGAGGGCGCCTGCGTACTCGAGGAAAAAGACCTCCTACCCCTCATGATTCGAGTCGCCCGACCTCAACAATTCCTA" \
           "TTATAGGTTGATCTACTCGAGTCACACTAACCCGTTACTCGAACACGGATTGGTCGGATCATCACGGTGCAAGTAGGACTAAAAAATAACAGTAGTCAATAATAT" \
           "CAACTTAGGGTTGATACGTTACCCATGCATGTTAGACCTACTAAGCTTTCCTGGCAGCCCGTTCTGTACGAAAGTGACTGAGCTTGCCTGCCCGAACATGTCCAC" \
           "GTTGAGCGAGTTTTTAACTCCGTCAGGTTGCCCGTATAGGTGGTATTGCTATCCGCAATAATGGGCGCCAGGTAATTAACGCATATAACCGACCCCCCGCGCTCCG"\
           "GTTGTTATAAGCCTCGGCAATGTACTACATTCGGTATCGCTGGAACCTTTATGCCACCAAGAGTCCGACTTCCCCGTGTTGCGAGACGCACGCATGACACCGGGT" \
           "ACTCGTTGGCAGGAATTTGCTGAGTTGTACTGTTTGTACGTTGTCCTCTTATAACTTTCAAGTTACAACTTTTCAAATGCTACGAGGGTTAAACTTATGGTGTAT" \
           "CAGAATAAACAATCAGGCGCCCTCGAAACAAGTGAACTTTTATAATACGATAAAACGATTCATGAGCCATTTACTGGTTCGCCCCTCAATACCAGCGCCGCTAAA" \
           "ATGTCTGATGACCCACGTGCTATTTTTCTAGCGCAGTCATAGTCACGGGACAAAGCTCCAACGTAGAGGACGGGCTAAGTCGCGTAAAACGCCAACCGAAGTCAA" \
           "GGACCCCGTCATGTTAGCATCCGTATACGCGCATGCGCGAAGGGCTGTCATTTCC"

    start = time.time()
    object1 = Sequence("SARS-CoV-2-isolate-Wuhan-Hu-1-complete-genome.fasta")

    trie = object1.build_trie(20)
    forward = find_primers(object1.frw_sequence, 60,)
    reverse = find_primers(object1.rvs_sequence, 60,  reverse=True)
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

    primer_pairs = []
    for frw_primer, frw_value in good_frw_primers.items():
        for rvs_primer, rvs_value in good_rvs_primers.items():
            fragment = (rvs_value["start"] - frw_value["start"] + 1)
            if (300 <= fragment <= 2000):
                pair_tuple = (frw_primer, rvs_primer, fragment)
                primer_pairs.append(pair_tuple)




