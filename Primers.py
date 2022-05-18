from Sequence import Sequence
from primer_algorithms import trie_primers, temp_calc
import re
from Bio import SeqIO
from Bio.SeqUtils import GC

class Primers:
    def __init__(self, sequence: Sequence, length):
        self.frw_primers = trie_primers(sequence.get_frw_sequence(), length)
        self.rvs_primers = trie_primers(sequence.get_rvs_sequence(), length)

    def temp_selection(self, temp):
        for primer, values in list(self.frw_primers.items()):
            if temp_calc(primer) != temp:
                self.frw_primers.pop(primer)

        for primer, values in list(self.rvs_primers.items()):
            if temp_calc(primer) != temp:
                self.rvs_primers.pop(primer)

    def GC_clamp(self):
        pattern = re.compile("(C|G)(C|G)(C|G)")
        for primer, values in list(self.frw_primers.items()):
            three_end = primer[-3:]
            one_end = primer[-1]
            m = pattern.match(str(three_end))
            conditions = (
                one_end in ("G", "C")
                and not m
            )
            if not conditions:
                self.frw_primers.pop(primer)

        for primer, values in list(self.rvs_primers.items()):
            three_end = primer[-3:]
            one_end = primer[-1]
            m = pattern.match(str(three_end))
            conditions = (
                    one_end in ("G", "C")
                    and not m
            )
            if not conditions:
                self.rvs_primers.pop(primer)

    def GC_content(self, upper_limit: int, lower_limit: int):
        for primer in list(self.frw_primers.keys()):
            if not (lower_limit <= GC(primer) <= upper_limit):
                self.frw_primers.pop(primer)

        for primer in list(self.rvs_primers.keys()):
            if not (lower_limit <= GC(primer) <= upper_limit):
                self.rvs_primers.pop(primer)

    def get_frw_primers(self):
        return self.frw_primers

    def get_rvs_primers(self):
        return self.rvs_primers


if __name__ == "__main__":
    object1 = Sequence("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta")
    primers = Primers(object1, 20)
    print(len(primers.get_frw_primers()))
    print(len(primers.get_rvs_primers()))
    primers.GC_clamp()
    print(len(primers.get_frw_primers()))
    print(len(primers.get_rvs_primers()))
    primers.GC_content(60, 40)
    print(len(primers.get_frw_primers()))
    print(len(primers.get_rvs_primers()))
    primers.temp_selection(60)
    print(len(primers.get_frw_primers()))
    print(len(primers.get_rvs_primers()))

