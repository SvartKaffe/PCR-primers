from PCR_primers import Sequence
from primer_algorithms import trie_primers


class Primers:
    def __init__(self, sequence: Sequence, length):
        self.frw_primers = trie_primers(sequence.get_frw_sequence(), length)
        self.rvs_primers = trie_primers(sequence.get_rvs_sequence(), length)
        self.trimmed_frw = {}
        self.trimmed_rvs = {}

    def temp_selection(self, temp):
        for primer, values in self.frw_primers:
            if values["length"] != temp:



