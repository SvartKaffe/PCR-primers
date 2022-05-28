from primer_algorithms import trie_primers, temp_calc
import re
from Bio.SeqUtils import GC
from Sequence import Sequence


class Primers:
    """
    The Primer class is used to store primers and to apply filters to them.
    """
    def __init__(self, sequence, length: int):
        """
        Generates all the potential primers of size length.
        :param sequence: Sequence object
        :param length: length of the primers to be generated
        """
        self.frw_primers = trie_primers(str(sequence.get_frw_sequence()), length)
        self.rvs_primers = trie_primers(str(sequence.get_rvs_sequence()), length, reverse=True)

    def temp_selection(self, temp: int):
        """
        Discards all primers that do not have Tm = temp.
        :param temp: the Tm temperature of the primers
        :return: nothing
        """
        for primer, values in list(self.frw_primers.items()):
            if temp_calc(primer) != temp:
                self.frw_primers.pop(primer)

        for primer, values in list(self.rvs_primers.items()):
            if temp_calc(primer) != temp:
                self.rvs_primers.pop(primer)

    def GC_clamp(self):
        """
        Discards all primers that do not have G/C clamp and the primers that have a three-peat of G/C,
        or any combination of the two bases and the last three positions of the primer.
        :return: nothing
        """
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
        """
        This method discards all primers that do not satisfy the G/C content criteria.
        :param upper_limit: upper limit of G/C content
        :param lower_limit: lower limit of G/C content
        :return: nothing
        """
        for primer in list(self.frw_primers.keys()):
            if not (lower_limit <= GC(primer) <= upper_limit):
                self.frw_primers.pop(primer)

        for primer in list(self.rvs_primers.keys()):
            if not (lower_limit <= GC(primer) <= upper_limit):
                self.rvs_primers.pop(primer)

    def get_frw_primers(self) -> dict:
        """
        This method is used to return the forward primers.
        :return: dictionary containing the current state of the forward primers
        """
        return self.frw_primers

    def get_rvs_primers(self) -> dict:
        """
        This method returns the reverse primers.
        :return: returns the current state of the reverse primers
        """
        return self.rvs_primers


if __name__ == "__main__":
    object1 = Sequence("test.fasta")
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

    frw_primers = primers.get_frw_primers()
    print(type(frw_primers))
    key_list = list(frw_primers.keys())
    print(key_list[0])
    print(type(key_list[0]))
    print(frw_primers)


