from Bio import SeqIO
from Bio.Restriction import *
import re
import time


"""
This file contains functions used in different files, some are not used in the main
file since they are part of the brute-force solution.
"""


def find_primers(sequence: str, temp: int, reverse=False) -> dict:
    """
    Initially used for the brute-force solution, but can be used in main program. Genereates primers based on temperature,
    not length.
    :param sequence: DNA string
    :param temp: Tm for the primers
    :param reverse: used to get the correct start/stop positions when called for the reverse strand.
    :return: a dictionary containing potential primers
    """
    i = 0
    primers = {}
    duplicate_primers = []
    pattern = re.compile("(C|G)(C|G)(C|G)")
    while i <= len(sequence):
        sum = 0
        primer = ""
        for base in sequence[i:]:
            if base == "A":
                if sum + 2 <= temp:
                    sum += 2
                    primer += base
                else:
                    break
            elif base == "T":
                if sum + 2 <= temp:
                    sum += 2
                    primer += base
                else:
                    break
            elif base == "G":
                if sum + 4 <= temp:
                    sum += 4
                    primer += base
                else:
                    break
            elif base == "C":
                if sum + 4 <= temp:
                    sum += 4
                    primer += base
                else:
                    break

        # Primer computations and conditions
        primer_length = len(primer)
        if reverse:
            start = len(sequence) - i
            stop = start + 1 - primer_length
        else:
            start = i+1
            stop = i+primer_length

        g_count = primer.count("G")
        c_count = primer.count("C")
        primer_end = primer[-3:]
        m = pattern.match(primer_end)

        if g_count > 0 and c_count > 0:
            gc_ratio = (g_count + c_count) / primer_length
        else:
            gc_ratio = 0

        primer_conditions = (
                primer_length == 20
                and (primer[-1] == "G" or primer[-1] == "C")
                and sum == temp
                and (0.4 <= gc_ratio <= 0.6)
                and not m
        )

        if primer_conditions:
            if (primer in primers.keys()) and (primer not in duplicate_primers):
                duplicate_primers.append(primer)
            else:
                primers[primer] = {"length": primer_length,
                                   "start": start,
                                   "stop": stop,
                                   "temp": sum
                                   }

        i += 1
    # Removes any primers found multiple times in the sequence
    for i in duplicate_primers:
        if i in primers.keys():
            primers.pop(i)

    return primers


def unique_primers(primers: dict, delta_t: int) -> dict:
    """
    aligns all primers to all other primers in the dictionary and only keeps the primers outside the specified
    delta_t value. This functions was really never used, was intended for the brute force solution.
    :param primers: dictionary containing primers
    :param delta_t: Ta
    :return: dictionary with primers
    """
    new_primers = {}
    for primer, values in primers.items():
        good_primer = True
        for primer2 in primers:
            if not primer == primer2:
                i = 0
                temp_difference = 0
                min_value = min(len(primer), len(primer2)) - 1
                while i <= min_value:
                    if primer[i] != primer2[i]:
                        if primer2[i] == ("A" or "T"):
                            temp_difference += 2
                        else:
                            temp_difference += 4

                    if temp_difference > delta_t:
                        break

                    i += 1
                if temp_difference < delta_t:
                    good_primer = False
                    break

        if good_primer:
            new_primers[primer] = values

    return new_primers


def primer_alignment(frw_sequence: str, rvs_sequence: str, primers: dict, delta_t: int, temp: int) -> dict:
    """
    This function is used for the brute force solution. It aligns all primers to the genome and discards the primer
     if it binds within the specified Ta (delta_t) value.
    :param frw_sequence: string of the forward sequence
    :param rvs_sequence: string of the reverse sequence
    :param primers: dictionary of primers
    :param delta_t: user specified Ta value
    :param temp: Tm for the primers in the dictionary
    :return: dictionary of possible primer pairs.
    """
    good_primers = {}
    melting_temp = temp - delta_t
    for primer, values in primers.items():
        appearances = 0
        i = 0
        window = len(primer)
        while i <= len(frw_sequence) and appearances <= 1:
            frw_binding_temp = 0
            rvs_binding_temp = 0
            frw_sequence_slice = frw_sequence[i:i + window]
            rvs_sequence_slice = rvs_sequence[i:i + window]
            if len(frw_sequence_slice) != window:
                break
            for j in range(len(primer)):
                nucleotide, score = complement(primer[j])
                if frw_sequence_slice[j] == nucleotide:
                    frw_binding_temp += score
                if rvs_sequence_slice[j] == nucleotide:
                    rvs_binding_temp += score
            if frw_binding_temp > melting_temp:
                appearances += 1
            if rvs_binding_temp > melting_temp:
                appearances += 1
            i += 1
        if appearances <= 1:
            good_primers[primer] = values

    return good_primers


def trie_primers(sequence, length: int, reverse=False) -> dict:
    """
    This function divides the genome into primers of size length.
    :param sequence: Sequence object
    :param length: length of the primers
    :param reverse: used to generate the correct start/stop positions if the reverse complement is used.
    :return: dictionary containing primers of size length
    """
    i = 0
    primers = {}
    duplicate_primers = []
    while i <= len(sequence):
        primer = sequence[i:i+length]
        primer_conditions = (
                len(primer) == length
        )
        primer_length = len(primer)
        if reverse:
            start = len(sequence) - i
            stop = start + 1 - primer_length
        else:
            start = i+1
            stop = i+primer_length

        if primer_conditions:
            if (primer in primers.keys()) and (primer not in duplicate_primers):
                duplicate_primers.append(primer)
            if (primer not in primers.keys()) and (primer not in duplicate_primers):
                primers[primer] = {"length": primer_length,
                                   "start": start,
                                   "stop": stop,
                                   }
        i += 1
    # Removes any primers found multiple times in the sequence
    for i in duplicate_primers:
        if i in primers.keys():
            primers.pop(i)

    return primers


def complement(nucleotide: str) -> str and int:
    """
    Simple function used in various different parts of the program. Returns the complement and the annealing temp
    cost.
    :param nucleotide: A DNA base (G, C, A, T)
    :return: The complement to the nucleotide as well as the Marmor Doty annealing cost.
    """
    if nucleotide == "A":
        return "T", 2
    if nucleotide == "T":
        return "A", 2
    if nucleotide == "G":
        return "C", 4
    if nucleotide == "C":
        return "G", 4


def temp_calc(string: str) -> int:
    """
    Simple function used to calculate the Tm value for a DNA sequence.
    :param string: DNA string
    :return: the calculated Tm value.
    """
    temp = 0
    for i in string:
        if i in ("A", "T"):
            temp += 2
        if i in ("C", "G"):
            temp += 4
    return temp


def search(trie, primers: dict, delta_t: int) -> dict:
    """
    This function calls the recursive search function hamming_distance.
    :param trie: Trie object
    :param primers: dictionary of primers
    :param delta_t: user specified Ta value
    :return: a dictionary containing primers outside the specified Ta (delta_t) value
    """
    good_primers = {}
    for primer, values in primers.items():
        result = trie.hamming_distance(primer, delta_t)
        if len(result) == 0:
            good_primers[primer] = values
    return good_primers


def sort_primers(frw_primers: dict, rvs_primers: dict, sequence) -> "two lists of tuples":
    """
    This function pairs forward with reverse primers, it really slows down the program if it is given many
    primers due to its O(n^2) complexity and other computations. This Function is not used anymore.
    :param frw_primers: dictionary of forward primers
    :param rvs_primers: dictionary of reverse primers
    :param sequence: Sequence object
    :return: Two lists of tuples, one containing "normal" primer pairs as well as other useful information. The other
    contains primer pairs that go over the beginning/end aka "circular primers" along with other useful information.
    """
    primer_pairs = []
    circular_pairs = []
    DNA = sequence.get_frw_sequence()
    DNA_length = len(DNA)

    for frw_primer, frw_value in frw_primers.items():
        for rvs_primer, rvs_value in rvs_primers.items():
            start = frw_value["start"]
            stop = rvs_value["start"]

            fragment_length = (stop - start + 1)
            circular_length = (DNA_length - start + stop + 1)
            fragments = EcoRI(DNA[start-1:stop])
            circular_DNA = DNA[start-1:]+DNA[:stop]
            circular_fragments = EcoRI(circular_DNA)

            conditions = (len(fragments) >= 2
                          and fragments[0] >= 300
                          and fragments[-1] >= 300
                          and fragment_length <= 3000
                          )

            circular_conditions = (len(circular_fragments) >= 2
                                   and circular_length <= 3000
                                   and circular_fragments[0] >= 300
                                   and circular_fragments[-1] >= 300
                                   )

            if conditions:
                pair_tuple = (frw_primer, rvs_primer, start, stop, fragment_length, fragments)
                primer_pairs.append(pair_tuple)

            if circular_conditions:
                circular_tuple = (frw_primer, rvs_primer, start, stop, circular_length, circular_fragments)
                circular_pairs.append(circular_tuple)

    return primer_pairs, circular_pairs


def forward(frw_primers: dict, rvs_primers: dict, max_size: int, min_size: int):
    """
    Pairs forward primers with reverse primers that form a fragment of suitable length.
    :param frw_primers: dictionary of primers
    :param rvs_primers: dictionary of primers
    :param max_size: max size of fragment
    :param min_size: min size of fragment
    :return: list of primer pairs along with other information
    """
    primer_pairs = []

    for frw_primer, frw_value in frw_primers.items():
        for rvs_primer, rvs_value in rvs_primers.items():
            start = frw_value["start"]
            stop = rvs_value["start"]
            fragment_length = (stop - start + 1)

            if min_size <= fragment_length <= max_size:
                primer_pairs.append([frw_primer, rvs_primer, start, stop, fragment_length])

    return primer_pairs


def circular(frw_primers: dict, rvs_primers: dict, sequence, max_size: int, min_size: int):
    """
    Pairs forward with reverse primers, this looks for circular primer pairs.
    :param frw_primers: dictionary of primers
    :param rvs_primers: dictionary of primers
    :param sequence: Sequence object used to get the sequence length
    :param max_size: max fragment size
    :param min_size: min fragment size
    :return: list of circular primer pairs along with other information
    """
    primer_pairs = []
    dna_length = sequence.sequence_length

    for frw_primer, frw_value in frw_primers.items():
        for rvs_primer, rvs_value in rvs_primers.items():
            start = frw_value["start"]
            stop = rvs_value["start"]
            circular_length = (dna_length - start + stop + 1)

            if min_size <= circular_length <= max_size:
                primer_pairs.append([frw_primer, rvs_primer, start, stop, circular_length])

    return primer_pairs


def EcoRI_digest(primer_pairs: list, sequence, circular=False):
    """
    This function digest a fragment, generated from a primer pairs, using EcoRI.
    :param primer_pairs: list of primer pairs
    :param sequence: sequence object
    :param circular: is set to True to get right fragment if the pair is circular
    :return: list of primer pairs along with the digestion map generated.
    """
    dna = sequence.get_frw_sequence()
    length = len(primer_pairs)

    if circular:
        for i in range(length):
            start = primer_pairs[i][2]
            stop = primer_pairs[i][3]
            dna_fragment = dna[start-1:] + dna[:stop]
            fragments = EcoRI(dna_fragment)
            primer_pairs[i].append(fragments)
    else:
        for i in range(length):
            start = primer_pairs[i][2]
            stop = primer_pairs[i][3]
            fragments = EcoRI(dna[start-1:stop])
            primer_pairs[i].append(fragments)

    return primer_pairs


def EcoRI(sequence) -> list[int]:
    """
    This function generates the digestion map for a DNA sequence, used in the main program.
    :param sequence: Sequence object
    :return: a list with the fragments from the EcoRI digestion.
    """
    sekvens = sequence
    enzymes = RestrictionBatch(["EcoRI"])
    result = enzymes.search(sekvens)
    locations = list(result.values())
    sekvens = str(sequence)
    location_single = [i-1 for sublist in locations for i in sublist]
    fragments = [len(sekvens[v1:v2]) for v1, v2 in zip([0]+location_single, location_single+[None])]

    return fragments


if __name__ == "__main__":
    for sequence1 in SeqIO.parse("SARS-CoV-2-isolate-Wuhan-Hu-1-complete-genome.fasta", "fasta"):
        sequence2 = sequence1.seq
        sequence3 = sequence1.reverse_complement().seq

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

    test2 = "AAGCCTCGGCAATGTACTACATTCGGTAC"

    start = time.time()

    forward = find_primers(sequence2, 60)
    reverse = find_primers(sequence3, 60, reverse=True)

    print(len(forward), "number of forward primers")
    print(len(reverse), "number of reverse primers")

    forward_alignment = primer_alignment(sequence2, sequence3, forward, 24, 60)
    reverse_alignment = primer_alignment(sequence2, sequence3, reverse, 24, 60)

    print(len(forward_alignment), "number of forward primers after alignment")
    print(len(reverse_alignment), "number of reverse primers after alignment")

    end = time.time()
    print(end - start)
    print("\n")

    print(forward_alignment)
    print(reverse_alignment)


