from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils.lcc import lcc_mult, lcc_simp
import re
import time

# TODO: look at this https://www.youtube.com/watch?v=QkwPf3fcxBs
# TODO: remove low complexity regions from sequence, lower amount of potential primers

def find_primers(sequence, temp, reverse=False):
    # TODO: make DNA circular
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
                primer_length > temp//4
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


def unique_primers(primers, delta_t):
    # checks all primers against all other primers and sees if they are within delta_t of each other
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


def primer_alignment(sequence, primers, delta_t, temp):
    good_primers = {}
    melting_temp = temp - delta_t
    for primer, values in primers.items():
        appearances = 0
        i = 0
        window = len(primer)
        while i <= len(sequence) and appearances <= 1:
            binding_temp = 0
            sequence_slice = sequence[i:i+window]
            if len(sequence_slice) != window:
                break
            for j in range(len(primer)):
                nucleotide, score = complement(primer[j])
                if sequence_slice[j] == nucleotide:
                    binding_temp += score
            if binding_temp > melting_temp:
                appearances += 1
            i += 1
        if appearances <= 1:
            good_primers[primer] = values

    return good_primers


def complement(nucleotide):
    if nucleotide == "A":
        return "T", 2
    if nucleotide == "T":
        return "A", 2
    if nucleotide == "G":
        return "C", 4
    if nucleotide == "C":
        return "G", 4


if __name__ == "__main__":
    for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
        sequence2 = sequence1.seq

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
    a = find_primers(sequence2, 60)
    print(len(a))
    b = primer_alignment(sequence2, a, 20, 60)
    print(len(b))
    end = time.time()
    print(end - start)



