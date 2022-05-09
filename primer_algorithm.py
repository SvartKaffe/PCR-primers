from Bio.Seq import Seq
from Bio import SeqIO
import re


def find_primers(sequence, temp, reverse=False):
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
    for item in duplicate_primers:
        if item in primers.keys():
            primers.pop(item)

    return primers, duplicate_primers


def unique_primers(sequence, primers, delta_t):
    print()


if __name__ == "__main__":
    for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
        sequence2 = sequence1.seq
    test = "ATGTCAATGGCTATCGTCACTATTGTACTTAGCCAATTAG"
    a,b = find_primers(sequence2, 60)
    print(len(a))
    print('ACCAACAGCACGTCTGAAAC' in a.keys())


