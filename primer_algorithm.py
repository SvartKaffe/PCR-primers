from Bio.Seq import Seq
from Bio import SeqIO
import re

# TODO: save primers in dictionary together with their start index, length etc
# TODO: should the index for the reverse strand be reversed?


def find_primers(sequence, temp):
    i = 0
    primers = []
    pattern = re.compile("(C|G)(C|G)(C|G)")
    while i <= len(sequence):
        sum = 0
        primer = ""
        for base in sequence[i:]:
            if base == "A":
                if sum + 2 <= temp:
                    sum += 2
                    primer = primer + base
                else:
                    break
            elif base == "T":
                if sum + 2 <= temp:
                    sum += 2
                    primer = primer + base
                else:
                    break
            elif base == "G":
                if sum + 4 <= temp:
                    sum += 4
                    primer = primer + base
                else:
                    break
            elif base == "C":
                if sum + 4 <= temp:
                    sum += 4
                    primer = primer + base
                else:
                    break

        primer_length = len(primer)
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
                and primer not in primers
                and (primer[-1] == "G"
                or primer[-1] == "C")
                and (0.4 <= gc_ratio <= 0.6)
                and not m
        )

        if primer_conditions:
            primers.append(primer)
        i += 1
    return primers


if __name__ == "__main__":
    for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
        sequence2 = sequence1.seq
    test = "ATGTCAATGGCTATCGTCACTATT"
    a = find_primers(sequence2, 60)
    print(a)
