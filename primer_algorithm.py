from Bio.Seq import Seq
from Bio import SeqIO

# TODO: save primers in dictionary together with their start index
# TODO: should the index for the reverse strand be reversed?
# random comment


def find_primers(sequence, delta_t):
    i = 0
    primers = []
    while i <= len(sequence):
        sum = 0
        primer = ""
        for base in sequence[i:]:
            if base == "A":
                if sum + 2 <= delta_t:
                    sum += 2
                    primer = primer + base
                else:
                    break
            elif base == "T":
                if sum + 2 <= delta_t:
                    sum += 2
                    primer = primer + base
                else:
                    break
            elif base == "G":
                if sum + 4 <= delta_t:
                    sum += 4
                    primer = primer + base
                else:
                    break
            elif base == "C":
                if sum + 4 <= delta_t:
                    sum += 4
                    primer = primer + base
                else:
                    break
        if primer not in primers:
            primers.append(primer)
        i += 1
    return primers


if __name__ == "__main__":
    for sequence1 in SeqIO.parse("Enterobacteria-phage-P2-NC_001895-complete-genome.fasta", "fasta"):
        sequence2 = sequence1.seq
    test = "ATGTCAATGGCTATCGTCACTATT"
    a = find_primers(sequence2, 10)
    print(a)
