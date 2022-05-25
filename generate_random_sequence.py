from random import choice
from Bio import SeqIO


def random_dna(length):
    dna = ""
    for count in range(length):
        dna += choice("CGTA")
    return dna


if __name__ == "__main__":
    sequence = random_dna(30000)
    fasta = ">123123123123123 \n" + sequence

    write_file = open(r"30k.fasta", "w+")
    write_file.write(fasta)
    write_file.close()




