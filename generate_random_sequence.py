from random import choice


def random_dna(length):
    dna = ""
    for count in range(length):
        dna += choice("CGTA")
    return dna


if __name__ == "__main__":
    print(random_dna(1000))
