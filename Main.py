from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers
import time

# test run to see if shit works as intended
# timer
start = time.time()

# read in genome
object1 = Sequence("SARS-CoV-2-isolate-Wuhan-Hu-1-complete-genome.fasta")

# build trie
trie = object1.build_trie(20)

# generate primers
primers = Primers(object1, 20)
primers.GC_clamp()
primers.GC_content(60, 40)
primers.temp_selection(60)

# go through the trie
forward_primers = search(trie, primers.get_frw_primers(), 18)
reverse_primers = search(trie, primers.get_rvs_primers(), 18)

print(len(forward_primers))
print(len(reverse_primers))

primer_pairs = sort_primers(forward_primers, reverse_primers, object1)

print(primer_pairs)
end = time.time()
print(f"It took {end-start} seconds to finish the program, from start to finish")