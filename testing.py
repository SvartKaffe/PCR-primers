from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers
import time

# test run to see if shit works as intended
# timer
start = time.time()

# read in genome
object1 = Sequence("10k.fasta")

# build trie
trie = object1.build_trie(20)

# generate primers
primers = Primers(object1, 20)
primers.GC_clamp()
primers.GC_content(60, 40)
primers.temp_selection(60)

# go through the trie
forward_primers = search(trie, primers.get_frw_primers(), 20)
reverse_primers = search(trie, primers.get_rvs_primers(), 20)

print(len(forward_primers))
print(len(reverse_primers))

primer_pairs, circular_pairs = sort_primers(forward_primers, reverse_primers, object1)

print(len(primer_pairs))
print(len(circular_pairs))
print(primer_pairs)
print(circular_pairs)
print("\n")
end = time.time()
print(f"It took {end-start} seconds to finish the program, from start to finish")

""" 
a.sort(key=lambda y: y[-1], reverse=True)
use this to sort the list in the end
"""