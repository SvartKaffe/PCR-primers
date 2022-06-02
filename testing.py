from Sequence import Sequence
from Primers import Primers
from primer_algorithms import search, sort_primers, forward, circular, EcoRI_digest
import time

"""
This file is used for testing purposes.
"""


# test run to see if shit works as intended
# timer
start_read_build = time.time()
# read in genome
object1 = Sequence("fasta_files\\40k.fasta")
# build trie
trie = object1.build_trie(20)
end_read_build = time.time()

start_filter = time.time()
# generate primers
primers = Primers(object1, 20)
primers.GC_clamp()
primers.GC_content(60, 40)
primers.temp_selection(60)
end_filter = time.time()

start_search = time.time()
# go through the trie
forward_primers = search(trie, primers.get_frw_primers(), 18)
reverse_primers = search(trie, primers.get_rvs_primers(), 18)
end_search = time.time()

print("\n")
print(len(forward_primers), len(reverse_primers))
print("\n")

start_pair = time.time()

primer_pairs = forward(forward_primers, reverse_primers)
circular_pairs = circular(forward_primers, reverse_primers, object1)
primer_pairs = EcoRI_digest(primer_pairs, object1)
circular_pairs = EcoRI_digest(circular_pairs, object1)

end_pair = time.time()

print(len(primer_pairs), len(circular_pairs))
print("\n")

end = time.time()

read_build = end_read_build - start_read_build
filters = end_filter - start_filter
searches = end_search - start_search
pairing = end_pair - start_pair
total_primers = int(len(primers.get_frw_primers()) + len(primers.get_rvs_primers()))
trie_len = trie.num_nodes
tot_time = end - start_read_build

print(f"It took {end-start_read_build} seconds to finish the program, from start to finish")
print([read_build, filters, searches, pairing, total_primers, trie_len, tot_time])
""" 
a.sort(key=lambda y: y[-1], reverse=True)
use this to sort the list in the end
"""