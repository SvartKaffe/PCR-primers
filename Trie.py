from primer_algorithms import complement
from Primers import Primers
# Trie taken from https://albertauyeung.github.io/2020/06/15/python-trie.html/
# The recursive function is heavily inspired by https://github.com/volpato30/hamming-d-search/blob/master/trie.py
# and https://albertauyeung.github.io/2020/06/15/python-trie.html/



class TrieNode:
    """A node in the trie structure"""

    def __init__(self, char):
        # the character stored in this node
        self.char = char

        # whether this can be the end of a word
        self.is_end = False

        # a counter indicating how many times a word is inserted
        # (if this node's is_end is True)
        self.counter = 0

        # a dictionary of child nodes
        # keys are characters, values are nodes
        self.children = {}


class Trie:
    """The trie object"""

    def __init__(self):
        """
        The trie has at least the root node.
        The root node does not store any character
        """
        self.root = TrieNode("")

    def insert(self, word):
        """Insert a word into the trie"""
        node = self.root

        # Loop through each character in the word
        # Check if there is no child containing the character, create a new child for the current node
        for char in word:
            if char in node.children:
                node = node.children[char]
            else:
                # If a character is not found,
                # create a new node in the trie
                new_node = TrieNode(char)
                node.children[char] = new_node
                node = new_node

        # Mark the end of a word
        node.is_end = True

        # Increment the counter to indicate that we see this word once more
        node.counter += 1

    def hamming_distance(self, primer: str, delta_t: int):

        result = []
        current_index = 0
        current_cost = 0

        root_node = self.root

        for nucleotide in root_node.children:
            self.recursive_search(root_node.children[nucleotide], result, current_index, current_cost, primer,
                                  nucleotide, delta_t, nucleotide)

        return result

    def recursive_search(self, node: TrieNode, result: list, current_index: int,
                         current_cost: int, primer: str, nucleotide: str, delta_t: int, trie_sequence: str):

        if len(result) > 0:
            return

        if nucleotide == primer[current_index]:
            if primer[current_index] in ("G", "C"):
                current_cost += 4
            if primer[current_index] in ("A", "T"):
                current_cost += 2
        if not (nucleotide == primer[current_index]):
            base, sum = complement(primer[current_index])
            if primer[current_index] == "A" and not nucleotide == base:
                current_cost += 2
            if primer[current_index] == "T" and not nucleotide == base:
                current_cost += 2
            if primer[current_index] == "G" and not nucleotide == base:
                current_cost += 4
            if primer[current_index] == "C" and not nucleotide == base:
                current_cost += 4

        if current_cost > delta_t:
            return

        current_index += 1

        # stopping criteria
        if node.is_end:
            if trie_sequence != primer:
                result.append(trie_sequence)
            else:
                return

        for nucleotide in node.children:
            self.recursive_search(node.children[nucleotide], result, current_index, current_cost, primer, nucleotide,
                                  delta_t, trie_sequence + nucleotide)


if __name__ == "__main__":
    t = Trie()
    t.insert("GATCA")
    t.insert("GATCG")
    t.insert("ATGCC")
    t.insert("AGCGT")

    debug = t.hamming_distance("GGGGG", 10)
    print(debug)