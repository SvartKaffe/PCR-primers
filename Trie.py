from primer_algorithms import complement, temp_calc
# Trie taken from https://albertauyeung.github.io/2020/06/15/python-trie.html/
# The recursive function is heavily inspired by https://github.com/volpato30/hamming-d-search/blob/master/trie.py
# and https://albertauyeung.github.io/2020/06/15/python-trie.html/


class TrieNode:
    """
    This class is used for the nodes in the Trie class.
    """

    def __init__(self, char: str):
        """
        Constructor for the TrieNode class, stores information about the node, such as which base it is, if it is an end
        node and its child nodes.
        :param char: A DNA base
        """
        # the character stored in this node
        self.char = char

        # whether this can be the end of a word
        self.is_end = False

        # remnant from source, not used
        self.counter = 0

        # a dictionary of child nodes
        # keys are characters, values are nodes
        self.children = {}


class Trie:
    """
    The trie datatype.
    """

    def __init__(self):
        """
        The trie class constructor, builds the first root node, which is an empty string, the number
        of nodes are also stored here.
        """
        self.root = TrieNode("")
        self.num_nodes = 0

    def insert(self, word: str):
        """
        Insert method, inserts a string into the trie. It first looks if the letter exists, if not,
        it creates a new node and assigns the letter to it.
        :param word: string
        :return: nothing
        """
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
                self.num_nodes += 1

        # Mark the end of a word
        node.is_end = True

        # Increment the counter to indicate that we see this word once more
        node.counter += 1

    def hamming_distance(self, primer: str, delta_t: int) -> list:
        """
        Initial function that calls the recursive search function. Still called hamming_distance despite not calculating
        the hamming distance, the name from the source repo stuck.
        :param primer: DNA primer that will be used to traverse the trie
        :param delta_t: Ta (delta_t) value.
        :return: a list, either empty or with one item in it.
        """

        result = []
        current_index = 0
        current_cost = 0

        root_node = self.root

        for nucleotide in root_node.children:
            self.recursive_search(root_node.children[nucleotide], result, current_index, current_cost, primer,
                                  nucleotide, delta_t)

        return result

    def recursive_search(self, node: TrieNode, result: list, current_index: int,
                         current_cost: int, primer: str, nucleotide: str, delta_t: int):
        """
        The recursive part of the search function. Iterates through all the nodes in the trie unless aborted early.

        :param node: current node in the trie
        :param result: a list, used to stop the search
        :param current_index: current index in the primer
        :param current_cost: the current delta_t (Ta) value for the primer
        :param primer: the primer which is being aligned to the genome
        :param nucleotide: Nucleotide of the current/next node
        :param delta_t: the Ta value (delta_t)
        :return: nothing
        """

        # needed to successfully terminate search
        if len(result) > 0:
            return

        # base at current_index in primer is being compared to the base in the node.
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

        # checks if the current primer has a bigger delta_t value than user specified one,
        # if true, it ends the search of the current branch
        if current_cost > delta_t:
            return

        current_index += 1

        # if at the end of a branch, current_cost = 0, a per
        # it means that the current primer has a location in the genome which is within the specified
        # delta_t value (Ta), it is added to the result list and search is aborted for the primer
        if node.is_end:
            if current_cost == 0:
                return
            else:
                result.append(1)

        # continuation of the recursive search
        for nucleotide in node.children:
            self.recursive_search(node.children[nucleotide], result, current_index, current_cost, primer, nucleotide,
                                  delta_t)


if __name__ == "__main__":
    t = Trie()
    t.insert("GATCA")
    t.insert("GATCG")
    t.insert("ATGCC")
    t.insert("AGCGT")
    t.insert("TCTGA")

    debug = t.hamming_distance("GTACT", 10)
    print(debug)
"""
        if node.is_end:
            if current_cost == 0:
                return
            else:
                result.append(1)
"""