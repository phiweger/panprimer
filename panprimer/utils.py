from collections import defaultdict
from itertools import product
import random
import re
from typing import Dict, List

import networkx as nx
import numpy as np
import screed
from screed import rc as revcomp
from tqdm import tqdm


def load_graph(path):
    '''Load graph from gfa file.

    Sequences are nodes, they are linked by an edge if they overlap:

    > Links are the primary mechanism to connect segments. Links connect oriented segments. A link from A to B means that the end of A overlaps with the start of B. If either is marked with -, we replace the sequence of the segment with its reverse complement, whereas a + indicates the segment sequence is used as-is. -- http://gfa-spec.github.io/GFA-spec/GFA1.html
    '''
    G = nx.Graph()
    strands = {}

    with open(path, 'r') as file:
        for line in tqdm(file):
            # S       10      GATGGAGA...ACG        DA:Z:1  CO:Z:1000000000000
            if line[0] == 'S':
                _, node, seq, _, weight = line.strip().split('\t')
                node = int(node)
                w = np.array([int(i) for i in weight.replace('CO:Z:', '')])
                G.add_node(node, **{'colors': w, 'sequence': seq})
            # L 53508   +   3225    -   30M
            if line[0] == 'L':
                _, node1, strand1, node2, strand2, _ = line.strip().split('\t')
                node1, node2 = [int(i) for i in (node1, node2)]
                G.add_edge(node1, node2)
                strands[(node1, node2)] = (strand1, strand2)
    
    cc = len(list(nx.connected_components(G)))
    print(f'The graph has {cc} connected component(s)')
    return G, strands


def mask_nodes(G: nx.classes.graph.Graph, pattern: List[int]) -> list:
    '''Apply a selection scheme to all node colors.

    Say a node is present in three genomes [1, 1, 1, 0, 0, 0]. The fn will
    return all nodes that present this pattern. It is thus possible to define
    which genomes to include and exclude, e.g. in the case where we only want
    to include a certain strain from a larger species.
    
    Usage:

    colors = np.array([1, 1, 1, 0, 0, 0])
    '''
    assert len(random.choice(G.nodes)['colors']) == len(pattern), 'The length of the pattern must match the number of colors in the graph.'

    candidates = set()
    for i in G:
        if all(G.nodes[i]['colors'] == pattern):
            candidates.add(i)
    print(f'Found {len(candidates)} candidate nodes ({round(len(candidates) / len(G), 4)}) that fit the specified pattern.')
    return candidates


def extract_color(G: nx.classes.graph.Graph, ix: int) -> nx.classes.graph.Graph:
    '''Subset graph by color index'''
    G_ = G.copy()
    nodes = list(G_.nodes)
    for i in nodes:
        if G_.nodes[i]['colors'][ix] == 0:
            G_.remove_node(i)
    return G_


def greedy_set_cover(universe: set, sets: Dict[int, set]) -> List[set]:
    '''
    https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect09-set-cover.pdf

    > Given a set of elements {1,2,...,n} (called the universe) and a collection S of m sets whose union equals the universe, the set cover problem is to identify the smallest sub-collection of S whose union equals the universe. -- https://en.wikipedia.org/wiki/Set_cover_problem

    Primer pairs are subsets, and they contain the colors. Minimum number of primer pairs to contain all colors.

    n .. size of universe .. number of colors

    - http://www.martinbroadhurst.com/greedy-set-cover-in-python.html
    - https://github.com/guangtunbenzhu/SetCoverPy

    > There is a well-known greedy approximation algorithm for set cover that is also easy to implement in whatever language of your choice. The algorithm itself is described here: https://en.wikipedia.org/wiki/Set_cover_problem#Greedy_algorithm -- https://stackoverflow.com/questions/7936037/any-good-implementation-of-greedy-set-cover-for-large-datasets

    > Unfortunately the Wikipedia entry doesn't actually cover weighted set cover, which is the case here. The extension is simple, and is described e.g. here: http://pages.cs.wisc.edu/~shuchi/courses/880-S07/scribe-notes/lecture03.pdf -- https://stackoverflow.com/questions/7936037/any-good-implementation-of-greedy-set-cover-for-large-datasets

    > Some more useful notes: http://www.cs.ucr.edu/~neal/non_arxiv/Young08SetCover.pdf http://www.cs.uiuc.edu/class/sp08/cs473/Lectures/lec20.pdf

    [Implementation](https://stackoverflow.com/questions/7942312/how-do-i-make-my-implementation-of-greedy-set-cover-faster)

    # Case where we don't mask unwanted genomes
    d = {}
    for n in G:
        pattern = G.nodes[n]['colors']
        d[n] = set([ix for ix, i in enumerate(pattern) if i == 1])

    greedy_set_cover({0, 1, 2, 3, 4, 5, 6, 7}, d)
    '''
    # TODO: evaluate here: using DFS, can we find a node x bases away
    # that has the same set members (and is valid as a primer base) IF
    # that segment is not already large enough to place two primers
    cover = []
    u = universe.copy()  # so we don't consume the original set

    while u:
        d = defaultdict(list)
        for k, v in sets.items():
            d[len(u.intersection(v))].append(k)
        
        # Randomly select a set that covers most items in the universe
        # TODO: filter in the next step, e.g. by length or whatnot
        n = random.choice(d[max(d)])
        cover.append(n)
        # Remove the covered items from the universe
        for i in sets[n]:
            try:
                u.remove(i)
            except KeyError:
                continue
    
    return cover


def dfs(G, node, visited, pathlen, limits, color):
    '''
    Do a depth first search but stop when the cumulative length of the node
    sequences is above some limit.
    '''
    minlen = limits[0]
    maxlen = limits[1]
    if node not in visited:
        # Make sure the region spanned between two nodes is not too short
        if minlen < pathlen[node]:
            visited.append(node)
        
        for i in nx.neighbors(G, node):
            pathlen[i] = pathlen[node] + len(G.nodes[i]['sequence'])
            try:
                # Some connected components seem to be colorless:
                # if all(G.nodes[i]['colors'] == pattern):
                # KeyError: 'colors'
                if (pathlen[i] < maxlen) and (G.nodes[i]['colors'][color] == 1):
                    dfs(G, i, visited, pathlen, limits, color)
                else:
                    continue
            except KeyError:
                continue
    return visited


def neighbors(G, node, color, limits):
    pathlen = {}
    pathlen[node] = 0
    visited = []
    return dfs(G, node, visited, pathlen, limits, color)


def orient_nodes_on_genome(seq, p1, p2, minlen):
    pairs = list(product([p1, revcomp(p1)], [p2, revcomp(p2)]))
    
    for pair in pairs:
        r1, r2 = [re.search(subseq, seq) for subseq in pair]
        if all([r1, r2]):
            '''
            This is the orientation we are looking for. Now arrange correct
            5'-3' sequence of the two subsequences (so the designed primers
            don't start transciption in opposing directions due to an
            artificial inversion).
            '''

            # coords = r1.start(), r1.end(), r2.start(), r2.end()
            invert = True if (r1.start() > r2.start()) else False
            # Ns = (r2.start() - r1.end()) * 'N'
            Ns = minlen * 'N'

            if not invert:
                return pair[0] + Ns + pair[1]
            else:
                return pair[1] + Ns + pair[0]
                
            # Only report the first instance (there could be multiple
            # bindings of primers on a query genome).
            # return oss, invert
            # mn, mx = min(coords), max(coords)
            # return o, mn, mx, seq[mn:mx]
        else:
            continue

    raise KeyError('Substrings not found in sequence')




