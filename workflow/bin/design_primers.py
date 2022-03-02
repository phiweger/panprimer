#!/usr/bin/env python3
'''
This script ...
'''


import argparse
# https://docs.python.org/3/library/argparse.html
import itertools
import os

import screed
import networkx as nx
from tqdm import tqdm

from panprimer.utils import \
    load_graph, \
    mask_nodes, \
    extract_color, \
    neighbors, \
    orient_nodes_on_genome
from panprimer.wrappers import design_primers


parser = argparse.ArgumentParser(
    description='Find primers in genome graphs')

parser.add_argument('--maxlen', default=4000, type=int, 
    help='Maximum len of PCR product')

parser.add_argument('--minlen', default=400, type=int, 
    help='Minimum len of PCR product')

parser.add_argument('-n', default=10, type=int, 
    help='Number of candidate primer pairs in result (best first)')

parser.add_argument('--outfile', default='primers.csv',
    help='Where to store primers')

parser.add_argument('--graph', required=True,
    help='Genome graph with color annotation [.gfa]')

parser.add_argument('--index', required=True,
    help='Index of the color annotation')

parser.add_argument('--pattern', required=True,
    help='Which genomes to include and which to exclude (0, 1) [.csv]')

parser.add_argument('--debug', action='store_true',
    help='Store some debugging logs')

parser.add_argument('--genome', required=True,
    help='Genome to validate primer positions')

parser.add_argument('--maxpairs', required=True,
    help='Maximum primer site pairs before we stop the recursion madness')

args = parser.parse_args()


maxlen = args.maxlen
minlen = args.minlen
limits = [minlen, maxlen]
npairs = args.n
outfile = args.outfile
fp_pattern = args.pattern
fp_gfa = args.graph
fp_ix = args.index
genome = args.genome
max_pairs = int(args.maxpairs)


'''
Each node in the genome graph has a set of colors, which means that the
represented sequence of DNA is present in this subset of organsims. We want
to be able to assing genomes into an include and exclude set (0 and 1, 
respectively). We do this through a pattern (mask) on the color index of the
genome graph. E.g.

[1, 0, 1, 0] in a genome graph w/ 4 genomes means that we only want to find
primers for the 1st and 3rd genome, that are NOT present on the 2nd and 4th.
'''
is_included_genome = {}
with open(fp_pattern, 'r') as patternfile:
    for line in patternfile:
        binary, path = line.strip().split(',')
        is_included_genome[os.path.basename(path)] = int(binary)

with open(fp_ix, 'r') as ixfile:
    files = [os.path.basename(i) for i in next(ixfile).split(' ')][:-1]
    # -1 .. last entry is a space so last item is '' (empty)
    pattern = [is_included_genome[i] for i in files]

if args.debug:
    with open('log', 'w+') as log:
        log.write(' '.join([str(i) for i in pattern]))


# Parse .gfa
print('Loading graph ...')
G, strands = load_graph(fp_gfa)
# Return candidate nodes, defined as nodes through which all genomes in the
# "include" set pass, but nodes in the "exclude" set don't.
candidates = mask_nodes(G, pattern)
# print(f'{len(candidates)} candidate nodes match pattern')


d = {}
print('Checking candidates against colors ...')
'''
This step looks through the candidates and checks whether they are suitable for
a PCR. For example: Are two nodes close enough so that a polymerase can span
the distance between them for all genomes in the include set?

We also make sure that we can traverse the graph from primer 1 to primer 2
using each color, ie that this path exists for all genomes in the include set.
'''
for color in tqdm(range(len(pattern))):
    if pattern[color] != 1:
        # Discard color bc/ it is not in the "include" set
        # TODO: We could generalize this to n primer sets, e.g. for 3 species
        continue

    # Extract all nodes that are traversed by a single genome
    # G_ = extract_color(G, i)
    pairs, singletons = set(), set()

    # TODO: Stop if a sufficient number of candidates have been found
    for node in candidates:
        # If there are few genomes to exclude, most nodes are potential
        # primer sites. So we cap them here to not search forever.
        if len(pairs) >= max_pairs:
            continue
        seq = G.nodes[node]['sequence']
        # If node is of sufficient length save. We can probably find 2 
        # primers on this segment. If not ...
        if len(seq) > minlen:
            # TODO: only include the singleton if found in all colors
            singletons.add(node)
            continue
        # ... traverse its neighborhood to see if any other candidate is
        # reachible. The polymerase can only span so many nucleotides.
        for n in neighbors(G, node, color, limits):
            if (n != node) and (n in candidates):
                pairs.add(tuple(sorted([node, n])))
                # sorted .. (1, 3) is the same as (3, 1) to us
                # tuple .. hashable, so we can call set() on it
    d[color] = pairs

# Make sure that 0'ed genomes do not appear in results
d = {k: v for k, v in d.items() if pattern[k] == 1}

# https://stackoverflow.com/questions/30773911/union-of-multiple-sets-in-python
valid_pairs = list(set.intersection(*map(set, d.values())))
print(f'Found {len(valid_pairs) + len(singletons)} potential primer region(s)')


sequences = [(G.nodes[n]['sequence'], 'singleton', [n]) for n in singletons]
# Now we make sure the sequence between valid pairs is not too short and on 
# the same strand (we'll glue them together later in a synthetic contig to
# feed into primer3, see below).
for i, j in valid_pairs:
    '''
    See how the path leaving from one primer unitig can be - and +,
    depending on the path:

    grep "2428690" pangenome.gfa
    S   2428690 GGATGTTAA...   DA:Z:2
    L   1686469 +   2428690 +   30M
    L   1686484 -   2428690 -   30M
    L   2428690 -   1686469 -   30M
    L   2428690 +   1686484 +   30M
    '''
    si = G.nodes[i]['sequence']
    sj = G.nodes[j]['sequence']


    # For primer3 design we need sequence fragments on the same strand,
    # see below.
    with screed.open(genome) as file:
        for read in file:
            try:
                seq = orient_nodes_on_genome(read.sequence, si, sj, minlen)
            except KeyError:
                continue
    
    # if not oss:
    #     sj = screed.rc(sj)

    # if invert:
    #     si, sj = sj, si

    # Create an artificial contig by connecting the two fragments w/ some Ns
    # TODO: Potential error here! si and sj need to be in the same 5'-3' order
    # as in the original sequence otherwise we create an artificial inversion
    # (and the PCR won't work bc/ the primers' 3' ends "point away from one
    # another").
    # seq = si + 'N' * minlen + sj
    sequences.append((seq, 'pair', [i, j]))  # i, j are nodes in the graph



# TODO: Check that the primers are unique in the genome. Or make this an option.


# TODO: separate script that shows insert size distribution and annotation and
# coordinates spanned by the primers.
print('Designing primers ...')
designs = []
for seq, type_, nodes in tqdm(sequences):
    singleton = True if type_ == 'singleton' else False
    try:
        u = design_primers(seq, singleton, minlen=minlen, maxlen=maxlen)
        v = '::'.join([str(n) for n in nodes])
        _ = [i.append(v) for i in u]
        designs.append(*u)
    except (KeyError, OSError) as e:  # no valid primers found
        continue
# [['ATCACTGATGGATTTGACGT', 'TACCCCAAAATGGCTAGAAC', 54.76, 55.01, 'NA', 0.251, '4324814::4324817']]
# [['AGGTTGTGTGGTTCGAATC', 'AAGCGGAGATCATACCCTTA', 55.45, 55.41, 'NA', 1.8569, '4324789::4324793']]


# Sort inplace by penalty value, smallest to largest
designs.sort(key=lambda x: x[5])  
# print(designs[:10])


'''
The case occurs where a singleton can also be reached from a smaller fragment, ie it ends up in a "valid pair". Then primer3 finds the best primer in the singleton, but bc/ it is part of a pair "NA" is written as the PCR product length.

['GCCTGTT...', 'CCCGAGC...', 54.99, 54.99, 465, 0.0146, '152624'],
['GCCTGTT...', 'CCCGAGC...', 54.99, 54.99, 'NA', 0.0146, '152624,152646'],

This can happen multiple times, if a singleton is part of more than one valid pair:
['GCCTGTT...', 'CCCGAGC...', 54.99, 54.99, 'NA', 0.0146, '152624,152646'],
['GCCTGTT...', 'CCCGAGC...', 54.99, 54.99, 'NA', 0.0146, '152624,154425'],
'''
# Deduplicate
result = {}
# Note that the results are sorted by penalty bc/ adding to dict since Python
# v3.6 preserves order.
for i, j, *rest, nodes in designs:
    p = tuple(sorted([i, j]))  # p .. pair
    if len(nodes.split('::')) == 1:  # singleton
        result[p] = [i, j, *rest, nodes]
    else:
        # If there is no entry for the primer pair from a singleton, then add
        try:
            _ = result[p]
        except KeyError:
            result[p] = [i, j, *rest, nodes]

print(f'Found {len(result)} primer pair(s), will save the best {npairs if npairs < len(result) else len(result)}.')
with open(outfile, 'w+') as out:
    out.write('fwd,rev,Tm fwd,Tm rev,product,penalty,nodes\n')  # header
    for v in itertools.islice(result.values(), npairs):
        out.write(','.join(str(i) for i in v) + '\n')




