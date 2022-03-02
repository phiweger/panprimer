from itertools import product

from panprimer.utils import orient_nodes_on_genome
from screed import rc as revcomp


def test_orient_nodes_on_genome():
    seq = '...................AAAAACCCCC.................CCTCCGGGGG.........'
    target = 'AAAAACCCCCNNNNNNNNNNNNNNNNNCCTCCGGGGG'
    
    p1 = 'AAAAACCCCC'
    p2 = 'CCCCCGGAGG'
    
    pairs = list(product([p1, revcomp(p1)], [p2, revcomp(p2)]))
    for i, j in pairs:
        assert orient_nodes_on_genome(seq, p1, p2, 17) == target