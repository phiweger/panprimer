import gffutils
import screed


prefix = 'PROKKA_07162020/PROKKA_07162020'
fp_annotation = f'{prefix}.gff'
fp_primer = 'results_bak_k31_n100/primers.csv'

db = gffutils.create_db(fp_annotation, ':memory:')

primers = []
with open(fp_primer, 'r') as file:
    _ = next(file)  # discard header
    for line in file:
        fwd, rev, *rest = line.strip().split(',')
        primers.append([fwd, rev])

for i in db.all_features():
    seq = i.sequence(f'{prefix}.fna')
    for fwd, rev in primers:
        if any([
            fwd in seq,
            rev in seq,
            screed.rc(fwd) in seq,
            screed.rc(rev) in seq]):
            (i.id)



!grep -A1 "AMAPBBHO_01368" PROKKA_07162020/PROKKA_07162020.faa

i.attributes['product']


# pyfaidx the protein, then search w/ mmseqs2

from collections import defaultdict


d = defaultdict(list)
fp = 'klebsiella/KPC-KP-UHL_nanopore_BK17487_whole-genome.fasta'
with screed.open(fp) as file:
    for read in file:
        for i, j in primers:
            
            p1 = re.search(i, read.sequence)
            p2r = re.search(screed.rc(j), read.sequence)

            p1r = re.search(screed.rc(i), read.sequence)
            p2 = re.search(j, read.sequence)

            if (p1 and p2r) or (p1r and p2):
                print('found', i, j)
                d[read.name].append((i, j))
            else:
                print('nope')

# TGTACAGTATCTTGGCGATG,ATAATTCTTCGTGGACGCAT,...,NA,0.2548,2428690::4253577
# TGTACAGTATCTTGGCGATG,GTAAGGTGGTTCTGACATCA,...,NA,0.3088,2428690::2453804

seq = 'GGATGTTAAAAACTTCGGACTTAGACGTATAAGTAACGGTGTACAGTATCTTGGCGATGCT'
p1 = 'TGTACAGTATCTTGGCGATG'
p2 = 'GTAAGGTGGTTCTGACATCA'



from panprimer.utils import on_same_strand


from glob import glob

for fp in glob('klebsiella/KPC-KP*.fasta'):
    with screed.open(fp) as file:
        for read in file:
            try:
                orientation = on_same_strand(read.sequence, p1, p2)
                print((read.name, orientation))
                # TODO: blastx sequence to find target
            except KeyError:
                pass




# TODO: inversions -- what happens? probably we now need rc of primer. lets
# exclude this special case for now.







