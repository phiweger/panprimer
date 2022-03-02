## Marker workflow

Automatic primer design from genome graphs.


### Install "panprimer" library and dependencies


```bash
conda create -y -n panprimer -c bioconda python=3.8 sourmash nextflow networkx tqdm gffutils numpy primer3-py
conda activate panprimer
git clone github.com/phiweger/panprimer
cd panprimer && pip install -e .
```


### Run workflow

We either have a folder with genomes against ALL of which we want to find primers, or we have a directory for those we wish to include and those we want to exclude. Optionally (see `nextflow.config`) we can remove genomes from the "exclude" set which are very similar to those in "include".


```bash
# Activate docker daemon
# ... continued from above
cd workflow

nextflow run main.nf \
    --outdir results \
    --exclude '../data/include/*.gz'
    --include '../data/exclude/*.fasta'
```


Results look like


```bash
fwd,rev,Tm fwd,Tm rev,product,penalty,nodes
TAGACAGGAGTACGCTGTTA,CATCGCCAAGATACTGTACA,55.05,54.99,NA,0.0612,4323::4325
TGATTCAGCTTCTTCCTGAC,CGTGGCTAAGGATTCAACTA,55.08,54.89,NA,0.1948,4323::4329
```

