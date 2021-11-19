# Comparative genomics of DIR

How conserved are the proteins near iodate reductase (IdrA) across the phylogeny of IdrA and IdrA/AioA-like proteins?

## Dependencies

The following must be installed and in your path:
- [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
- [MMSeqs](https://github.com/soedinglab/MMseqs2)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [MUSCLE](http://www.drive5.com/muscle/muscle.html)

## Methods

### 1. Download genomes
Download genomes from a list of RefSeq/GenBank accessions

### 2. Define phylogeny
Identify HMM hit in each genome
Generate phylogenetic tree of HMM hits
Define groups within the tree ("clades")

### 3. Define gene neighborhood composition
Obtain genes within +/- 10 positions of HMM hits ("gene neigborhoods")
Group proteins from gene neighborhoods by sequence similarity ("subfamilies")

## Analysis
- Which subfamilies are conserved which clades?
- What the functions of those subfamilies?



