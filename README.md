# Neighborhood Analyses of Genomes in Genbank and RefSeq

This script replicates the subfamilies analysis used to ascertain the gene neighborhood around iodate reductases (IdrA). This procedure was used to demonstrate that the idrA gene is accompanied by the idrB, IdrP1, and IdrP2. The IdrP1 and IdrP2 genes are absent from the very similar AioA. This is the basis of the proposed functional difference between IdrA and AioA, which have similar primary amino acid sequences.

## Dependencies

The following must be installed and in your path:
- [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
- [MMSeqs](https://github.com/soedinglab/MMseqs2)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [MUSCLE](http://www.drive5.com/muscle/muscle.html)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [ete3](http://etetoolkit.org/)
- [HMMER](http://hmmer.org/download.html)

This script was developed and executed on Ubuntu 18.04.5 LTS 64-Bit, Intel® Core™ i5-6200U CPU @ 2.30GHz × 4, and 8 GB of Memory.

## Using this for your protein of interest

While this script is primarily made available for anyone to follow along the workflow, this script can be adapted to search for the genomic context of your gene of interest. To achieve this for your own purposes you will need:

Once you have HMMER installed, you will need to develop an HMM that searches for your protein of interest. Generally, a good way to go about it is to use a multifasta seed set for the protein motif you're interested in (e.g., [PF00384](https://pfam.xfam.org/family/PF00384#tabview=tab3)

You will want to do a simple MUSCLE alignment
```
$ muscle -in fasta.faa -out fasta.aln
```
You will then want to build your HMM
```
$ hmmbuild fasta.hmm fasta.aln
```

You should then test your HMM to set an appropriate threshold. There are many considerations in determining threshold beyond the scope of this tutorial, but one great place to start is testing your HMM against the reference protein dataset on [EMBL-EBI's HMMER Search](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch). 

Last thing you're going to want is genomes! This is covered in the methods below; however, a good start would be to take the Genbank or RefSeq accessions from the EMBL-EBI search and place them into the appropriate `refseq-accessions.txt` or `genbank-accessions.txt` lists in your `genomes` folder.


If this script turns out being useful for your research, I'd greatly appreciate the citation
Please cite: __Reyes-Umana, V., Henning, Z., Lee, K. et al. Genetic and phylogenetic analysis of dissimilatory iodate-reducing bacteria identifies potential niches across the world’s oceans. ISME J (2021). https://doi.org/10.1038/s41396-021-01034-5__

## Comparative genomics of DIR Tutorial

How conserved are the proteins near iodate reductase (IdrA) across the phylogeny of IdrA and IdrA/AioA-like proteins?

## Methods

It is recommended to run this script as a python notebook. The script can be run as an stand along script, but there will be a fair amount of editing needed beforehand. Nonetheless, if you choose to run it as a script, be sure to make all the edits provided in the instructions before hand. 
Additionally, you need to structure your directories appropriately or the code may break. The directory structure of the files provided here will be valid for the code within the code itself.

### 1. Download genomes 
The first two cells of the notebook will allow you to download genomes to your directory. Make sure you use the directory structure in the repository.
Download genomes from a list of RefSeq/GenBank accessions. This list is named either:

- refseq-accessions.txt
- genbank-accessions.txt

Good practice here is to save old accession lists with the date as a prefix, and keep the above lists as your running active lists.
If you already downloaded genomes, toggle the download option to  False!
If you add new genomes delete old genome files and re-download.
```
download = False 
```
Have a duplicate problem?
Open a terminal window in the directory with you accessions list, and deduplicate your entries.  
```
$ sort genbank-accessions.txt | uniq -u > genbank-accessions_dedup.txt
``` 

## Define phylogeny

### 2. Identify HMM hit in each genome
The next cell will then search your dataset of downloaded genomes for hits matching your HMM. There are several considerations to take into account when creating an HMM, such as gene level indicators and meaningful cutoffs. For the sake of this tutorial, the HMM provided is for a combined AioA/IdrA HMM that was used to initially train an HMM that includes all AioA/IdrA proteins. A subsequent HMM was used to specifically distinguish IdrA and AioA using the gene level indicators produced by this script. 
For additional information on this, I recommend looking at the methods in our [publication](https://www.nature.com/articles/s41396-021-01034-5)

It is important to define the inputs here, as they will determine what you are searching for and how stringent your search is:
```
hmm = './data/hmm/combined_iriA_aioA.hmm'
threshold = 640
```

It is essential that you harmonize the dictionaries to your hits, otherwise the code breaks
```
for k in list(hmmhits.keys()):
    if k not in path_to_hits:
        del hmmhits[k]
with open('iodate_reducing_genomes.txt', 'w') as fh:
    for item in hmmhits: 
        fh.write(item + '\n')
```
This step will put out the total number of hits identified during the search and comapre it to those in the path. The number should be the same. 

### 3. Generate phylogenetic tree of HMM hits
This next step generates a phylogenetic tree as a newick file. 

Inputs here are important as well. 
- __faa_ingroup__: This are the protein sequences from your HMM hits placed into a single multifasta file
- __faa_outgroup__: Depending on your analysis, you are going to want to have sequences that you have verified as part of your outgroup. Since here we want to place the arsenite oxidase AioA as our outgroup, we use a multifasta containing known arsenite oxidases
- __temp__: This is a temporary file that the code generates while going through this cell. You should delete it before running the code if the file is present.  
- __faa__: This file is generated by the code and combines the ingroup and outgroup. Delete this file if present.
```
faa_ingroup = hits
faa_outgroup = './data/tree/aioA.faa' 
temp = './data/tree/temp.faa'
faa = './data/tree/iriA-all.faa'
```
A couple things to note here:
There is a filtering option for you to use:
```
filt = 80
f_aln = './data/tree/iriA-'+ str(filt) + '.aln'
filter_gaps = ' '.join(['python3', './scripts/remove-gapped-positions.py', '-p', str(filt), '-i', aln, '-o', f_aln])
sp.call(filter_gaps, shell=True)
```
The tree alignment doesn't use this option, but you can if you want. Just make sure to change:
```
fasttree = ' '.join(['fasttree','-boot 10000', aln, ">", tree])
```
to
```
fasttree = ' '.join(['fasttree','-boot 10000', f_aln, ">", tree])
```

Also note that the `-boot 10000` option is turned on. You can change this in case you want to change the number of bootstraps. I recommend looking at the fasttree documentation. 

After running this cell, you should have files with the extensions:
- .aln
- .nwk

The subsequent cell will do some pythonic magic to add human readable names to the accession numbers.

Finally, you will get to the point where you will draw your tree.
This will use the ete3 package and draw the tree inline. 

A few considerations at this point:
- Make sure to [set your outgroup](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#tree-rooting) so that you can meaningfully infer phylogeny on your tree
- Be sure to look at the ete3 documentation. There are a lot of things you can do with this package that allows you to color clades, and change tree styles.
- The support values are drawn onto the tree as black (>0.99), dark gray (>0.90), and light gray (>0.80). I don't recommend you change these, but you can. 

You will get two trees out of the cells that draw trees. Each has its utility. The expanded tree tells you the phylogeny of each hit. The collapsed tree shows you the different IdrA clades:

![image](https://user-images.githubusercontent.com/27031932/142703922-4c6e75fa-efa8-4ad4-be38-56bea9a9d248.png)

The default settings give you some good stuff, but I recommend changing it to your needs. 

### 4. Define gene neighborhood composition
This step has two different functions, the first part will obtain genes within +/- 10 positions of HMM hits to search for "gene neigborhoods". This data is extremely useful when trying to understand the genomic context of your hits, and it's crucial for testing the stringency of your HMM model. You can expand the neighborhood size, but that is _not recommended_ doing so may be taxing on your computer and often provides you with marginally more useful information. You will be able to export these data into a csv. This data is searchable, but only useful once you've grouped your proteins into subfamilies.

The second part allows you to group proteins from gene neighborhoods by sequence similarity ("subfamilies"). 
It is imperative that the subfamilies directory is empty at this point except for the following files:
- neighborhoods.faa
- subfamilies.py

Move any other file to an archive folder for your later use and analysis.

First step is to run the cell containing:
```
seqs ="/home/victor/Downloads/dir-comparative-genomics/data/subfamilies/neighborhoods.faa"
mmseq_output = './data/subfamilies/mmseq' 
```

Next, open your shell in the `./data/subfamilies` directory, change the directory, and run the following commands in sequence:
```
$ awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' /path/to/data/subfamilies/neighborhoods.faa | awk '/^>/{f=!d[$1];d[$1]=1}f' > /path/to/data/subfamilies/neighborhoods_clean.faa
$ ./subfamilies.py neighborhoods_clean.faa
```
Running this will put out a file called `neighborhoods_clean.faa` which the `subfamilies.py` script uses to produce a folder called `neighborhoods_clean.faa_proteinClustering`.

These files will then be used by the following cell in the python notebook called __Create a readable tab separated file with protein subfamilies__ to do just that. Run the cell and it will produce:
- orf2subfamily.tsv
- orf2subfamily_clean.tsv

You can read these files, along with everything else in the proteinClustering folder, but the most useful data comes up in the next section.

## Analysis
- What are the functions of these subfamilies?
Running the first analysis cell will provide you with the answer you've been looking for all along! You will get an excel file titled `subfamily_analysis.xlsx` saved onto your ./data directory. This file is very useful and tells you all the genes near your gene of interest, what subfamily they belong to, and a bunch of additional information. You can use your favorite excel functions to filter, group, and pivot table the data into oblivion. 


- Which subfamilies are conserved which clades?

This next output is by far the most useful summary of your data, especially when paired with the information from the above analysis. You will get a graphical output, organized in the tree order from your earlier phylogenetic tree. All these pieces of data together visualizes the conservation of protein subfamilies within each organism, and the clade. Examples of the output are below (the raw output is much larger!)

### The top:
![top](https://user-images.githubusercontent.com/27031932/142739150-c50a00ee-d0ed-4688-b597-63519e41749b.png)

### Genomes with the IRI:
![DIR](https://user-images.githubusercontent.com/27031932/142739155-a0357ba1-0653-4b2b-9afb-e07c57b1d3d9.png)

### Colorbar:
![legend](https://user-images.githubusercontent.com/27031932/142739156-8c6db0f3-8069-460c-aba5-d5b1bb441862.png)

These images can be manipulated further by saving the resulting output as an SVG and displaying the data in a meaningful way.

Hopefully this tutorial was useful, and you can apply this analysis notebook to  whatever you deem necessary!

