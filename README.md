# Motif-mark


## The Problem ##


Given:  
- a text file indicating desired motifs to identify
- a fasta file indicating sequences with introns (lower-case) and exons (upper-case) 

Output:
- a png file depicting each sequence to scale with the locations of found motifs 



## Data Organization ##


This should be run in an environment with python and pycairo packages installed.
>  Mine is called 'motif'. (hint:  mamba activate)

This __needs__ to be done with __object-oriented programming__. 
- Classes for:

    - each Motif object (to be populated in for loop reading input text file)
    - each Gene object populated from a fasta file that has all sequences for a gene on oneline (to be searched for presence of motifs)
    - creation of exons and of introns for every gene
    - groups of exons and introns for each gene to simply downstream processing
    - information gathered from exons and introns for each gene to be used to illustrate the output
    - pycairo image creation


## Classes ##

- class Motif:
    - slots: ['_seq', '_length', '_type']
    - methods:
        - create(requires:  motif from line read in from text file -m in argparse):  used to populate information after object is created

- class Gene: 
    - slots: ['name','seq','chr','gstart','gend','eistart','eiend','length','_type','_owner_length']
    - methods:
        - create(requires key as tuple from the gene dictionary created from oneline_fasta function):  used to populate information after object is created
        - splice_me:  creates exons and introns as gene objects from a gene object with type as "transcript" (straight from fasta file)

- class GeneGroupInfo:
    - slots: ['name','parts']
    - methods:
        - populate(requires:  list of exon and intron objects):  will take a list of gene_objects that are exons and introns and group them by gene  

- class ForIllustration:
    - slots: ['info','name','positions','canvassize']
    - methods: 
        - gather_info(requires: list of exons and introns which are grouped by gene held in parts(list) attribute of GeneGroupInfo object):  will populate a dictionary of all the information related to the exons and introns needed for downstream output
        - gather_positions(requires: dictionary with regular expressions to search and info must be populated already): creates dictionary to be used for drawing found motifs at respective positions of each gene

- class DrawMe:
    - slots: ['illustrationgroups','canvaswidth','canvasheight','surface','con']
    - methods: 
        - canvas(requires: the set of Gene objects with type transcript from fasta file): prepares pycairo canvas
        - draw_legend(requires: regular expression dictionary): draws legend on canvas
        - draw_genes_with_motifs:  uses positions attribute for each ForIllustration object in illustrationsgroups to draw found motifs in genes


## Needed Functions ##

- __One_line_fasta__ to make each sequence entry in fasta a single line and output a dictionary of information for each gene in file

- __ambiguity__ to populate dictionary for each motif to assign color values and the regular expression to use in searching the genes in gather_positions ForIllustration method

## Basic Workflow ##

1. read in files to create Motif and Gene objects
    - use for loop to read -m motif text file to initialize Motif object and then use .create method
    - call oneline_fasta function on the input fasta file and *must* save variable as "__genes__" for the Gene.create method to work

2. create Exon/Intron objects (gene_objects) from gene_objects
    - initialize empty set to be filled with gene_objects 
    - populate in for loop over set of gene_objects with .splice_me method

3. create GeneGroupInfo objects for each gene's exon/intron it has
    - populate list of these grouped exon and intron gene_objects for each gene in for loop with .populate method

4. create ForIllustration objects from GeneGroupInfo and identify positions of motifs per exon/intron per gene
    - call ambiguity function on set of motif_objects
    - populate list with illustration objects in for loop over list of gene_group_info objects with .gather methods

5. create output png with DrawMe object
    - initialize draw_me object with list of illustration objects
    - use DrawMe methods in order of 
        - canvas
        - draw_legend
        - draw_genes_with_motifs







