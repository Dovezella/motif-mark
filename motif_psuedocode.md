# Motif-mark psuedocode:  


## The Problem ##


Given:  
- a file indicating desired motifs to identify
- a fasta file indicating sequences with introns (lower-case) and exons (upper-case) 

Output:
- a png file depicting each sequence to scale with the locations of found motifs 



## Data Organization ##


This should be run in an environment with python and pycairo packages installed.
>  Mine is called 'motif'. (hint:  mamba activate)

This __needs__ to be done with __object-oriented programming__. 
- Classes for:

    - parsing motif to create Motif objects (or object with all Motifs?) 
    - parsing fasta file to populate objects with sequences to identify motifs
    - creation of exons
    - creation of introns (both associated with specific sequence)
    - pycairo image creation


## Classes ##

- class Motif:
    - slots: 'seq', 'length', 'intron/exon_type'
    - def __repr__(self) -> str: 
    - methods:
        - ??

- class Gene: 
    - slots: 'name', 'seq', 'length', 'intron', 'exon'
    - methods:
        - ??

- class Fasta:
    - slots: 'file', 'Genes'
    - methods:
        - ??

- do I want one to hold all Motifs?

- class Intron:
    - slots: 'Gene_owner', 'start, 'end'
    - methods: 
        - ??

- class Exon:
    - slots: 'Gene_owner', 'start, 'end'
    - methods: 
        - ??

- class image:
    - slots: 'Genes', 'Motifs', 'colors', 'Gene-Motif_locations'
    - methods:
        - create_image
        - update_image

## Needed Functions ##

- __One-line fasta__ to make each sequence entry in fasta a single line

-  

## Psuedocode ##

read in files

create Motifs objects from file

create Gene and Fasta objects from fasta

create related exon and intron objects from Gene objects

identify positions of motifs per exon/intron per gene
- input to image object

create output png with image object









