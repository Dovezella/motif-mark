#!/usr/bin/env python

## Be sure to run in environment with pycairo and python ##
## Imports ##

from __future__ import annotations

import argparse
import re
import math
import cairo


def get_args():
    parser = argparse.ArgumentParser(description="A program to parse FASTA file looking for motifs and output a png visualization")
    parser.add_argument("-f", "--fasta", help="fasta filename", required=True, type=str)
    parser.add_argument("-m", "--motif", help="motif file name that ends in _unique.sam", required= True, type=str)

    return parser.parse_args()

args=get_args()
file=args.fasta
motif=args.motif
outfilename = re.match("^[^.]*",file).group()

## Functions ##

fasta_dictionary = {}   #this is important for the oneline_fasta function to output the sequences in a dictionary 
                        #keys: tuple of header info     #values: oneline sequence
def oneline_fasta(filein: str) -> dict[tuple,str]:
    '''This will convert a multiline fasta file to a dictionary (sequences are all on one line wihtout new line characters).
    The key is the information from the header line, but this may need to be modified,
    depending on the information you want from the header.
    Don't forget to specify filein.'''
    with open (filein, "r") as fin:
        header = re.split(" |:|>|-",fin.readline().strip())
        header = tuple(header[1:5])
        nt: str = " "
        for line in fin: 
                if line[0] == ">":
                    fasta_dictionary[header] = nt
                    nt: str = " "
                    header = re.split(" |:|>|-",line.strip())
                    header = tuple(header[1:5])
                else:
                    nt += line.strip()
        fasta_dictionary[header] = nt            
    return fasta_dictionary

colors = [("green",(0.2,0.9,0.6)),("blue",(0,0,1)),("yellow",(1,0.75,0)),("purple",(0.9,0.3,0.9)),("red",(0.95,0.15,0.1))]
defined = {}    #this empty dictionary will be used in combination with the colors list of rgb values above to populate the ouput of the ambiguity function
def ambiguity(motifs: set) -> dict[tuple,str]:
    '''This function will take a motif sequence and populate a dictionary of regular expressions to be used to search for all possibilites the amibiguities could represent.
    The dictionary will have the tuple(str[motif sequence], tuple[color,tuple[rgbvalues]]) as key and the str[regex] as the value.'''
    cycle = 0   #this will control the motif chosen, by increasing the index of the set of motifs by one as each regular expression is finished for a motif
    for motif in motifs:
        find = ""   #this will hold the regex being created
        sequence = str(motif._seq).upper()  #ensure that the motif becomes upper case to properly match the gene sequence which will also be uppercased
        color = colors[cycle][0]
        color_values = colors[cycle][1]
        color_key = (sequence, color, color_values) #assign each motif a color with its rgb values in a tuple 
        for letter in sequence:
            if letter == "Y":
                find += "[CT]"
            elif letter == "U":
                find += "[T]"
            elif letter == "R":
                find += "[AG]"
            elif letter == "K":
                find += "[GT]"
            elif letter == "M":
                find += "[AC]"
            elif letter == "S":
                find += "[CG]"
            elif letter == "W":
                find += "[AT]"
            elif letter == "B":
                find += "[CGT]"
            elif letter == "D":
                find += "[TAG]"
            elif letter == "H":
                find += "[TAC]"
            elif letter == "V":
                find += "[ACG]"
            elif letter == "N":
                find += "[TACG]"    
            else:
                find += f'[{letter}]'
            
        defined[color_key] = find
        cycle += 1
    return defined  #use this defined motif dictionary when searching for motifs in sequences 

## Classes ##

class Motif:

    __slots__ = ['_seq', '_length', '_type']

    def __init__(self):
        '''This object holds the information about motif of interest.'''
        ## Data ##
        self._seq = None
        self._length = None
        self._type = None

    def __repr__(self): 
        return f"motif_{self._seq}"

    ## Methods ##
    def create(self, motif_line_from_file: str):       #call this function in a for loop to iterate over input motif text file 
        self._length = len(motif_line_from_file)
        self._seq = motif_line_from_file
        if motif_line_from_file.isupper():
            self._type = "exon"
        else:
            self._type = "intron"
    

class Gene:

    __slots__ = ['name','seq','chr','gstart','gend','eistart','eiend','length','_type','_owner_length']

    def __init__(self):
        '''This object holds the information about a gene from the FASTA file.'''
        ## Data ##
        self.name = None
        self.seq = None
        self.chr = None
        self.gstart = None
        self.gend = None
        self.eistart = None
        self.eiend = None
        self.length = None
        self._type = None
        self._owner_length = None

    def __repr__(self):
        if self._type == "transcript":
            return f"gene_{self.name}_{self._type}_{self.gstart}"
        else:
            return f"gene_{self.name}_{self._type}_{self.eistart}"
    
    ## Methods ##
    def create(self, gene_dict_key: tuple):        #use the fasta_dictionary output from the oneline_fasta function which must be called genes as input in a for loop for each item in the dict
        self.name = gene_dict_key[0]                            #^be sure to store output in a set/ or modify script if using list or other structure
        self.seq = genes[gene_dict_key]
        self.gstart = gene_dict_key[2]
        self.chr = gene_dict_key[1][3:]
        self.gend = gene_dict_key[3]
        self.length = len(genes[gene_dict_key])         
        self._type = "transcript"                         
    
    def splice_me(self):    #remember to save as a variable to retreive ouput set
        '''This function will be used to populate a list of the introns and exons as Gene Class objects for a specific gene_object with ._type = transcript.
        This should be called in a for loop over the set of gene objects, and the output should be saved in a list which was initialized prior to for loop.
        These can then be separated into GeneGroupInfo Class objects.'''
        introns_and_exons = []
        introns = re.findall("[tuagc]+",self.seq)
        exons = re.findall("[TUAGC]+",self.seq)     #even though these are underlined, it is confirmed to be working
        compare = str(self.seq)
        for i in introns:
            intron = Gene()
            intron.chr = self.chr
            intron.name = self.name
            intron.seq = i
            intron.length = len(i)
            intron._type = "intron"
            intron.gstart = self.gstart
            intron.gend = self.gend      
            intron.eistart = int(compare.find(i)) - 1              
            intron.eiend = int(intron.eistart) + int(intron.length) - 1
            intron._owner_length = self.length
            introns_and_exons.append(intron)
        for e in exons:
            exon = Gene()
            exon.chr = self.chr
            exon.name = self.name
            exon.seq = e
            exon.length = len(e)
            exon._type = "exon"
            exon.gstart = self.gstart
            exon.gend = self.gend 
            exon.eistart = int(compare.find(e)) - 1
            exon.eiend = int(exon.eistart) + int(exon.length) - 1
            exon._owner_length = self.length
            introns_and_exons.append(exon)
        
        return introns_and_exons

class GeneGroupInfo:

    __slots__ = ['name','parts']

    def __init__(self):
        '''This class will hold the intron-exon lists and all their info from a specific gene. Info will be passed to an Illustration object.'''
        self.name = str
        self.parts = list
    
    def __repr__(self):
        return f'group_{self.name}'
    
    ## Methods ##
    def populate(self, exon_intron_list: list):
        '''This should be called in a loop over a list of all Gene exon/introns where each item is a list per related Gene.
        This will separate the entire list into a single Gene_Group_info object per each Gene.'''
        self.parts = exon_intron_list
        self.name = exon_intron_list[0].name
            

class ForIllustration:

    __slots__ = ['info','name','positions','canvassize']

    def __init__(self):
        '''This object holds the information about where motifs are located in genes from the FASTA file, and also performs the illustration.'''
        ## Data ##
        self.info = {}  #see function gather_info for more details
        self.name = str
        self.positions = {} #see function gather_positions for more details
        # self.canvassize = int
        
    def __repr__(self): 
        return f"{self.name}_to_illustrate"
    
    ## Methods ##   
    def gather_positions(self,defined_regex_dict: dict[tuple, str]):
        '''Call this function while looping over a list or set of Illustration objects to use their self.info to populate the self.positions for drawing.
        The output will be a dictionary of keys = str[Motif_name/seq] and values = list[match starting positions in entire Gene sequence].'''

        for key,value in defined_regex_dict.items():    #populate positions dictionary with an empty list value for every motif key
            self.positions[key] = []
        for k,v in self.info.items():
            gene_info = k
            seq = v
            gene_start = gene_info[6]
            gene_start_pos = int(re.search("[0-9]+",gene_start).group())
            for key,value in defined_regex_dict.items():
                find_me = f"(?={value})"
                search = seq
                starts = list(re.finditer(find_me,search))
                if starts == []:
                    continue
                else:
                    for pos in starts:
                        motif_starts_here = int(pos.start())
                        corrected = gene_start_pos + motif_starts_here
                        self.positions[key].append(corrected)
        return


    def gather_info(self, intron_exon_list: list[Gene]):
        '''This function will populate the Illustration object's info attribute with the info from a group of gene objects.
        It should be called in a for loop to iterate over a list of Gene_Group_Info objects, which are a list of Gene Objects for the same gene. 
        The self.info attribute will be a dictionary of keys(tuple[Gene_name,Gene_Length,Gene_exon/intron, Chromosome, gene start on chromosome, gene end on chromosome, Ex/it start, Ex/it end positions]
        and values = sequence of exon/intron.'''

        for exit in intron_exon_list:
            self.name = exit.name
            gene_name = exit.name
            gene_len = exit._owner_length
            gchr = exit.chr
            gene_info = tuple([gene_name, gene_len, exit._type, gchr, exit.gstart, exit.gend, f'start: {exit.eistart}', f'end:{exit.eiend}']) 
            seq = str(exit.seq)
            search = seq.upper()   
            self.info[gene_info]=search

        return 

class DrawMe:

    __slots__ = ['illustrationgroups','canvaswidth','canvasheight','surface','con']

    def __init__(self, illustrationgroups: list):
        '''This object will take the list of Illustration objects and conduct pycairo to draw the desired output.'''
        ## Data ##
        self.illustrationgroups = illustrationgroups
        self.canvaswidth = int
        self.canvasheight = int
        self.surface = None
        self.con = None #surface and context must be populated with canvas function
        
    def __repr__(self): 
        return f"draw_motifs_for_genes"
    
    ## Methods ##  

    def canvas(self, set_of_gene_transcript_objects: set):
        '''This function will create the canvas and background for the image.'''
        i = 0
        width = int   
        height = 90     #ensures room for legend at the top of the image            
        for g in set_of_gene_transcript_objects:
            height += 40 
            if i == 0:
                width = g.length
                i += 1
            else:
                if g.length == width:
                    continue
                if g.length > width:
                    width = g.length
        self.canvaswidth = width + 20 #20 gives left and right margins
        self.canvasheight = height
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.canvaswidth, self.canvasheight)
        self.con = cairo.Context(self.surface)

        self.con.set_source_rgb(1,1,1)
        self.con.rectangle(0,0,self.canvaswidth,self.canvasheight)   #make the background white
        self.con.fill()
        self.con.stroke()
    
    def draw_legend(self, defined_regex_dict: dict[tuple, str]):
        self.con.set_source_rgb(0,0,0)
        self.con.move_to(10,40)
        self.con.set_font_size(10)
        self.con.select_font_face("HelveticaNeueLT Std Lt", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)   #Arial
        self.con.show_text("MOTIFS LEGEND:")
        self.con.move_to(10,70)
        self.con.show_text("GENES:")

        legend_x= 120
        self.con.move_to(legend_x, 40)
        legend_rec_x = 105
        for motif_name,express in defined_regex_dict.items():
            text = motif_name[0]
            r = float(motif_name[2][0])
            g = float(motif_name[2][1])
            b = float(motif_name[2][2])
            self.con.show_text(text)
            legend_x += 100
            self.con.rectangle(legend_rec_x,32,10,10)
            self.con.set_source_rgb(r,g,b)
            self.con.fill()
            self.con.stroke()
            self.con.set_source_rgb(0,0,0)
            self.con.move_to(legend_x,40)
            legend_rec_x += 100
        
    def draw_genes_with_motifs(self):
        self.con.set_source_rgb(0,0,0)   #set color to black 
        self.con.set_line_width(2)
        startoflinex = 10
        startofliney = 100
        self.con.set_source_rgb(0,0,0)
        rec_start_y = 95   #just add 10 each time
        self.con.move_to(startoflinex,startofliney)
        for gobj in self.illustrationgroups:
            for key, values in gobj.info.items():
                if key[2] == "exon":
                    # print(f'{key=}')
                    self.con.set_source_rgb(0,0,0)
                    exstart = key[6]
                    ex_x_axis_start = int(re.search("[0-9]+",exstart).group()) + 11   #add one to move from 0-based numbers to matrix-math based numbers and 10 for margin space
                    line_len = int(key[1]) + 10
                    self.con.line_to(line_len,startofliney)
                    self.con.stroke()
                    exend = key[7]
                    gname = key[0]
                    chrname = key[3]
                    chrstart = key[4]
                    chrend = key[5]
                    gtext = f'{gname} chr{chrname}:{chrstart}-{chrend}'
                    gname_y = startofliney - 10
                    self.con.move_to(startoflinex, gname_y)
                    self.con.show_text(gtext)
                    ex_x_axis_end = (int(re.search("[0-9]+",exend).group())+11) - ex_x_axis_start
                    self.con.rectangle(ex_x_axis_start,rec_start_y,ex_x_axis_end,10)  #(x,y of top left corner, width, height)
                    self.con.fill()
                    self.con.stroke()
                    startofliney += 40
                    for mkey,sval in gobj.positions.items():
                        w = len(mkey[0])
                        # print(f'{w=},{mkey=}')
                        r = float(mkey[2][0])
                        g = float(mkey[2][1])
                        b = float(mkey[2][2])
                        self.con.set_source_rgb(r,g,b)
                        for xaxism in sval:
                            x = xaxism + 11
                            self.con.rectangle(x,rec_start_y,w,10)
                            self.con.fill()
                            self.con.stroke()
                    self.con.move_to(startoflinex,startofliney)
                    rec_start_y += 40
        self.surface.write_to_png(f'{outfilename}.png')

## Now the script to achieve the goal of producing a png output based on the input Fasta file and motif text file for where each motif is found in each fasta sequence ## 
   
readmotif = open(motif, "r")    #begin by reading in the motifs
motif_set = set()               #initialize empty set to hold Motif objects
for line in readmotif:
    line = line.strip()
    motif = line.strip()
    motif = Motif()             #assign the current line as Motif class
    motif.create(line)          #populate attributes
    motif_set.add(motif)        #add to set

readmotif.close()               

genes = oneline_fasta(file)     #populate dictionary of fasta information for each sequence/gene, it MUST be called genes

gene_objects = set()    #create empty set to populate with gene objects created by for loop over genes: dict from fasta file
for item in genes:
    gene = Gene()
    gene.create(item)       #populate attributes
    gene_objects.add(gene)  #add to set of gene_objects 

is_and_es = []      #initialize an empty list to hold lists of spliced genes (introns/exons as gene objects) for each gene transcript object
for object in gene_objects:
    addit = object.splice_me()  #output of splice_me is a list of introns and exons for one gene
    is_and_es.append(addit) 

groups = []                 #This part is important to make it easier to separate Genes objects to group related introns and exons. 
for exit in is_and_es:
    name = GeneGroupInfo()  
    name.populate(exit)
    groups.append(name)

expressions = ambiguity(motif_set)  #populate a dictionary of regular expressions to be used to search for each motif in a sequence

draw_these = []         #this for loop populates a list of illustration objects for each gene
for glist in groups:        #groups is a list of gene_group_info objects used to create illustration objects
    draw = ForIllustration()
    draw.gather_info(glist.parts)
    draw.gather_positions(expressions)
    draw_these.append(draw) 


## Now to draw ##
to_draw = DrawMe(draw_these)
to_draw.canvas(gene_objects)
to_draw.draw_legend(expressions)
to_draw.draw_genes_with_motifs()    
           