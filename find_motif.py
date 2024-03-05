#!/usr/bin/env python

from __future__ import annotations

import argparse
import re
import math
import cairo

## Be sure to run in environment with pycairo and python ##

def get_args():
    parser = argparse.ArgumentParser(description="A program to parse FASTA file looking for motifs and output a png visualization")
    parser.add_argument("-f", "--fasta", help="fasta filename", required=True, type=str)
    parser.add_argument("-m", "--motif", help="motif file name that ends in _unique.sam", required= True, type=str)
    parser.add_argument("-o", "--output", help="name of output png file", required = False, type=str)
    
    return parser.parse_args()

args=get_args()
file=args.fasta
motif=args.motif        
png=args.output

fasta_dictionary = {}   #this is important for the oneline_fasta function to output the sequences in a dictionary 
                        #keys: tuple of header info     #values: oneline sequence
def oneline_fasta(filein: str) -> dict:
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

colors = [("magenta",(0.9,0.3,0.6)),("purple",(0.7,0.3,0.9)),("cyan",(0.3,0.9,0.9)),("yellow",(1,1,0)),("orange",(0.9,0.6,0.3))]
# print(f'{colors[0]=}')
# cycle = 0
defined = {}
find = ""
def ambiguity(motifs: set) -> dict:
    '''This function will take a motif sequence with ambiguous Y and populate a dictionary of regular expressions to be used to search for all possibilites the amibiguous Y could represent.
    The dictionary will have the tuple(motif sequence, color) as key and the regex as the value.'''
    cycle = 0
    for motif in motifs:
        find = ""
        sequence = str(motif._seq).upper()
        color = colors[cycle][0]
        color_values = colors[cycle][1]
        color_key = (sequence, color, color_values)
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
    return defined



class Motif:

    __slots__ = ['_seq', '_length', '_type']

    def __init__(self):
        '''This object holds the information about motif of interest.'''
        ## Data ##
        self._seq = None
        self._length = None
        self._type = None

    def __repr__(self): 
        return f"Motif_{self._seq}"

    ## Methods ##
    def create(self, motif_from_set):
        self._length = len(motif_from_set)
        self._seq = motif_from_set
        if motif_from_set.isupper():
            self._type = "exon"
        else:
            self._type = "intron"
    

class Gene:

    __slots__ = ['name','seq','chr','start','end','length','_type','_owner_length']

    def __init__(self):
        '''This object holds the information about a gene from the FASTA file.'''
        ## Data ##
        self.name = None
        self.seq = None
        self.chr = None
        self.start = None
        self.end = None
        self.length = None
        self._type = None
        self._owner_length = None

    def __repr__(self): 
        return f"Gene_{self.name}_{self._type}_{self.start}"
    
    ## Methods ##
    def create(self, gene_dict_key):
        self.name = gene_dict_key[0]
        self.seq = genes[gene_dict_key]
        self.start = gene_dict_key[2]
        self.chr = gene_dict_key[1][3:]
        self.end = gene_dict_key[3]
        self.length = len(genes[gene_dict_key])
        self._type = "transcript"
    
    def splice_me(self):    #remember to save as a variable to retreive ouput set
        # sequence = print(self.seq)
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
            intron._type = "intron"                     ###if this is zero-based bc python, why does .start have the first nucleotide positon output as 1 and not zero?
            intron.start = int(compare.find(i)) - 1              #if the exon or intron had multiple hits it would not find the right positions potentially 
            intron.end = int(intron.start) + int(intron.length) - 1
            intron._owner_length = self.length
            introns_and_exons.append(intron)
        for e in exons:
            exon = Gene()
            exon.chr = self.chr
            exon.name = self.name
            exon.seq = e
            exon.length = len(e)
            exon._type = "exon"
            exon.start = int(compare.find(e)) - 1
            exon.end = int(exon.start) + int(exon.length) - 1
            exon._owner_length = self.length
            introns_and_exons.append(exon)
        
        return introns_and_exons

class Gene_Group_Info:

    __slots__ = ['name','parts']

    def __init__(self, name):
        self.name = str
        self.parts = list
    
    def __repr__(self):
        return f'Group:  {self.name}'
    
    ## Methods ##
    def populate(self, exon_intron_list):
        '''This should be called in a loop over a list of all Gene exon/introns where each item is a list per related Gene.
        This will separate the entire list into a single Gene_Group_info object per each Gene.'''
        self.parts = exon_intron_list
        self.name = exon_intron_list[0].name
            

class Illustration:

    __slots__ = ['info','name','positions','canvassize']

    def __init__(self):
        '''This object holds the information about where motifs are located in genes from the FASTA file, and also performs the illustration.'''
        ## Data ##
        self.info = {}
        self.name = str
        self.positions = {}
        self.canvassize = int
        
    def __repr__(self): 
        return f"{self.name}_to_illustrate"
    
    ## Methods ##   
    def gather_positions(self,defined_regex_dict: dict[str, str]):
        '''Call this function while looping over a list or set of Illustration objects to use their .info to populate the .positions for drawing.
        The output will be a dictionary of keys = Motif_name/seq and values = match starting positions in entire Gene sequence.'''

        for key,value in defined_regex_dict.items():    #populate positions dictionary with an empty list value for every motif key
            self.positions[key] = []
        for k,v in self.info.items():
            gene_info = k
            seq = v
            gene_start = gene_info[3]
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


    def gather_info(self, intron_exon_set: list[Gene]):    #output is a dictionary of keys = tuple(gene name and total nt length : values = list of lists with exon/intron gene start pos and end pos, with motif, and with motif matches start positions on entire gene length
        '''This function will populate the Illustration object's info attribute with the info from a group of gene objects.
        It should be called in a for loop to iterate over a list of Gene_Group_Info objects, which are a list of Gene Objects for the same gene. 
        The .info attribute will be a dictionary of keys(tuple[Gene_name,Gene_Length,Gene_exon/intron, Ex/it start, Ex/it end positions]
        and values = sequence of exon/intron.'''

        for exit in intron_exon_set:
            self.name = exit.name
            gene_name = exit.name
            gene_len = exit._owner_length
            gene_info = tuple([gene_name, gene_len, exit._type, f'start: {exit.start}', f'end:{exit.end}']) 
            seq = str(exit.seq)
            search = seq.upper()   
            self.info[gene_info]=search

        return 
    
    def canvas(self, gene_objects_set):
        '''If I choose to do it this way, a for loop over all Illustration objects must be done to properly account for the longest gene length.
        This will take that set of Gene transcript objects and make .canvassize the same for each Illustration object.
        Canvas size will be 20 extra to the count of the longest gene to have a buffer of 10 counts on each side of the figure.'''
        i = 0
        length = int    
        for g in gene_objects_set:
            if i == 0:
                length = g.length
                i += 1
            else:
                if g.length == length:
                    continue
                if g.length > length:
                    length = g.length
        self.canvassize = length + 20
        


   

readmotif = open(motif, "r")
motif_set = set()
for line in readmotif:
    line = line.strip()
    motif = line.strip()
    motif = Motif()
    motif.create(line)
    motif_set.add(motif)

readmotif.close()
# print(motif_set)

genes = oneline_fasta(file)
# print(genes)

gene_objects = set()    #this needs to match exactly the format of Gene Class so the .create method works
for item in genes:
    gene = Gene()
    gene.create(item)
    gene_objects.add(gene)

is_and_es = []
for object in gene_objects:
    # print(object)
    addit = object.splice_me()
    is_and_es.append(addit) 

groups = []                 #This part is important to make it easier to separate Genes objects to group related introns and exons. 
for exit in is_and_es:
    name = Gene_Group_Info(exit[0].name)
    name.populate(exit)
    groups.append(name)
    # print(type(name.parts[0]))
    # print(name.parts[0].chr)
# print(groups)
    
# print(gene_objects)
# print(is_and_es)
    
# for item in is_and_es:
#     for thing in item:
#         print(thing.seq)
#         print(thing.start)
#         print(thing.end)
#         print(thing.length)
    
expressions = ambiguity(motif_set)
# print(expressions)

draw_these = []
for glist in groups:        #groups is a list of gene_group_info objects used to create illustration objects
    draw = Illustration()
    draw.gather_info(glist.parts)
    draw_these.append(draw)
# print(draw_these[0].info)
# print(draw_these)

for gobj in draw_these:
    gobj.gather_positions(expressions)
    gobj.canvas(gene_objects)
    # print(f'{gobj.name=}{gobj.canvassize=}')

## Find longest Gene to prepare canvas ##
i = 0
width = int    
for g in gene_objects:
    if i == 0:
        width = g.length
        i += 1
    else:
        if g.length == width:
            continue
        if g.length > width:
            width = g.length
## Get height of canvas by adding 10 for every number of Gene transcript objects in gene_obj set
height = 15
for g in gene_objects:
    height += 25
# print(f'{height=}')
# print(length)

## Make the context of the canvas equal to length plus 20 for ten extra counts on each end ##

canvaswidth = width + 20 #this is confirmed to work

## begin drawing ## if it works, turn it into a function maybe??

sf = cairo.ImageSurface(cairo.FORMAT_ARGB32, canvaswidth, height)
con = cairo.Context(sf)

con.set_source_rgb(1,1,1)
con.rectangle(0,0,canvaswidth,height)
con.fill()
con.stroke()

con.set_source_rgb(0,0,0)
con.set_line_width(1)

con.move_to(10,10)
con.set_font_size(10)
con.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
con.show_text("MOTIFS LEGEND:")


legend_x= 120
con.move_to(legend_x, 10)
legend_rec_x = 110
# cycle = 0
# print(f'{colors[0][1]=}')
for motif_name,express in expressions.items():
    text = motif_name[0]
    r = float(motif_name[2][0])
    g = float(motif_name[2][1])
    b = float(motif_name[2][2])
    # print(f'{r=}{g=}{b=}')
    con.show_text(text)
    legend_x += 100
    con.rectangle(legend_rec_x,2.5,10,10)
    # con.set_source_rgb(colors[cycle][1][0],colors[cycle][1][1],colors[cycle][1][2])
    con.set_source_rgb(r,g,b)
    con.fill()
    con.stroke()
    con.set_source_rgb(0,0,0)
    con.move_to(legend_x,10)
    # cycle += 1
    legend_rec_x += 100


startoflinex = 10
startofliney = 30
# motif_recs_y = 
con.set_source_rgb(0,0,0)
rec_start_y = 27.5   #just add 10 each time
con.move_to(startoflinex,startofliney)
for gobj in draw_these:
    for key, values in gobj.info.items():
        if key[2] == "exon":
            con.set_source_rgb(0,0,0)
            exstart = key[3]
            ex_x_axis_start = int(re.search("[0-9]+",exstart).group()) + 1   
            line_len = int(key[1]) + 10
            con.line_to(line_len,startofliney)
            con.stroke()
            exend = key[4]
            gname = key[0]
            gname_y = startofliney - 5
            con.move_to(startoflinex, gname_y)
            con.show_text(gname)
            ex_x_axis_end = int(re.search("[0-9]+",exend).group()) - ex_x_axis_start
            con.rectangle(ex_x_axis_start,rec_start_y,ex_x_axis_end,5)  #(x,y of top left corner, width, height)
            con.fill()
            con.stroke()
            # rec_start_y += 25
            startofliney += 25
            # con.move_to(startoflinex,startofliney)
            for mkey,sval in gobj.positions.items():
                w = len(mkey[0])
                r = float(mkey[2][0])
                g = float(mkey[2][1])
                b = float(mkey[2][2])
                con.set_source_rgb(r,g,b)
                for xaxism in sval:
                    x = xaxism + 10
                    con.rectangle(x,rec_start_y,w,5)
                    con.fill()
                    con.stroke()
            con.move_to(startoflinex,startofliney)
            rec_start_y += 25


sf.write_to_png(f'{png}.png')
            

# print(f'{ex_x_axis_start=}')

# for x in ex_x_axis_start:
#     con.rectangle = 

# for 
