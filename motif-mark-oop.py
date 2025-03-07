import cairo
import argparse as arp
import os
import re
import textwrap

### GLOBALS ###
BP_PER_PIXEL = 1                    # how many bps in each horizontal pixel
ROW_HEIGHT = 50                     # how many pixels in each row
EXON_HEIGHT = ROW_HEIGHT//2         # how tall to draw the exon boxes
MOTIF_STAGGER = 1                   # how much to stagger the motifs
MARK_HEIGHT = EXON_HEIGHT - \
                MOTIF_STAGGER*2     # how tall to draw the marks
MAX_WIDTH = 1000                    # maximum width (in pixels) to draw each record
INDENT = 5                          # how many pixels (from left) to indent the image
FONT = "Arial"                      # text font
FONTSIZE = 14                       # text font size
HEADING_OFFSET = ROW_HEIGHT//5      # vertical offset for row headings
ROW_OFFSET = -ROW_HEIGHT//2         # vertical offset for row contents

COLORS = [
    (220/255, 20/255, 60/255),      # Crimson
    (34/255, 139/255, 34/255),      # Forest Green
    (0, 0, 139/255),                # Dark Blue
    (255/255, 223/255, 0),          # Gold
    (0, 255/255, 255/255),          # Cyan
    (255/255, 0, 255/255),          # Magenta
    (255/255, 165/255, 0),          # Orange
    (128/255, 0, 128/255),          # Purple
    (0, 128/255, 0),                # Dark Green
    (186/255, 85/255, 211/255),     # Medium Orchid
    (255/255, 105/255, 180/255),    # Hot Pink
    (255/255, 99/255, 71/255),      # Tomato Red
    (0, 128/255, 128/255),          # Teal
    (255/255, 192/255, 203/255),    # Light Pink
    (255/255, 255/255, 0),          # Yellow
    (240/255, 128/255, 128/255),    # Light Coral
    (0, 0, 255/255),                # Blue
    (0, 255/255, 0),                # Green
    (255/255, 0, 0),                # Red
    (255/255, 228/255, 181/255),    # Moccasin
]


DEGEN_BASES = {
    'A': '[A]',     'a': '[a]',   
    'T': '[T]',     't': '[t]',
    'C': '[C]',     'c': '[c]',
    'G': '[G]',     'g': '[g]',
    'U': '[U]',     'u': '[u]',
    'R': '[AG]',    'r': '[ag]',       # Purine
    'Y': '[CT]',    'y': '[ct]',       # Pyrimidine
    'S': '[GC]',    's': '[gc]',       # Strong 
    'W': '[AT]',    'w': '[at]',       # Weak
    'K': '[GT]',    'k': '[gt]',       # Keto
    'M': '[AC]',    'm': '[ac]',       # Amino
    'B': '[CGT]',   'b': '[cgt]',      # not A
    'D': '[AGT]',   'd': '[agt]',      # not C
    'H': '[ACT]',   'h': '[act]',      # not G
    'V': '[ACG]',   'v': '[acg]',      # not T
    'N': '[ACGT]',  'n': '[acgt]'      # Any one base
}

DEGEN_COMPS = {
    'A': 'T',   'a': 't',   # Adenine -> Thymine
    'T': 'A',   't': 'A',   # Thymine -> Adenine
    'C': 'G',   'c': 'g',   # Cytosine -> Guanine
    'G': 'C',   'g': 'c',   # Guanine -> Cytosine
    'R': 'Y',   'r': 'y',   # Purine (A or G) -> Pyrimidine (C or T)
    'Y': 'R',   'y': 'r',   # Pyrimidine (C or T) -> Purine (A or G)
    'S': 'S',   's': 's',   # G or C -> G or C
    'W': 'W',   'w': 'w',   # A or T -> A or T
    'K': 'M',   'k': 'm',   # G or T -> A or C
    'M': 'K',   'm': 'k',   # A or C -> G or T
    'B': 'V',   'b': 'v',   # C or G or T -> A or G or C
    'D': 'H',   'd': 'h',   # A or G or T -> A or C or T
    'H': 'D',   'h': 'd',   # A or C or T -> A or G or T
    'V': 'B',   'v': 'b',   # C or G or A -> C or G or T
    'N': 'N',   'n': 'n'    # Any base -> Any base
}

DEGEN_TAB = str.maketrans(DEGEN_COMPS)


def get_args():
    '''Parses user inputs from command line'''
    parser = arp.ArgumentParser(description='''###
                                This is a python program for visualizing motifs in sequencing data. 
                                It accepts two files: a fasta file containing the sequencing data, with introns
                                and exons designated by lower and upper case respectively, and a text file containing the 
                                motifs that you would like to mark (one motif per line). Each sequence in the fasta
                                file must be <= 1000 bp, and each motif must be <= 10 bp.
                                ###''') 
    parser.add_argument("-f", "--fasta", help="path to input fasta file",
                        required=True)
    parser.add_argument("-m", "--motifs", help="path to input motifs file",
                        required=True)
    return parser.parse_args()

def revcomp(seq: str) -> str:
    return seq[::-1].translate(DEGEN_TAB)

def get_reg(seq: str) -> str:
    '''
    Takes string input and converts it into a regular expression

    Args:
        seq (str): sequence of characters

    Returns:
        str: regular expression
    '''

    reg = ""

    for char in seq:
        reg += "+" + DEGEN_BASES[char]

    return reg[1:]


def get_mots(mFile: str) -> dict:
    '''
    takes motifs file and outputs dictionary of motifs and their 
    corresponding regular expressions

    Args:
        mFile (str): motifs text file   

    Returns:
        dict: dictionary of motifs of the form <motif>:<regex>
    '''

    motifs={}
    with open(mFile, 'r') as m:
        for mot in m:
            mot = mot.strip()
            motifs[mot] = get_reg(mot)

    return motifs

def get_recs(faFile: str, motifs: dict) -> dict:
    '''
    generates dictionary of Record objects from a fasta file and motif dictionaries

    Args:
        faFile (str): destination of fasta file
        motifs (dict): dictionary of motifs and their regex patterns

    Returns:
        dict: Records dictionary with gene names as the keys
    '''
    currGene=''
    currSeq=''
    recs = {}
    with open(faFile, 'r') as f:

        line = f.readline().strip()

        while line:
            header = line
            line = line[1:].split()
            currGene = line[0]
            loc = line[1]
            
            line = f.readline().strip()

            while line[0] != ">":
                currSeq += line

                line = f.readline().strip()

                if not line: 
                    break

            chrom,pos = loc.split(":")
            length = len(currSeq)
            

            recs[currGene] = Record(header,currSeq, chrom,pos,length,currGene, motifs)

            # print(recs[currGene])
            currSeq = ''

    return recs

class Record():
    def __init__(self, 
                 header: str,
                 seq: str, 
                 chrom: str,
                 pos: str, 
                 length: int, 
                 gname: str,
                 motifs: dict) -> None:
    
        self.header = header
        self.pos = pos                              # mapping position
        self.chrom = chrom
        self.seq = seq                              # sequence
        self.length = length                        # sequence length
        self.gname = gname                          # name of gene
        self.lines = self.get_lines()
        self.motifs = self.find_motifs(motifs)      # list of motif objects
        self.title = f'{gname} ({chrom}:{pos}):'
        # self.title = f'{header}'

    def __str__(self):

        words = f'{self.header}:'

        for mot in self.motifs:
            words += f'\n{mot}\t{len(mot)}'

        return words
        
                
    def get_lines(self) -> list:
        '''
        splits record into lines based on drawn image width, then splits each line
        into intron/exon segments

        Returns:
            list: list of tuples containing intron/exon segments for each line in the record
        '''

        wid = (MAX_WIDTH-2*INDENT)*BP_PER_PIXEL

        lines = textwrap.wrap(self.seq, width=wid)
        seglines = []
        for line in lines:
            segments = re.split(r"([a-z]+)" , line)

            if line[0].islower():
                segments = segments[1:]
            
            if line[-1].islower():
                segments = segments[:len(segments)-1]

            seglines.append(tuple(segments))
            
        return seglines

    def draw(self, row: int, ctx):
        '''
        Draws record to an existing context object, returns the new row value

        Args:
            row (int): counter for current row in image
            ctx (_type_): cairo.Context object

        Returns:
            new row value after drawing record
        '''

        # Write label
        x = INDENT
        y = row*ROW_HEIGHT - HEADING_OFFSET + ROW_OFFSET 
        ctx.move_to(x,y)
        ctx.show_text(self.title)

        for i,line in enumerate(self.lines):
            x2 = INDENT

            for seg in line:
                if seg[0].islower():
                    x1 = x2
                    x2 += len(seg)
                    y = (row+i)*ROW_HEIGHT + EXON_HEIGHT//2 + ROW_OFFSET

                    # draw line
                    ctx.move_to(x1,y)
                    ctx.line_to(x2,y)
                    ctx.stroke()

                else:
                    # draw rectangle
                    x1 = x2
                    x2 += len(seg)
                    y = (row+i)*ROW_HEIGHT + ROW_OFFSET
                    ctx.rectangle(x1,y,len(seg),EXON_HEIGHT)
        
        # draw motifs
        for i,mot in enumerate(self.motifs):
            ctx.set_source_rgba(*markcolors[mot.seq],0.8)
            stag = MOTIF_STAGGER if  i%2==0 else -MOTIF_STAGGER
            ctx.rectangle(mot.x,mot.y+row*ROW_HEIGHT + stag,len(mot.seq),MARK_HEIGHT)
            ctx.fill()
            ctx.set_source_rgba(0,0,0,1)
                
        return row + len(self.lines)
        

    def find_motifs(self, motifs: dict) -> list:
        '''
        finds loci of each motif in the record

        Args:
            motifs (dict): regex dictionary for each motif

        Returns:
            list: Motif objects sorted by starting position
        '''
        mots = []
        for mot in motifs:
            mot_locs = [match.start() for match in re.finditer(motifs[mot],self.seq)]
            if len(mot_locs) != 0:
                mots += [Motif(mot,mot_locs[i]) for i in range(len(mot_locs))]
       
        # return motifs sorted by starting position
        return sorted(mots,key = lambda x: x.start) 

class Motif():
    def __init__(self, seq: str, start: int):
        self.seq = seq
        self.length = len(seq)
        self.start = start
        self.x = start%(MAX_WIDTH*BP_PER_PIXEL) + INDENT
        self.y = (start//(MAX_WIDTH*BP_PER_PIXEL))*ROW_HEIGHT +\
              (EXON_HEIGHT-MARK_HEIGHT)//2 + ROW_OFFSET
        self.color = markcolors[seq]


    def __str__(self):
        return f'''{self.seq}, {self.length}, {self.start}'''


### Main ###

args = get_args()

faFile, mFile = args.fasta, args.motifs

basename = os.path.basename(faFile).split('.')[0]

motifs = get_mots(mFile)

# designate colors for each motif 
markcolors = {seq:color for seq,color in zip(list(motifs), COLORS[:len(motifs)])}

# generate Records
recs = get_recs(faFile,motifs)


# calculate total number of lines to draw
numlines = sum([len(recs[rec].lines) for rec in recs])
recwid = min(max([recs[rec].length for rec in recs]),MAX_WIDTH)


# start drawing
with cairo.ImageSurface(cairo.FORMAT_ARGB32,
                        recwid + 10*FONTSIZE + 2*INDENT, 
                        ROW_HEIGHT*(max(len(motifs),numlines)) + 5) as surface:

    # Set up context
    ctx = cairo.Context(surface)

    ctx.set_source_rgba(0,0,0,1)

    # for blending colors in case of motif overlap
    ctx.set_operator(cairo.OPERATOR_ADD)

    ctx.select_font_face(FONT, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(FONTSIZE)

    # Start at row #1
    row = 1

    # iterate through each record
    for i,rec in enumerate(recs): 
        row = recs[rec].draw(row,ctx)

    # draw legend
    row = 1

    # show legend title text
    x = recwid 
    y = row*ROW_HEIGHT - HEADING_OFFSET + ROW_OFFSET 
    ctx.move_to(x,y)
    ctx.show_text("LEGEND:")
    
    # iterate through motifs
    for i,motif in enumerate(motifs):
        # set x,y coords
        x = recwid + 5*INDENT - len(motif)
        y = (row+i)*ROW_HEIGHT//2 - EXON_HEIGHT//4 

        # set color and width based on the motif
        ctx.set_source_rgba(*markcolors[motif])
        ctx.rectangle(x,y,len(motif),EXON_HEIGHT//2)

        # draw motif box
        ctx.fill()
        ctx.set_source_rgba(0,0,0,1)

        # show motif sequence
        x = recwid + 7*INDENT
        y = (row+i)*ROW_HEIGHT//2 + 5
        ctx.move_to(x,y)
        ctx.show_text(f'{motif}')

    # write to file
    surface.write_to_png(f"{basename}.png")


