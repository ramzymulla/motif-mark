import cairo
import argparse as arp
import os
import re
import textwrap

### GLOBALS ###
BP_PER_PIXEL = 1
ROW_HEIGHT = 50
EXON_HEIGHT = ROW_HEIGHT//2
IMAGE_WIDTH = 1000
INDENT = 5
FONTSIZE = 14
HEADING_OFFSET = ROW_HEIGHT//5
ROW_OFFSET = -ROW_HEIGHT//2 
MARKDIMS = (4, EXON_HEIGHT)

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


def get_mots(mFile: str) -> tuple:
    '''
    takes motifs file and outputs dictionary of motifs and their 
    corresponding regular expressions

    Args:
        mFile (str): motifs text file   

    Returns:
        dict: dictionary of motifs of the form <motif>:<regex>
    '''

    motifs={}
    revmots = {}
    with open(mFile, 'r') as m:
        for mot in m:
            mot = mot.strip()
            revmot = revcomp(mot)
            motifs[mot] = get_reg(mot)
            revmots[revmot] = get_reg(revmot)

    return (motifs, revmots)

def get_recs(faFile,motifs,revmots) -> dict:
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
            

            if "rev" in header:
                currMots = revmots
            else:
                currMots = motifs

            recs[currGene] = Record(header,currSeq, chrom,pos,length,currGene, currMots)

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
        self.pos = pos                              # leftmost mapping position
        self.strand = "-" if "rev" in header else "+"
        self.chrom = chrom
        self.seq = seq                              # sequence
        self.length = length                        # sequence length
        self.gname = gname                          # name of gene
        self.lines = self.get_lines()
        self.motifs = self.find_motifs(motifs)      # list of motif objects
        # self.title = f'>{gname} ({length} bp, ({self.strand}) strand):'
        self.title = f'{header}'

    def __str__(self):

        words = f'{self.header}:'

        for mot in self.motifs:
            words += f'\n{mot}\t{len(mot)}'

        return words
        
                
    def get_lines(self):

        wid = (IMAGE_WIDTH-2*INDENT)*BP_PER_PIXEL

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

        for motif in self.motifs:
            for mot in self.motifs[motif]:
                ctx.set_source_rgba(*markcolors[mot.seq],0.8)
                ctx.rectangle(mot.x,mot.y+row*ROW_HEIGHT,len(mot.seq),MARKDIMS[1])
                ctx.fill()
                ctx.set_source_rgba(0,0,0,1)
                
        return row + len(self.lines)
        

    def find_motifs(self, motifs: dict) -> dict:
        '''
        _summary_

        Args:
            motifs (dict): _description_

        Returns:
            dict: _description_
        '''
        mots = {}
        for mot in motifs:
            mot_locs = [match.start() for match in re.finditer(motifs[mot],self.seq)]
            if len(mot_locs) != 0:
                mots[mot] = [Motif(mot,mot_locs[i]) for i in range(len(mot_locs))]

        return mots

class Motif():
    def __init__(self, seq: str, start: int):
        self.seq = seq
        self.length = len(seq)
        self.start = start
        self.x = start%(IMAGE_WIDTH*BP_PER_PIXEL) + INDENT
        self.y = (start//(IMAGE_WIDTH*BP_PER_PIXEL))*ROW_HEIGHT +\
              (EXON_HEIGHT-MARKDIMS[1])//2 + ROW_OFFSET
        self.color = markcolors[seq]


    def __str__(self):
        return f'''{self.seq}, {self.length}, {self.start}'''


### Main ###

args = get_args()

faFile, mFile = args.fasta, args.motifs

basename = os.path.basename(faFile).split('.')[0]

motifs,revmots = get_mots(mFile)

markcolors = {seq:color for seq,color in zip(list(motifs), COLORS[:len(motifs)])}
markcolors.update({seq:color for seq,color in zip(list(revmots), COLORS[:len(motifs)])})


recs = get_recs(faFile,motifs,revmots)


numlines = sum([len(recs[rec].lines) for rec in recs])


with cairo.ImageSurface(cairo.FORMAT_ARGB32,
                        IMAGE_WIDTH + 10*FONTSIZE, 
                        ROW_HEIGHT*(max(len(motifs),numlines)) + 5) as surface:

    # Set up context
    ctx = cairo.Context(surface)

    ctx.set_source_rgba(0,0,0,1)
    ctx.set_operator(cairo.OPERATOR_ADD)

    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(FONTSIZE)

    # Start at row #1
    row = 1

    # iterate through each record
    for i,rec in enumerate(recs): 
        row = recs[rec].draw(row,ctx)

    # draw legend
    row = 1
    x = IMAGE_WIDTH+1
    y = row*ROW_HEIGHT - HEADING_OFFSET + ROW_OFFSET 
    ctx.move_to(x,y)
    ctx.show_text("LEGEND:")
    
    for i,motif in enumerate(motifs):
        x = IMAGE_WIDTH + 15 - len(motif)
        y = (row+i)*ROW_HEIGHT//2 - EXON_HEIGHT//4 
        ctx.set_source_rgba(*markcolors[motif])
        ctx.rectangle(x,y,len(motif),EXON_HEIGHT//2)
        ctx.fill()
        ctx.set_source_rgba(0,0,0,1)

        x = IMAGE_WIDTH + 25
        y = (row+i)*ROW_HEIGHT//2 + 5
        ctx.move_to(x,y)
        ctx.show_text(f'{motif}')

    # write to file
    surface.write_to_png(f"{basename}.png")


