# motif-mark

```motif-mark.py``` is a simple python program that visualizes motifs in nucleotide sequences. It takes two inputs:
1. a FASTA file with read lengths $\leq$ 1000 bp (introns in lowercase, exons in uppercase)
2. a text file with the motif sequences of $\leq$ 10 bp per motif (one motif per line)


### Requirements

This program requires the ```pycairo``` package and ```python3.12```


### Usage

This program can be called from the command line using -f and -m arguments containing the location of the the FASTA and motifs files respectively. The output, a PNG image, goes to your current directory.

```
python motif-mark.py -f <path_to_fasta_file> -m <path_to_motifs_file>

```

Example Input:
```
python motif-mark.py -f ./Figure_1.fasta -m ./Fig_1_motifs.txt
```

Example Output: (better if viewed in light mode :) )

![Example Output](./Figure_1.png)




