# author : Honey Bhardwaj
# date   : 1'st april 2021

from Bio import SeqIO
import sys

#this is a function to check whether input file is a fasta file or not
def checkfasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

# this is a function to translate dna sequence into protien sequence
def translate(seq):     
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    val=len(seq)
    if len(seq)%3 != 0:
    	val=len(seq)-(len(seq)%3)
    
    try:
        for i in range(0, val, 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    except:
    	print("input file does not contain DNA sequence data")
    	sys.exit(0);

    return protein

# this function is used to read and clean the dna sequence
def read_seq(inputfile):
    with open(inputfile, "r") as f:
        seq = f.read()
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    seq = seq.replace("</pre></div><div><pre>", "")
    seq = seq.replace("N", "")
    return seq

# this function is for counting the protein found
def count_protein(protein):
    seq = protein
    seq = seq = seq.replace("__", "_")
    seq = seq.replace("___", "_")
    lsp = seq.split("_")
    final=[]
    for i in lsp:
        if (len(i)>44):
            final.append(i)
    count=len(final)
    return count



# driver code

inputfile =input("Enter the name of file:")
if(checkfasta(inputfile)):
    sequence = read_seq(inputfile)[150:]
    print("total protien found with min length 44")
    protein = translate(sequence)
    numberp = count_protein(protein)
    print("The number of protein found are",numberp)
    print()
    print("number of protien found under sequence length 44-100")
    protein = translate(sequence[43:100])
    numberp = count_protein(protein)
    print("The number of protein found are",numberp)
    print()
    print("number of protien found under sequence length 100-500")
    protein = translate(sequence[99:500])
    numberp = count_protein(protein)
    print("The number of protein found are",numberp)
    print()
    print("number of protien found under sequence length 500-beyond")
    protein = translate(sequence[499:])
    numberp = count_protein(protein)
    print("The number of protein found are",numberp)
else:
	print("input file is not in .fasta format")
	sys.exit(0);