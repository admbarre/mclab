from Bio import SeqIO
from Bio.Seq import Seq
import re
import os

# NOTE: Josiahs code
# TODO: clean this up later
def read_barcode(st):
    # create dictionaries to turn each "bit" into an integer
    bitsdic = {'CGCA':0,'GTCG':1,'TCCT':2,'GAAA':3,'TTTA':4,'AAGT':5,'CATG':6}
    # again, special dictionary for the first bit
    bitsdic_old_a = {'TAAT':0, 'ATCC':1, 'CGCA':2, 'TCTC':3, 'GGTG':4, 'GCGT':5, 'AAGA':6}
    m = re.search('ctacc([agtc]{4})CAAC([agtc]{4})aggagccctcaagtca([agtc]{4})GCTC([agtc]{4})ccatcc',st,re.IGNORECASE)
    if m:
        bcarray = [bitsdic_old_a[m.groups()[0]]]+[bitsdic[bit] for bit in m.groups()[1:]]
        return 'bc'+''.join([str(i) for i in bcarray])
    
def getseq(abifile):    
    with open(abifile,'rb') as f:
        rec = next(SeqIO.parse(f,'abi'))
    seq = Seq("".join([chr(b) for b in rec.annotations["abif_raw"]["PBAS1"]])).reverse_complement()
    return str(seq)

def get_barcode(abi):
    seq = getseq(abi)
    bcode = read_barcode(seq)
    return bcode

def readabis(root):
    abis = [f for f in os.listdir(root) if f.endswith(".ab1")]
    barcodes = [read_barcode(getseq(root+abi)) for abi in abis]
    return (abis,barcodes)
