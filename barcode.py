from Bio import SeqIO
from Bio.Seq import Seq
import re
import os

# NOTE: Josiahs code
# TODO: ask Josiah about theory
def read_barcode(st):
    # create dictionaries to turn each "bit" into an integer
    bitsdic = {'CGCA':0,
               'GTCG':1,
               'TCCT':2,
               'GAAA':3,
               'TTTA':4,
               'AAGT':5,
               'CATG':6}
    # again, special dictionary for the first bit
    bitsdic_old_a = {'TAAT':0,
                     'ATCC':1,
                     'CGCA':2,
                     'TCTC':3,
                     'GGTG':4,
                     'GCGT':5, 
                     'AAGA':6}

    # Searches for barcodes (4 bits) spaced across fixed sequences
    regex_pattern = "ctacc([AGTC]{4})caac([AGTC]{4})aggagccctcaagtca([AGTC]{4})gctc([AGTC]{4})ccatcc"

    # What happens if not all sections are found...?
    # TODO: key errors need to be handled! if barcode is wrong it will throw an
    # error
    m = re.search(regex_pattern,st,re.IGNORECASE)
    if m:
        first,*rest = m.groups()
        # NOTE: adding an X if we have an unknown barcode in a position
        # not sure if this is the best way to do it
        # but prevents keyerrors from halting
        try:
            first_bit = bitsdic_old_a[first]
        except(KeyError):
            first_bit = "x"

        rest_bits = []
        for bit in rest:
            try:
                rest_bits.append(bitsdic_old_a[bit])
            except(KeyError):
                rest_bits.append("x")
            
        barcode_array = [first_bit] + rest_bits
        return 'bc'+''.join([str(i) for i in barcode_array])
    else:
        return "barcode not found"
    
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
