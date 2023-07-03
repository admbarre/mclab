from Bio import SeqIO
bad_read = "/Users/adrian/Documents/research/sequencing/adrian.barrera-velasquez@ucsf.edu_SEQ2154089A_061323T/B1_3_1_p85_C04.ab1"
good_read = "/Users/adrian/Documents/research/sequencing/adrian.barrera-velasquez@ucsf.edu_SEQ2154089A_061323T/B1_4_1_p85_D04.ab1"

def get_qual(abi):
    record = SeqIO.read(abi,"abi")
    qual = record.letter_annotations["phred_quality"]
    return qual

print("bad")
print(get_qual(bad_read))
print("---"*15)
print("good")
print(get_qual(good_read))
print("---"*15)
