#this is the main code for the dna project
#Names: Shaun Lee, __
dna = 'AAAACCCGGT'
# s(dna)
def s(dna):
    count = {}
    count['A'] = dna.count('A')
    count['C'] = dna.count('C')
    count['G'] = dna.count('G')
    count['T'] = dna.count('T')
    return count

#dna2rna
def dna2rna(dna):
    rna = ''
    for symbol in dna:
        if symbol == 'T':
            rna = rna + 'U'
        else:
            rna = rna + symbol
    return rna
#reverse complement
def reverse_complement(dna):
    rna = ''
    dna = dna[::-1]
    for symbol in dna:
        if symbol == 'A':
            rna = rna + 'T'
        elif symbol == 'C':
            rna = rna + 'G'
        elif symbol == 'G':
            rna = rna + 'C'
        elif symbol == 'T':
            rna = rna + 'A'
    return rna
print (reverse_complement(dna))

#rna2codon
def rna2codon(triplet):
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    allowed_codons = set('ACGU')
    if triplet in genetic_code:
        return (genetic_code[triplet])
    else:
        return ("Invalid")

def locate_substring(dna_snippet, dna):
    indexes=[]
    for i in range(len(dna_snippet)):
        if dna_snippet.startswith(dna, i):
            indexes.append(i)
    print(indexes)
    
def  hamming_dist(dna1, dna2):
    hammingdistance = 0
    for i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            hammingdistance += 1
        else:
            continue
    return hammingdistance
