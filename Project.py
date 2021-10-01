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
