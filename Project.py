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

#=======================================================================================================================
#Start of Joshua's work
#=======================================================================================================================

def count_dom_phenotype(genotypes):
    if len(genotypes) != 6:
        print("nope")
        exit()
    dom_chance=[1, 1, 1, 0.75, 0.5, 0] #chance any one child displays the dominant trait
    dom_kids=0
    for i in range(6):
        dom_kids+= 2*dom_chance[i]*genotypes[i] #two children per pair, multiplied by the No. of pairs
    return dom_kids

amino_freq={'F': 2, 'L': 6, 'I': 3, 'M': 1, 'V': 4, 'S': 6, 'P': 4, 'T': 4, 'A': 4, 'Y': 2, 'Stop': 3, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2, 'C': 2, 'W': 1, 'R': 6, 'G': 4}
def source_rna(protein):
    protein = protein.upper()
    possible_sorces=3 #there are three possible stops
    for i in protein:
        possible_sorces = possible_sorces*amino_freq[i] #multiply by number of possible amino acids
    return possible_sorces%1000000
