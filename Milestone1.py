#this is the main code for the dna project

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

def mendels_law(hom, het, rec):
    def comb(n,k):
        a = 1
        b = 1
        c = 1
        for i in range(1,n+1):
            a = a*i
        for i in range(1,k+1):
            b = b*i
        for i in range(1,n-k+1):
            c = c*i
        result = a/(b*c)
        return result
    p_total = 0 #initial value for total probability
    p_total += (comb(hom,2) + (3/4)*comb(het,2))
    p_total += hom*het + hom*rec
    p_total += (1/2)*het*rec
    counter = comb(hom,2) + comb(het,2) + comb(rec,2) + hom*het + hom*rec + het*rec
    return p_total/counter


#fibonacci_rabbits
def fibonacci_rabbits(n,k):
    if n == 1:
        return 1
    elif n == 2:
        return 1
    else:
        return fibonacci_rabbits(n-1,k) + k * fibonacci_rabbits(n-2,k)

def GC_content(dna_list):
    n = ' '
    n_count = 0
    for i in dna_list:
        i_count = i.count('G')+i.count('C')
        if i_count/len(i) > n_count/len(n):
            n = i
            n_count = n.count('G')+n.count('C')
    percentage = n_count/len(n)*100
    return (dna_list.index(n), percentage)

def rna2codon(rna):  
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
    protein=''
    for i in range(0,len(rna),3):
        codon= rna[i:i+3]
        if codon in genetic_code:
            protein+= (genetic_code[codon])
        else:
            continue
    return protein

def locate_substring(dna_snippet, dna):
    indexes=[]
    for i in range(0,len(dna)-len(dna_snippet)):
        if dna_snippet == dna[i:i+len(dna_snippet)]:
            indexes.append(i)
    return indexes


def  hamming_dist(dna1, dna2):
    hammingdistance = 0
    for i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            hammingdistance += 1
        else:
            continue
    return hammingdistance

#============================================================================

def splice_rna(dna, intron_list): #this has now been tested
    for intron in intron_list:
        dna=dna.replace(intron,"#")
    exons=dna.split("#")
    protein=[]
    protein_str=''
    for exon in exons:
        exon_rna= dna2rna(exon)
        protein.append(rna2codon(exon_rna))
    for part in protein:
        protein_str += part
    return protein_str

amino_freq={'F': 2, 'L': 6, 'I': 3, 'M': 1, 'V': 4, 'S': 6, 'P': 4, 'T': 4, 'A': 4, 'Y': 2, 'Stop': 3, 'H': 2, 'Q': 2, 'N': 2, 'K': 2, 'D': 2, 'E': 2, 'C': 2, 'W': 1, 'R': 6, 'G': 4}
def source_rna(protein):
    protein = protein.upper()
    possible_sorces=3 #there are three possible stops
    for i in protein:
        possible_sorces = possible_sorces*amino_freq[i] #multiply by number of possible amino acids
    return possible_sorces%1000000


def count_dom_phenotype(genotypes):
    if len(genotypes) != 6:
        print("nope")
        exit()
    dom_chance=[1, 1, 1, 0.75, 0.5, 0] #chance any one child displays the dominant trait
    dom_kids=0
    for i in range(6):
        dom_kids+= 2*dom_chance[i]*genotypes[i] #two children per pair, multiplied by the No. of pairs
    return dom_kids


