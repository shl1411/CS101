#milestone2
def find_splice(dna_motif, dna):
    position = []
    count = 0
    for l in dna_motif:
        for i in dna:
            if l == i:
                position.append(dna.index(i)+count)
                count = dna.index(i) + count+1
                dna = dna[dna.index(i)+1:]
                break
            else:
                position = []
                break
    return position

print (find_splice(dna_motif,dna))
