def load_file(filename):
    file=open(filename,'r')
    data=file.read()
    file.close()
    dna_dict=data.split('\n')
    return dna_dict

def get_edges(dna_dict):
    adj_list = []
    for key in dna_dict:
        for n in dna_dict:
            if dna_dict[dna_dict.index(key)][-8:] == dna_dict[dna_dict.index(n)][0:8]:
                adj_list.append((dna_dict.index(key),dna_dict.index(n)))
                #print(dna_dict[key][-3:],dna_dict[n][0:3])
    return adj_list
adj_dict = {}
def get_dict(dna_dict):
    
    for key in dna_dict:
        for n in dna_dict:
            if dna_dict[dna_dict.index(key)][-8:] == dna_dict[dna_dict.index(n)][0:8]:
                adj_dict[key] = n
    return adj_dict

count = 0
dna = ''
def assemble_genome2(dna_dict):
    for str in dna_dict:
        while count <  len(dna_dict):
            if str not in adj_dict and count+1 < len(dna_dict):
                dna = ''
                count = 0
                break
            elif str not in adj_dict:
                real = dna
                break
            elif count == 0:
                dna += str
            #print(str)
            dna += adj_dict[str][8:]
            count += 1
            str = adj_dict[str]
    return real
