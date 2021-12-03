def load_file(filename):  #open and read from file
    file=open(filename,'r')
    data=file.read()
    file.close()
    dna_dict=data.split('\n')
    return dna_dict

adj_dict = {}    #set up and initialize dictionary

def get_dict(dna_dict):    #create a dictionary of adjacent DNA fragments, where the suffix of the key matches with the prefix of the corresponding value
    for key in dna_dict:
        for n in dna_dict:
            if dna_dict[dna_dict.index(key)][-8:] == dna_dict[dna_dict.index(n)][0:8]:
                adj_dict[key] = n
    return adj_dict

def assemble_genome2(dna_dict):
    adj_dict = get_dict(dna_dict)
    count = 0
    dna = ''
    real = ''
    break_out_flag = False
    for str in dna_dict:   #str is a dna fragment, loop through all the str in the list
        while count <  len(dna_dict):
            if str not in adj_dict and count+1 < len(dna_dict):    #if str does not have a corresponding pair, reset the dna string and break
                dna = ''
                count = 0
                break
            elif str not in adj_dict:    #if str does not have a corresponding pair and we have looped through all the fragments, stop loop and return dna
                break_out_flag = True
                break
            elif count == 0:    #for the first one in the dna chain, add the whole string
                dna += str
            #print(str)
            dna += adj_dict[str][8:]  #add the fragment, from the 8th index till the end to the 'dna' string
            count += 1
            str = adj_dict[str]
        if break_out_flag:
            break
    return dna
