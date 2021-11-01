#milestone2
def find_splice(dna_motif, dna):
    position = []
    count = 0
    c = 0
    for l in dna_motif:
        for i in dna:
            if l == i:
                position.append(dna.index(i)+count)
                count = dna.index(i) + count+1
                dna = dna[dna.index(i)+1:]
                c += 1
                break
            else:
                continue
    if c != len(dna_motif):
        position = []
    return position

def func(i):
    return len(i)
def shared_motif(dna_list):
    dna_list.sort(key=func)     #sort list by length of element
    short=dna_list[0]           #shortest string
    motif=''                    #aggregator
    a=0                         #start index
    size=1                      #size of test substring
    while size<= len(short):
        test=short[a:a+size]    #creating the test substring
        if a+size > len(short):       #resets a and increments size when we go past the last letter in short
            a=0
            size+=1
            continue
        count=0
        for dna in dna_list:      # Checking that every DNA has this substring
            if test in dna:
                count+=1
        if count != len(dna_list):
            a+=1
            continue
        else:
            motif=test             #sets new longest substring, increments size
            size+=1                #we don't need to keep checking for substrings of the same length
            a=0
    #print(motif)
    return motif

def get_edges(dna_dict):
    adj_list = []
    for key in dna_dict:
        for n in dna_dict:
            if n == key:
                continue
            elif dna_dict[key][-3:] == dna_dict[n][0:3]:
                adj_list.append((key,n))
                #print(dna_dict[key][-3:],dna_dict[n][0:3])
    return adj_list

import itertools
def assemble_genome(dna_list):
    all_perms=list(itertools.permutations(dna_list))
    all_results=[]
    for perm in all_perms:
        string=perm[0]
        for i in range(len(perm)-1):
            add=perm[i+1]
            overlap=''
            bign=0
            for n in range(len(perm[0])):
                if string[-n:] == add[:n]:
                    if len(overlap) < len(add[:n]):
                        overlap=add[:n]
                        bign=n
            string+= add[bign:]
        all_results.append(string)
    return min(all_results, key=len) 

import math
def perfect_match(rna):
    rna=rna.upper()
    A_count=rna.count('A')
    C_count=rna.count('C')
    G_count=rna.count('G')
    U_count=rna.count('U')
    if U_count == A_count and C_count == G_count:
        return math.factorial(A_count) * math.factorial(C_count)
    else:
        return 0

import math
def random_genome(dna, gc_content):
    ans_list=[]
    dna=dna.upper()
    for val in gc_content:
        CG_freq=val/2
        AT_freq=(1-val)/2
        CG_count=dna.count('C')+dna.count('G')
        AT_count=dna.count('A')+dna.count('T')
        ans=(CG_freq ** CG_count) * (AT_freq ** AT_count)
        log_ans=math.log(ans,10)        #not sure if we need to round or not
        ans_list.append(log_ans)
    return ans_list 
def rev_palindrome(dna):
    return 0
