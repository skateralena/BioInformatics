import sys

def read_file(file_name):
    '''
    Reading data from file
    '''
    F = open(file_name)
    defline = ''
    secuence = ''
    linecounter = 0
    for line in F:
        linecounter += 1
        if line[0] == '>':
            defline = line.strip()
        else:
            sequence = line.strip()
    if linecounter != 2:
        print(file_name + 'should have exactly 2 lines')
    return(sequence.lower())


def alphabet_indexes_dict(alphabet=[]):
    '''
    args:
    alphabet - list of letters that cound be in the T(text)  (list of one letter strings)
    
    return: 
    d - dict for alphabet: d['a'] returns index of 'a' in alphabet list
    alphabet - dublicated or new one if was not exist
    '''
    d = dict()
    if alphabet == []:  # if alphabet is unknown => creating dict for all English letters (using only low ones)
        ord_a = ord('a')
        for letter_ind in range(ord('a'), ord('z')+1):
            d[chr(letter_ind)] = letter_ind - ord_a
            aplhabet.append(chr(letter_ind))
    else:
        for letter_ind in range(len(alphabet)):
            d[alphabet[letter_ind]] = letter_ind 
    return d, alphabet
    

def BCR_table(P, d, alphabet):
    '''
    Creating helper table for Bad Character Rule (BCR)
    
    args:
    P - pattern(type=str)
    d - dict of alphabet letters (key - letter(str), value - index in alphabet list (int))
    alphabet - list of alphabet letters

    returns: table, where first_ind - alphabet ind, second - pattern ind
    '''
    table = [[0] * len(P) for _ in range(len(alphabet))]
    for i in range(len(alphabet)):
        for j in range(len(P)):
            if P[j] == alphabet[i]:
                table[i][j] = 0
            elif j == 0:
                table[i][j] = 1
            else:
                table[i][j] = table[i][j-1] + 1
    return table


def GSR_list(P):
    '''
    Creating list for Good Suffix rule

    arg: P - pattern
    return: list shift_gsr: shift_gsr[i] = k means if we have suffix of len i, so the shift showld be k

    helper list: gsr
        if gsr[i] = 0 means that P has no substrings (not suffix!) equal to suffix of len i (i > 0)    

        gsr[i] = k <=> k - max: k < len(P) and P[m-i::] = P[k:k+i] 
        i - len of suffix
        k - index of appering substring of P equal to suffix of len i
    complexity O(n**3)
    '''
    m = len(P)
    shift_gsr = [m] * m
    gsr = [0] * m   # !
    for i in range(1, m):
        suffix_i = P[m-i::]
        for k in range(m-i-1, -1, -1):
            if P[k:k+i] == suffix_i:
                gsr[i] = k
                shift_gsr[i] = m - k - i
                break
    if shift_gsr[1] != m and m > 1:
        for i in range(1, m): #!
            if shift_gsr[i] == m:
                shift_gsr[i] = shift_gsr[i-1] 
    return shift_gsr


def max_match_suffix(T, P, i):
    '''
    returns max len of matching suffix OR '-1' if strings equal  
    comparing T[i:i+m] with P (starting in the last element)
    '''
    m = len(P)
    suf_len = -1
    for j in range(0, m):
        if P[m-j-1] != T[i+m-j-1]:
            suf_len = j
            break
    return suf_len


'''Format of command in cmd(Windows)'''
usage = "Usage:" + sys.argv[0] + "<FASTA file> <FASTA file>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

'''main part'''
P = read_file(sys.argv[1])
T = read_file(sys.argv[2])

a = ['a', 't', 'g', 'c']  # if alphabet not given: a = []
d, alphabet = alphabet_indexes_dict(a)
table_bcr = BCR_table(P, d, alphabet)
list_gsr = GSR_list(P)
m = len(P)
locus_list = []  # list of locuses
i = 0  # index of compearing T[i:i+m] with P[i:i+m] 
while i <= len(T) - m:

    #print(T)         
    #print('-'*i + P)   # for visualisation

    suf_len = max_match_suffix(T, P, i)
    if suf_len == -1:
        locus_list.append(i+1)  # indexing should be from 1 (in program - from 0)
        i += list_gsr[1]
        continue

    shift_gsr = list_gsr[suf_len]
    ind = d[T[i + m - suf_len - 1]]  # index of mismach letter in alphabet
    shift_bcr = table_bcr[ind][m - suf_len - 1]

    if suf_len == 0:
        i += shift_bcr
    else:
        i += max(shift_gsr, shift_bcr)

print(*locus_list, sep=', ')



