# Problem BA2B

#Given: An integer k and a collection of strings Dna.
#Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern. (If multiple answers exist, you may return any one.)

# 1. Find the minimum distance (Mismatches)

def MinHammingDistance(p,s):
    k= len (p) #length of pattern 
    min_d= k   # Initial minimum distance (the highiest possible value, to have all of them mutated)
    for i in range(len(s) - len(p) + 1):
        Distance = sum([1 for j in range(len(p)) if p[j] != s[i:i+len(p)][j]])
        if Distance < min_d:
            min_d= Distance 
    return min_d

# 2. Generate all possible k-mers

def Possible (k):
    Bases= ['A','C','G','T']
    A= Bases
    for t in range (k-1):
        A = [i+j for i in A for j in Bases] #Create an array with all possible mismatches of length k
    return A

# 3. Find the pattern with the lowest value of d 

def MedianString (k, dna):
    pattern = Possible (k) #Call the array with all possible mismatches of length k
    d_p_dna= {} #Create a set of patterns and their No. of mismatches
    min_distance= len (pattern) * len (dna) #Assign the minimum distance to be the highest possible value 
    for i in pattern:
        Distances= 0 #Counting
        for j in range (len(dna)):
            Distances= Distances + MinHammingDistance(i,dna[j])
        d_p_dna[i]= Distances 
        if Distances < min_distance:
            min_distance = Distances #Finding the minimum distance
    for t in d_p_dna.keys():
        if d_p_dna [t]== min_distance:
            print (t) #Printing the pattern with the lowest value of d
           

print (MedianString (6, ['TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT', 'CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA', 'TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT', 'TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA', 'ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG', 'TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA', 'TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC', 'GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA' ,'CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG' ,'CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG']))



