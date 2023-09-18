#Problem BA2D
#Given: Integers k and t, followed by a collection of strings Dna.
#Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t). If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

def CountMatrix(Dna):
    Dna_Matrix = [[Dna[j][i]for j in range (len(Dna))]for i in range (len(Dna[0]))]
    count_A=0
    count_C=0
    count_G=0
    count_T=0
    list =[]
    for row in Dna_Matrix:
        for base in row:
            if base == 'A':
                count_A+=1
            if base == 'C':
                count_C+=1
            if base == 'G':
                count_G+=1
            if base == 'T':
                count_T+=1
        list.append(str(count_A)+ str(count_C)+str(count_G)+str(count_T) )
        count_A=0
        count_C=0
        count_G=0
        count_T=0
        
    count_Matrix= [[list[j][i]for j in range (len(list))]for i in range(len(list[0]))]
    
    
    return count_Matrix

def ProfileMatrix(count_Matrix):
    list = []
    total=0
    n= len(count_Matrix)
    a= [int(item [0]) for item in count_Matrix]
    total= sum(a)
         
    profile_Matrix= [[int(x)/total for x in lst]for lst in count_Matrix]
    
    return profile_Matrix

def KmerProbability(kmer,profile_Matrix):
    
    probability= 1
    count = 0
    
    for i in kmer:
        
        global nuc
        if i=='A':
            nuc= profile_Matrix[0][count]
        if i=='C':
            nuc= profile_Matrix[1][count]
        if i=='G':
            nuc= profile_Matrix[2][count]
        if i=='T':
            nuc= profile_Matrix[3][count]
        count +=1
        
        probability= probability*nuc
        
    return probability

def SelectFirstKmer (k, Dna):
    list= []
    for i in Dna:
        list.append(i[0:k])
    return list

def AllPossibleKmers(Dna, k):
   
    kmer_list=[]
    for dna in Dna:
        for x in range (0, len(dna)-k-1):
            for i in range(len(dna)-k+1):
                row= dna[i:i+k]
                kmer_list.append(row) 
    kmer_list= list(dict.fromkeys(kmer_list))
            
    return kmer_list

def MostProbableK_mer(AllPossiblekmers, profile_Matrix):
    
    max_probability=0
    
    dict= {}
    
    for i in AllPossiblekmers:
        probability= KmerProbability(i,profile_Matrix)
        dict[probability]=i
        
        if probability> max_probability:
            max_probability = probability
    
    return dict[max_probability]

def Score(Motifs):
    score = 0
    for i in range(len(Motifs[0])):
        j = [motif[i] for motif in Motifs]
        score += (len(j) - max(j.count("A"), j.count("C"), j.count("T"), j.count("G")))
    return score

def GreedyMotifSearch(Dna, k, t):
    
    BestMotifs= SelectFirstKmer (k, Dna)
    BestMotifsScore= Score(BestMotifs)
    
    for i in Dna:
        Motifs= []
        count= CountMatrix(Dna)
        profile= ProfileMatrix(count)
        text= [i]
        rowAllKmers= AllPossibleKmers(text, k)
        mostprobablerowkmer= MostProbableK_mer(rowAllKmers, profile)
        Motifs.append(mostprobablerowkmer)
        if Score(Motifs) < BestMotifsScore:
            BestMotifs = Motifs
        print (BestMotifs)

print (GreedyMotifSearch (['GGCGTTCAGGCA',
'AAGAATCAGTCA',
'CAAGGAGTTCGC',
'CACGTCAATCAC',
'CAATAATATTCG'], 3, 5))
