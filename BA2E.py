#Problem BA2E
#Given: Integers k and t, followed by a collection of strings Dna.
#Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t) with pseudocounts. If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

def KmerProbability( string, matrix):

    Probable= 1 #Assign the probability to 1
    for i in range (len(string)):
        if string [i]== 'A':
            Probable= Probable * matrix [0][i] #Multiply the value with probability of A in the 1st row and i column
        if string [i]== 'C':
            Probable= Probable * matrix [1][i] #Multiply the value with probability of C in the 2nd row and i column
        if string [i]== 'G':
            Probable= Probable * matrix [2][i] #Multiply the value with probability of G in the 3rd row and i column
        if string [i]== 'T':
            Probable= Probable * matrix [3][i] #Multiply the value with probability of T in the 4th row and i column
    return Probable


def MostProbableKmer(string, k, matrix):
    sequence = {}
    for i in range(len(string) - k + 1):
        sequence[string[i:i + k]] = KmerProbability(string[i:i + k], matrix)
    max_key = sorted(sequence.items(), key=lambda x:x[1], reverse=True)[0][0]
    return max_key


def Score(Motifs):
    score = 0
    for i in range(len(Motifs[0])):
        j = [motif[i] for motif in Motifs]
        score += (len(j) - max(j.count("A"), j.count("C"), j.count("T"), j.count("G")))
    return score


def GreedyMotifSearch(Dna, k, t):

    
    BestMotifs = [dna[:k] for dna in Dna]
    
    for k_mer in [Dna[0][i:i+k] for i in range(len(Dna[0])-k+1)]:
        
        Motifs = [k_mer]
        
        for i in range(1, t):
            
            motifs = Motifs[:i]
            
            matrix = []
            for n in ["A", "C", "G", "T"]:
                mat = []
                for j in range(k):
                    mm = [m[j] for m in motifs]
                    mat.append((mm.count(n)+1)/(t+4))
                matrix.append(mat)
               
            Motifs.append(MostProbableKmer(Dna[i], k, matrix))
        
        
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    
    return BestMotifs

print (GreedyMotifSearch (['GGCGTTCAGGCA',
'AAGAATCAGTCA',
'CAAGGAGTTCGC',
'CACGTCAATCAC',
'CAATAATATTCG'], 3, 5))
