#Problem BA2F 
#Given: Positive integers k and t, followed by a collection of strings Dna.
#Return: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times. Remember to use pseudocounts!


import random

def Score(Motifs):
    Score = 0
    for i in range(len(Motifs[0])):
        j = [motif[i] for motif in Motifs]
        Score += (len(j) - max(j.count("A"), j.count("C"), j.count("T"), j.count("G")))
    return Score

def Motifs_Profile(Motifs):
    d = {}
    n = float(len(Motifs))
    z = list(zip(*Motifs))
    for i in range(len(z)):
        d.setdefault('A', []).append((z[i].count('A')+1)/n/2)
        d.setdefault('C', []).append((z[i].count('C')+1)/n/2)
        d.setdefault('G', []).append((z[i].count('G')+1)/n/2)
        d.setdefault('T', []).append((z[i].count('T')+1)/n/2)
    return d

def Most_Prob_kmer(text, k , matrix):
    maxp = None
    Probable_kmer = None
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        pt = 1
        for j in range(k):
            p = matrix[kmer[j]][j]
            pt *=p
        if maxp == None or pt > maxp:
            maxp = pt
            Probable_kmer = kmer
    return Probable_kmer


def Random_kmers(DNA,k):
    Motifs = []
    for i in DNA:
        Position = random.randrange(0,len(DNA[0])-k+1)
        Motifs.append(i[Position:Position+k])
    return Motifs

def Randomized_Motif_Search(DNA,k,t):
    Motifs = Random_kmers(DNA,k)
    best = Motifs
    while True:
        matrix = Motifs_Profile(Motifs)
        Motifs = []
        for m in range(t):
            Motifs.append(Most_Prob_kmer(DNA[m],k,matrix))
        if Score(Motifs) < Score(best):
            best = Motifs
        else:
            return best


def Random_Times(dna,k,t):
    Best_Motifs = []
    High_Score = None
    for i in range(int(1000)):
        Temp_Motifs = Randomized_Motif_Search(dna, k, t)
        Temp_Score = Score(Temp_Motifs)
        if High_Score == None or High_Score > Temp_Score:
            High_Score = Temp_Score
            Best_Motifs = Temp_Motifs
    return Best_Motifs



print ( Random_Times(['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'] ,8, 5))
