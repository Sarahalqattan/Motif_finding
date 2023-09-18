#Problem BA2A
#Given: Integers k and d, followed by a collection of strings Dna.
#Return: All (k, d)-motifs in Dna.

def findkmers(dna, k):
    n = len(dna)
    kmers = []
    for i in range(n - k + 1):
        kmer = dna[i:i + k]
        kmers.append(kmer)
    return kmers

def mismatch(text1, text2):
    n1 = len(text1)
    count = 0
    for i in range(n1):
        if text1[i] != text2[i]:
            count += 1
    return count

def Neighbors(pattern, d):
    base = ['A', 'T', 'G', 'C']
    Neighbourhood = set()
    n2 = len(pattern)
    n3 = len(base)
    for i in range(n2):
        for j in range(n3):
            Neighbour = pattern[:i] + base[j] + pattern[i+1:]
            if d <= 1:
                Neighbourhood.add(Neighbour)
            else:
                Neighbours(Neighbour, d - 1, Neighbourhood)  
    return Neighbourhood

def MotifEnumeration(dna, k, d):
    patterns = set()
    for sequence in dna:
        kmers = findkmers(sequence, k)  
        for kmer in kmers:
            kmerpattern = Neighbors(kmer, d)
            patterns.update(kmerpattern)
    motifpatterns = []
    for pattern in patterns:
        count = 0
        for sequence in dna:
            kmers = findkmers(sequence, k)  
            for kmer in kmers:
                if mismatch(pattern, kmer) <= d:
                    count += 1
                    break 
        if count == len(dna):
            motifpatterns.append(pattern)

    return motifpatterns
