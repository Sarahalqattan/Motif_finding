#Problem BA2G:
#Given: Integers k, t, and N, followed by a collection of strings Dna.
#Return: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!

# Need some random number functionality and the operator.add ability.
import random
import operator

def accumulate(iterable, func=operator.add):

    it = iter(iterable)
    try:
        total = next(it)
    except StopIteration:
        return
    yield total
    for element in it:
        total = func(total, element)
        yield total

def hammingDistance(str1, str2):

        diffs = 0
        for ch1, ch2 in zip(str1, str2):
                if ch1 != ch2:
                        diffs += 1
        return diffs

def d(kmer, dna):

    k = len(kmer)
    motif = []
    totDist = 0
    for phrase in dna:
        localDist = len(phrase)+len(kmer)
        word = ""
        for i in range(len(phrase)-k+1):
            subPattern = phrase[i:i+k]
            if localDist > hammingDistance(kmer, subPattern):
                localDist = hammingDistance(kmer, subPattern)
                word = subPattern
        motif.append(word)
        totDist += localDist
    return totDist

def score(motifs):

    k = len(motifs[0])
    pattern = []
    for i in range(k):
        A=0
        C=0
        G=0
        T=0
        for string in motifs:
            if string[i]=='A':
                A+=1
            elif string[i]=='C':
                C+=1
            elif string[i]=='G':
                G+=1
            elif string[i]=='T':
                T+=1


        if A >= C and A >= G and A >= T:
            pattern.append('A')
        elif C >= G and C >= T:
            pattern.append('C')
        elif G >= T:
            pattern.append('G')
        else:
            pattern.append('T')

    pattern = "".join(pattern)

    score = 0
    for string in motifs:
        score += hammingDistance(string, pattern)
    return score

def BuildMotifs(profile, dna, k):
    
    motif = []
    for string in dna:
        bestSubStr = ''
        for i in range(len(string)+1-k):
            substr = string[i:i+k]
            prob = 1
            bestProb = -1
            for j in range(k):
                if substr[j] == 'A':
                    prob *= profile[j][0]
                elif substr[j] == 'C':
                    prob *= profile[j][1]
                elif substr[j] == 'G':
                    prob *= profile[j][2]
                elif substr[j] == 'T':
                    prob *= profile[j][3]
            if prob > bestProb:
                bestProb = prob
                bestSubStr = substr
        motif.append(bestSubStr)
    return motif

def BuildProfile(motif):

    k = len(motif[0])
    profile = [[0 for y in range(4)] for x in range(k)]
    for count in range(k):
        A = 1
        C = 1
        G = 1
        T = 1
        for string in motif:
            if string[count]=='A':
                A+=1
            elif string[count]=='C':
                C+=1
            elif string[count]=='G':
                G+=1
            elif string[count]=='T':
                T+=1
        profile[count][0] = float(A)/(A+C+G+T)
        profile[count][1] = float(C)/(A+C+G+T)
        profile[count][2] = float(G)/(A+C+G+T)
        profile[count][3] = float(T)/(A+C+G+T)
    return profile

def singleReplacementMotif(motifs, dna_i):

    k = len(motifs[0])
    profile = BuildProfile(motifs)
    
    
    kmerDensities = [0 for x in range(len(dna_i)-k+1)]
    for i in range(len(dna_i)-k+1):
        prob = 1
        for j in range(k):
            if dna_i[i+j] == 'A':
                prob *= profile[j][0]
            elif dna_i[i+j] == 'C':
                prob *= profile[j][1]
            elif dna_i[i+j] == 'G':
                prob *= profile[j][2]
            elif dna_i[i+j] == 'T':
                prob *= profile[j][3]
        kmerDensities[i] = prob

    
    normalizationTot = sum(kmerDensities)
    for i in range(len(dna_i)-k+1):
        kmerDensities[i] = kmerDensities[i]/normalizationTot

    
    kmerDensities = list(accumulate(kmerDensities))

    
    randVal = random.random()
    for i in range(len(dna_i)-k+1):
        if randVal < kmerDensities[i]:
            break
    replacementKmer = dna_i[i:i+k]

    return replacementKmer

def gibbsSampler(dna, k, N):

    t = len(dna)
    # randomly select k-mers Motifs = (Motif_1, , Motif_t) in each string from Dna
    motifs = []
    for strand in dna:
        i = random.randrange(len(strand)-k+1)
        substr = strand[i:i+k]
        motifs.append(substr)

    bestMotifs = list(motifs)
    bestMotifsScore = score(bestMotifs)
    for j in range(1,N):
        i = random.randrange(t)
        subsetMotifs = motifs[0:i]+motifs[i+1:t]
        replacementMotif = singleReplacementMotif(subsetMotifs, dna[i])
        motifs[i] = replacementMotif

        if score(motifs) < bestMotifsScore:
            bestMotifs = list(motifs)
            bestMotifsScore = score(bestMotifs)
    return bestMotifs

def multipleSeedsGibbsSampling(dna, k, t, N):

    results = gibbsSampler(dna, k, N)
    bestScore = score(results)
    bestMotifs = list(results)
    for i in range(1, t):
        results = gibbsSampler(dna, k, N)
        if score(results) < bestScore:
            bestScore = score(results)
            bestMotifs = list(results)
    return bestMotifs


with open('rosalind_ba2g-2.txt') as f:
    
    k = 15
    t = 20
    N=2000
    strings = [st.rstrip() for st in f.readlines()]
print('\n'.join(multipleSeedsGibbsSampling(strings,k,t,N)))
