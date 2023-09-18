#Problem BA2H:
#Given: A DNA string Pattern and a collection of DNA strings Dna.
#Return: DistanceBetweenPatternAndStrings(Pattern, Dna).

#Here, the CalculateHammingDistance function compares the nucleotides of two different strings to scan for any mismatches. 
#If any mismatch is found, it adds a count of 1 into a storing variable "hd". This reiterates until all nucleotides
#are scanned and hd is returned.



import math
def CalculateHammingDistance(a, b):
    hd = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            hd += 1
    return hd


# Here, the function uses the previous function to find the hamming distance between a given input 'pattern' and all possible k-mers of similar
#length in each DNA string from input 'strings'. It essentially would tally the minimum number of mismatches between 'pattern' and possible k-mers in each index of 'strings' then returns that tally.

def DistanceBetweenPatternAndStrings(pattern, strings):
    distance = 0
    k = len(pattern)
    for string in strings:
        HammingDistance = math.inf
        for i in range(len(string)-k+1):
            k_pattern = string[i:i+k]
            hd = CalculateHammingDistance(pattern, k_pattern)
            if HammingDistance > hd:
                HammingDistance = hd
        distance += HammingDistance
    return distance

print (DistanceBetweenPatternAndStrings ('AAA', ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']

