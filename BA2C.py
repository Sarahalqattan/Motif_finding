#Problem BA2C
#Given: A string Text, an integer k, and a 4 × k matrix Profile.
#Return: A Profile-most probable k-mer in Text. (If multiple answers exist, you may return any one.)

# 1. Find the porobability of each k-mer in the string 

# 2. Find the most probable k-mer

def MostProbableK_mer(text, k, profile):
    max_probability = 0.0
    most_probable_kmer = ""

    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probability = 1.0
        for j in range(k):
            nucleotide = kmer[j]
            for n in nucleotide:
                if n == 'A':
                    probability = probability * profile[0][j]
                if n == 'C':
                    probability = probability * profile[1][j]
                if n == 'G':
                    probability = probability * profile[2][j]
                if n == 'T':
                    probability = probability * profile[3][j]

        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer
            
   

print (profile_most_probable_kmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, [[0.2, 0.2, 0.3, 0.2, 0.3],
[0.4, 0.3, 0.1, 0.5, 0.1],
[0.3, 0.3, 0.5, 0.2, 0.4],
[0.1, 0.2, 0.1, 0.1, 0.2]]))
