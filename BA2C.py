#Problem BA2C
#Given: A string Text, an integer k, and a 4 × k matrix Profile.
#Return: A Profile-most probable k-mer in Text. (If multiple answers exist, you may return any one.)

# 1. Find the porobability of each k-mer in the string 

def profile_most_probable_kmer(text, k, profile):
    max_probability = 0.0
    most_probable_kmer = ""

    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probability = 1.0

        for j in range(k):
            nucleotide = kmer[j]
            probability = probability * profile[nucleotide][j]

        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer
            
   

print (profile_most_probable_kmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, [[0.2, 0.2, 0.3, 0.2, 0.3],
[0.4, 0.3, 0.1, 0.5, 0.1],
[0.3, 0.3, 0.5, 0.2, 0.4],
[0.1, 0.2, 0.1, 0.1, 0.2]]))
