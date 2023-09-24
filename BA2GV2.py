import random

# Function to randomly select a k-mer index from a sequence based on probabilities
def random_kmer_index(sequence, k, profile_matrix):
    probabilities = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        probabilities.append(calculate_kmer_probability(kmer, profile_matrix))
    total_prob = sum(probabilities)
    normalized_probs = [prob / total_prob for prob in probabilities]
    return random.choices(range(len(sequence) - k + 1), weights=normalized_probs, k=1)[0]

# Function to calculate the probability of a k-mer given a profile matrix
def calculate_kmer_probability(kmer, profile_matrix):
    prob = 1
    for i in range(len(kmer)):
        prob *= profile_matrix[kmer[i]][i]
    return prob
    
# Function to calculate the profile matrix with pseudocounts
def calculate_profile_matrix(motifs, pseudocount=1):
    k = len(motifs[0]) # index 0
    profile_matrix = {'A': [pseudocount] * k, 'C': [pseudocount] * k, 'G': [pseudocount] * k, 'T': [pseudocount] * k}
    for motif in motifs:
        for i in range(k):
            profile_matrix[motif[i]][i] += 1
            print(profile_matrix[motif[i]][i])
    return profile_matrix

# Gibbs sampling function
def gibbs_sampler(Dna, k, t, N):
    best_motifs = [sequence[:k] for sequence in Dna]
    current_motifs = best_motifs.copy()
    for _ in range(N):
        i = random.randint(0, t - 1)  # Generate a random integer between 0 and t-1
        motifs_except_i = current_motifs[:i] + current_motifs[i+1:]
        profile = calculate_profile_matrix(motifs_except_i)
        selected_index = random_kmer_index(Dna[i], k, profile)
        current_motifs[i] = Dna[i][selected_index:selected_index + k]
        
        current_score = score_motifs(current_motifs)
        best_score = score_motifs(best_motifs)
        
        if current_score < best_score:
            best_motifs = current_motifs.copy()
    
    return best_motifs

# Main function to run GibbsSampler with 20 random starts
def main(Dna, k, t, N, num_starts=20):
    best_score = float('inf')
    best_motifs = None
    for _ in range(num_starts):
        random_starts = [random.randint(0, len(Dna[0]) - k) for _ in range(t)]
        motifs = [Dna[i][start:start+k] for i, start in enumerate(random_starts)]
        current_motifs = gibbs_sampler(motifs, k, t, N)
        current_score = score_motifs(current_motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs
    return best_motifs

# Function to calculate the score of a set of motifs
def score_motifs(motifs):
    t = len(motifs)
    k = len(motifs[0])
    score = 0
    for j in range(k):
        column = [motifs[i][j] for i in range(t)]
        most_common = max(set(column), key=column.count)
        score += sum(1 for nucleotide in column if nucleotide != most_common)
    return score

# Example usage
k = 8
t = 5
N = 100
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

best_motifs = main(Dna, k, t, N)
for motif in best_motifs:
    print(motif)
