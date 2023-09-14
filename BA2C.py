#Problem BA2C
#Given: A string Text, an integer k, and a 4 × k matrix Profile.
#Return: A Profile-most probable k-mer in Text. (If multiple answers exist, you may return any one.)

# 1. Find the porobability of each k-mer in the string 

def Probability (string, matrix):
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

# 2. Find the most probable k-mer

def MostProbableK_mer(string, k, matrix):
    sequence= {} #Make a set of k-mers and their probabilities
    max_value= 0 #Assign the maximum value to zero 
    for i in range (len(string)-k+1):
        sequence [string [i:i+k]]= Probability (string[i:i+k], matrix) #Call Probability function
    for key, value in sequence.items():
        if value > max_value:
            max_value = value #Find the maximum probability 
    sequences =[]
    for key, value in sequence.items():
        if value == max_value:
            sequences.append(key) #Return the k-mer with maximum probability
    return sequences
   

print (profile_most_probable_kmer('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, [[0.2, 0.2, 0.3, 0.2, 0.3],
[0.4, 0.3, 0.1, 0.5, 0.1],
[0.3, 0.3, 0.5, 0.2, 0.4],
[0.1, 0.2, 0.1, 0.1, 0.2]]))
