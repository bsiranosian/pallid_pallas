import numpy as np

n = 4
folder = "YAY"

#TODO Finish this method. Take in folder, output list of fasta names
def listOfFastasInFolder(folderLocation):
    return []

def palToInt(seq, nucs):
    int1 = 0
    for i in range(0, len(seq)//2):
        int1 = int1 + 4**i *nucs.index(seq[i])
    return int1
    
#n is length of a sequence
def intToPal(int1, nucs, n):
    seq = ""
    for i in range(n//2):
        num = int1//4**(n//2 - i - 1)
        int1 = int1 - num*4**(n//2 - i - 1)
        seq = seq + nucs[num]
    otherHalf = ""
    for i in seq:
        otherHalf = otherHalf + nucs[3-nucs.index(i)]
    return seq[::-1] + otherHalf

#TODO not yet tested
def getSequenceFromFasta(file):
    with open(file, 'r') as f:
        line = f.readline()
        while line[0]=='>':
            line=f.readline()
        sequence = line.replace("\n", "")
        for line in f:
            sequence = sequence + line.replace("\n", "")
    return sequence
    
def isPal(seq, nucs):
    n = len(seq)
    if n%2 == 1:
        return False
    for i in range(n//2):
        if nucs.index(seq[i]) + nucs.index(seq[n - i - 1]) != 3:
            return False
    return True
    
def palNucCount(seq, n, nucs):
    counts = np.zeros(4)
    palCount = np.zeros(4**(n//2))
    for i in range(n - 1):
        counts[nucs.index(seq[i])] = counts[nucs.index(seq[i])] + 1
    for i in range(n , len(seq) + 1):
        counts[nucs.index(seq[i-1])] = counts[nucs.index(seq[i-1])] + 1
        ntide = seq[i - n: i]
        if isPal(ntide, nucs):
            palCount[palToInt(ntide, nucs)] = palCount[palToInt(ntide, nucs)] + 1
    return (counts, palCount)

#TODO not yet tested
def zScores(nucCount, palCount, nucs, n, dnaLength):
    totalNucs = np.sum(nucCount)
    zScores = np.zeros(4**(n//2))
    for i in range(4**(n//2)):
        pal = intToPal(i, nucs, n)
        prob = np.product([nucCount[nucs.index(j)] for j in pal]/totalNucs)
        sigma = np.sqrt((dnaLength - n + 1)*prob*(1- prob))
        zScores[i] = (palCount[i] - prob * (dnaLength - n + 1))/sigma
    return zScores
    

    
#KEEP THIS ORDERING
nucs = ['A', 'C', 'G', 'T']
fastas = listOfFastasInFolder(folder)
for fasta in fastas:
    dna = getSequenceFromFasta(fasta)
    data = palNucCount(dna, n, nucs)
    zscores = zScores(data[0], data[1], nucs, n, len(dna))
    #Now output in some nice way...
    

