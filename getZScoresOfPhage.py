import numpy as np
from os import listdir
import json

n = 4
#Sets the path to find the phage fastas. This folder contains all genomes on phagesDB
folder = "C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\phageFastas"

#Take in folder, output list of fasta names
def listOfFastasInFolder(folderLocation):
    return [folderLocation +"\\" + f for f in listdir(folderLocation)]

#This takes the four nucleotides (nucs) and converts them to integerrs so that ultimately, all palindromes can be found and defined as integers. 
def palToInt(seq, nucs):
    int1 = 0
    for i in range(0, len(seq)//2):
        int1 = int1 + 4**i *nucs.index(seq[i])
    return int1
    
#n is length of a sequence. This converts integer palindrome representations back into strings
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

#inputs a .fasta file for a single genome and outputs the genome sequence as a string
def getSequenceFromFasta(file):
    with open(file, 'r') as f:
        line = f.readline()
        while line[0]=='>':
            line=f.readline()
        sequence = line.replace("\n", "")
        for line in f:
            sequence = sequence + line.replace("\n", "")
    return sequence.replace("\t", '')

#tests if an input kmer is a palindrome and outputs the appropriate truth value 
def isPal(seq, nucs):
    n = len(seq)
    if n%2 == 1:
        return False
    for i in range(n//2):
        if nucs.index(seq[i]) + nucs.index(seq[n - i - 1]) != 3:
            return False
    return True

#takes in a genome as a string and outputs the number of each palindrome   
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

def tud(nucCount, palCount, nucs, n):
    totalNucs = np.sum(nucCount)
    tud = np.zeros(4**(n//2))
    for i in range(4**(n//2)):
        pal = intToPal(i, nucs, n)
        prob = np.product([nucCount[nucs.index(j)] for j in pal]/totalNucs)
        tud[i] = palCount[i]/(prob*(totalNucs - n + 1))
    return tud

def zScores(tuds):
    for i in range(len(tuds[0])):
        sigma = np.std([tuds[j][i] for j in range(len(tuds))])
        for j in range(len(tuds)):
            tuds[j][i] = (tuds[j][i] - 1)/sigma
    return tuds

#This runs the code properly in order, finding Z scores for palindrome use
#Finds scores for each palindrome for all genomes, using fastas downloaded from phages DB
def mainThing(n, folder):
    #KEEP THIS ORDERING
    nucs = ['A', 'C', 'G', 'T']
    fastas = listOfFastasInFolder(folder)
    #fastas = ["C:\Users\Ray\Documents\GitHub\pallid_pallas\phageFastas\Yahalom-B3.fasta"]
    allTUD = []
    for i in range(len(fastas)):
        fasta = fastas[i]
        print fasta
        dna = getSequenceFromFasta(fasta)
        data = palNucCount(dna, n, nucs)
        tuds = tud(data[0], data[1], nucs, n)
        allTUD.append(list(tuds))
    z = zScores(allTUD)
    print(z)
    with open("zscores.txt", "w") as f:
        json.dump(z, f)
    with open("phages.txt", "w") as f:
        json.dump(fastas, f)
    with open("pals.txt", "w") as f:
        json.dump([intToPal(int1, nucs, n) for int1 in range(4**(n//2))], f)
    
    
    #Now output in some nice way...

#==============================================================================
# nucs = ['A', 'C', 'G', 'T']
# fasta = "C:\Users\Ray\Documents\GitHub\pallid_pallas\phageFastas\Conspiracy-A5.fasta"
# dna = getSequenceFromFasta(fasta)
# data = palNucCount(dna, n, nucs)
# zscores = zScores(data[0], data[1], nucs, n, len(dna))
#==============================================================================

mainThing(n, folder)

