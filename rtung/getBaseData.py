# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 17:09:00 2015

@author: Ray
"""

import numpy as np
from os import listdir
import json
import pandas as pd

import sys

#Take in folder, output list of fasta names
def listOfFastasInFolder(folderLocation):
    fastas = [folderLocation +"\\" + f for f in listdir(folderLocation)]
    fastas.remove(folderLocation +"\\" + "badNames.txt")
    return fastas

#This takes the four nucleotides (nucs) and converts them to integerrs so that ultimately, all palindromes can be found and defined as integers. 
def palToInt(seq, nucs):
    int1 = 0
    for i in range(0, len(seq)):
        int1 = int1 + 4**i *nucs.index(seq[i])
    return int1
    
#n is length of a sequence. This converts integer palindrome representations back into strings
def intToPal(int1, nucs, n):
    seq = ""
    for i in range(n):
        num = int1//4**(n - i - 1)
        int1 = int1 - num*4**(n - i - 1)
        seq = seq + nucs[num]
    return seq[::-1]

#inputs a .fasta file for a single genome and outputs the genome sequence as a string
def getSequenceFromFasta(file1):
    with open(file1, 'r') as f:
        line = f.readline()
        while line[0]=='>':
            line=f.readline()
        sequence = line.replace("\n", "")
        for line in f:
            sequence = sequence + line.replace("\n", "")
    if 'N' in sequence:
        return ""
    return sequence.replace("\t", '')

#tests if an input kmer is a palindrome and outputs the appropriate truth value 
def isPal(seq, nucs):
    return True

#takes in a genome as a string and outputs the number of each palindrome   
def palNucCount(seq, n, nucs):
    counts = np.zeros(4)
    palCount = np.zeros(4**(n))
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
    tud = np.zeros(4**(n))
    for i in range(4**(n)):
        pal = intToPal(i, nucs, n)
        prob = np.product([nucCount[nucs.index(j)] for j in pal]/totalNucs)
        tud[i] = palCount[i]/(prob*(totalNucs - n + 1))
    return tud

def clean(str1):
    if str1 == "None" or str1 =="Singleton":
        return "Unclustered"
    return str1


folder = "C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\phageFastas"
n = 4
myco = True


#KEEP THIS ORDERING
badFastas = []
nucs = ['A', 'C', 'G', 'T']
fastas = listOfFastasInFolder(folder)
if myco:
    fastas = ["C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\MycobacteriumSmegmatis.fasta"]
allnMer = []
allCounts = []
for i in range(len(fastas)):
    fasta = fastas[i]
    print fasta
    dna = getSequenceFromFasta(fasta)
    if dna == "":
        badFastas.append(fasta)
    else:
        data = palNucCount(dna, n, nucs)
        allnMer.append(list(data[1]))
        allCounts.append(list(data[0]))
print(badFastas)
for bad in badFastas:
    fastas.remove(bad)
nicePhageNames = [f1.replace(folder + "\\", "").replace(".fasta", "") for f1 in fastas]

nmerlist = [intToPal(int1, nucs, n) for int1 in range(4**n)]

nblah = []
for i in range(len(nicePhageNames)):
    for j in range(len(nmerlist)):
        if myco:
            nblah.append([nicePhageNames[i], nmerlist[j], allnMer[i][j]])
        else:
            nblah.append([nicePhageNames[i].split("-")[0], clean(nicePhageNames[i].split("-")[1]), nmerlist[j], allnMer[i][j]])

if myco:
    nblah = pd.DataFrame(nblah, columns = ["Name", "Mer", "MerCount in Sequence"])
else:
    nblah = pd.DataFrame(nblah, columns = ["Name", "SubCluster", "Mer", "MerCount in Sequence"])

nucCount = [[i[0].split("-")[0]] + i[1] for i in zip(nicePhageNames, allCounts)]
nucCount = pd.DataFrame(nucCount, columns = ["Name", 'A', 'C', 'G', 'T'])

res = pd.merge(nblah, nucCount, left_on = "Name", right_on = "Name", how = 'left')
if myco:
    res.to_json("Myco"+str(n)+".txt")
else:
    res.to_json("Data" + str(n)+ ".txt")
#==============================================================================
#     with open("phages"+ str(n)+"All.txt", "w") as f:
#         json.dump(nicePhageNames, f)
#     with open("pals"+ str(n)+"All.txt", "w") as f:
#         json.dump([intToPal(int1, nucs, n) for int1 in range(4**n)], f)
#     with open("tuds"+ str(n)+"All.txt", "w") as f:
#         json.dump(allTUD, f)
#     with open("gc"+ str(n)+"All.txt", "w") as f:
#         json.dump(gc, f)
#==============================================================================