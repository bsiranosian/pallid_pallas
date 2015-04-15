# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:25:22 2015

@author: juliagross
"""
#Imports, and loads in the master 6-mer TUD array:

import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('multipage.pdf')
matplotlib.get_backend()
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import json

#Get 6mer TUD array, assign it to tuds
file_name = '/Users/juliagross/Desktop/Phage_Informatics_Project/kmer_analysis/tuds6.txt'
with open(file_name, 'r') as f:
    tuds = json.load(f)

#Create pals dict:
ps = ["AAATTT", "CAATTG", "GAATTC", "TAATTA", "ACATGT", "CCATGG", "GCATGC", "TCATGA", "AGATCT", "CGATCG", "GGATCC", "TGATCA", "ATATAT", "CTATAG", "GTATAC", "TTATAA", "AACGTT", "CACGTG", "GACGTC", "TACGTA", "ACCGGT", "CCCGGG", "GCCGGC", "TCCGGA", "AGCGCT", "CGCGCG", "GGCGCC", "TGCGCA", "ATCGAT", "CTCGAG", "GTCGAC", "TTCGAA", "AAGCTT", "CAGCTG", "GAGCTC", "TAGCTA", "ACGCGT", "CCGCGG", "GCGCGC", "TCGCGA", "AGGCCT", "CGGCCG", "GGGCCC", "TGGCCA", "ATGCAT", "CTGCAG", "GTGCAC", "TTGCAA", "AATATT", "CATATG", "GATATC", "TATATA", "ACTAGT", "CCTAGG", "GCTAGC", "TCTAGA", "AGTACT", "CGTACG", "GGTACC", "TGTACA", "ATTAAT", "CTTAAG", "GTTAAC", "TTTAAA"]
master_pals = dict()
for i in range(64):
    master_pals[i] = ps[i]
#Useful Test: print master_pals which will produce a pals dict, #(0-63):PAL

#This transposes the 6-mer TUDs array, so that each list is a list of TUDS for every given phage w/in a Pal (as opposed to the other way around)
mp3 = np.transpose(tuds)


#Lines 19-34 sort the 6-mers into bins based on GC content:
# 0 GC, 2 GC, 4 GC, or 6 GC.
gccounter = 0
gc0 = []
gc2 = []
gc4 = []
gc6 = []
for pal in ps:
    for letter in pal:
        if letter == "G" or letter == "C":
            gccounter += 1
    if gccounter == 0:
        gc0.append(pal)
    if gccounter == 2:
        gc2.append(pal)
    if gccounter == 4:
        gc4.append(pal)
    if gccounter == 6:
        gc6.append(pal)
    gccounter = 0

#Initialize counters for each GC bin
gcs0 = 0
gcs2 = 0
gcs4 = 0
gcs6 = 0

#Sets up the correct lengths
tuds_len = len(tuds)
pals_len =  len(tuds[0])

#Iterates through the transposed TUDS array, by pal lists on the outside and individual phage TUDS within a given pal on the inside
for pal in range (pals_len): #iterates through the list of pals, one pal at a time for as many pals as the list contains
    name = master_pals[pal] #sets name = to the appropriate name from the master dict
    for i in range(tuds_len): #uses an arbitrary i as a placeholder to iterate through 
        phage = mp3[pal][i]
        if (name in gc0) and (phage == 0):
            gcs0 += 1
        if (name in gc2) and (phage == 0):
            gcs2 += 1
        if (name in gc4) and (phage == 0):
            gcs4 += 1
        if (name in gc6) and (phage == 0):
            gcs6 += 1
            
total0 = gcs0 + gcs2 + gcs4 + gcs6 # Sums all the 0s across all categories
total_data = (pals_len * tuds_len) # Convenient name for total # of datapoints
x = ["GC0", "GC2", "GC4", "GC6"] #X-labels for eventual graph
raw_y = [gcs0, gcs2, gcs4, gcs6] #Raw data, the # of 0 occurences in each bin
normed_y = [(gcs0/len(gc0)), (gcs2/len(gc2)), (gcs4/len(gc4)), (gcs6/len(gc6))] # Normed data, raw # / size of each bin

#Super unsophisticated output generation becase Julia doesn't want to deal w/ matplotlib and file handling, until we know exactly what we want in the figures
print total0
print total0 / float(total_data)
print x
print raw_y
print normed_y






