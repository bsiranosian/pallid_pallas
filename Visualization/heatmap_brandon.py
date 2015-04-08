# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:10:48 2015

@author: Brandon
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import pylab

#------------------------------------------------------------------------------

def getKey(item):
    return item[2]

def createHeatmap(phages, pals, zscores):
    phagesList = json.load(open(phages, "r"))
    x = np.array(json.load(open(pals, "r")))
    scoresList = json.load(open(zscores, "r"))
    
    # sort data by phage clusters
    dataList = []
    for phage, scores in zip(phagesList, scoresList):
        cluster = phage.split("-")[1]
        dataList.append((cluster, phage, scores))
    dataList = sorted(dataList, key=getKey) 
    
    # convert scores into numpy array
    tempY = []
    for el in dataList:
        tempY.append(el[2])
    y = np.array(tempY)
    
    plt.pcolor(y,cmap=plt.cm.Reds)
    plt.xticks(np.arange(0, len(x))+0.5, x)
    plt.xlim(0, len(x))
    plt.ylim(0, len(y))
    plt.title("Palindromic TUD Z-Scores")
    plt.xlabel("Palindrome")
    plt.ylabel("Bacteriophage")
    plt.colorbar()
    
    pylab.savefig('PTUD_heatmap.pdf', bbox_inches='tight')
    
#------------------------------------------------------------------------------

def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python heatmap_brandon.py phages pals zscores \n \
                  phages is the file path of the phages \n \
                  pals is the file path of the palindromes \n \
                  zscores is the file path of the zscores')
    else:
        phages = sys.argv[0]
        pals = sys.argv[1]
        zscores = sys.argv[2]

    plt.hist(createHeatmap(phages, pals, zscores))
    
if __name__ == '__main__':
    main()
    
#------------------------------------------------------------------------------

phages = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/phages4.txt"
pals = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/pals4.txt"
zscores = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/zscores4.txt"

plt.hist(createHeatmap(phages, pals, zscores))