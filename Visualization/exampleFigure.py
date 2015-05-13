# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:10:48 2015

@author: Brandon

Creates an example heatmap that displays our data.
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------

def getMajorCluster(item):
    return item[0]
def getSubCluster(item):
    return item[1]

def createHeatmap(phages, pals, scores):
    phagesList = json.load(open(phages, "r"))
    x = np.array(json.load(open(pals, "r")))
    scoresList = json.load(open(scores, "r"))
    
    # sorts data by phage clusters
    dataList = []
    singleton = []
    unclustered = []
    none = []
    for phage, scores in zip(phagesList, scoresList):
        subCluster = phage.split("-")[1][1:]
        majorCluster = phage.split("-")[1][0]

        if subCluster == "ingleton":
            singleton.append(("Singleton", "N/A", phage, scores))
        elif subCluster == "nclustered":
            unclustered.append(("Unclustered", "N/A", phage, scores))
        elif subCluster == "one":
            none.append(("None", "N/A", phage, scores))
        else:
            dataList.append((majorCluster, subCluster, phage, scores))
    dataList = dataList + singleton + unclustered + none
    dataList.sort(key=getMajorCluster)
    
    # gets number of clusters in dataList
    count = 0
    prev = None
    for l in dataList:
        if l[0] != prev:
            count += 1
        prev = l[0]
    
    # creates list of cluster lists
    clustersList = [[] for i in range(count)]
    for l in dataList:
        
        if l[0] == "Unclustered":
            clustersList[len(clustersList) - 3].append(l)
        elif l[0] == "Singleton":
            clustersList[len(clustersList) - 2].append(l)
        elif l[0] == "None":
            clustersList[len(clustersList) - 1].append(l)
        else:
            clustersList[ord(l[0])-65].append(l)
            
    # sorts by subcluster
    for cluster in clustersList:
        cluster.sort(key=getSubCluster)
    
    # filters out all non B clusters, where clustersList[1] represents all B clusters
    bCluster = clustersList[1]
    
    # Attaches 3 of each B-subcluster page to finalList
    finalList = []
    usedClusters = []
    ylabels = []
    for phage in bCluster:
        if usedClusters.count(phage[1]) < 3:
            usedClusters.append(phage[1])
            finalList.append(phage)
    
    for item in usedClusters:
        toAppend = "B" + item
        if toAppend not in ylabels:
            ylabels.append("")
            ylabels.append("")
            ylabels.append(toAppend)
        
#==============================================================================    
# Creates a single heatmap containing all clusters  
    
    # writes phages to csv file
    file = open("heatmap_visuals/example_heatmap_scores.csv", "w")
    file.write(",")
    for pal in x:
        file.write(pal + ",")
    file.write("\n")
    
    for phage in finalList:
        file.write(phage[2] + ",")
        for score in phage[3]:
            file.write(str(score) + ",")
        file.write("\n")
    file.close
    
    # convert scores into numpy array
    tempY = []
    for l in finalList:
        tempY.append(l[3])
    y = np.array(tempY)
    
    # creates the heatmap
    heatmap = plt.figure(1)
    plot = heatmap.add_subplot(1,1,1)
    plt.pcolor(y,cmap=plt.cm.Reds)
    
    # x labels
    plot.tick_params(axis='x', which='major', labelsize=5) 
    plt.xticks(np.arange(0, len(x))+0.5, x)    
    plt.xlim(0, len(x))
    
    # y labels
    print(len(y))
    plt.yticks(np.arange(0, len(y))+1, ylabels)
    plt.ylim(0, len(y))
    
    plt.title("Palindromic TUD Scores")
    plt.xlabel("Palindrome")
    plt.ylabel("Bacteriophages by Cluster")
    plt.colorbar()

    fig = plt.gcf()
    xsize = 8.5
    ysize = 6.0
    fig.set_size_inches(xsize, ysize)
    
    fig.savefig("heatmap_visuals/example_heatmap.pdf", bbox_inches='tight')
    plt.close(fig)
    
#==============================================================================

phages = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/phages4.txt"
pals = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/pals4.txt"
zscores = "C:/Users/Brandon/Documents/GitHub/pallid_pallas/datas/tuds4.txt"

plt.hist(createHeatmap(phages, pals, zscores), bins=1)

#==============================================================================

def main():
    if len(sys.argv) != 3:
        sys.exit('USAGE: python heatmap_brandon.py phages pals zscores \n \
                  phages is the file path of the phages \n \
                  pals is the file path of the palindromes \n \
                  zscores is the file path of the zscores')
    else:
        phages = sys.argv[0]
        pals = sys.argv[1]
        zscores = sys.argv[2]

    plt.hist(createHeatmap(phages, pals, zscores), bins=1)
    
if __name__ == '__main__':
    main()