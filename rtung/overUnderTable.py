# -*- coding: utf-8 -*-
"""
Created on Wed May 06 18:35:21 2015

@author: Ray
"""

# -*- coding: utf-8 -*-
"""
Created on Tue May 05 21:21:02 2015

@author: Ray
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

n = 6

def mainCluster(str1):
    if str1 == 'Unclustered' or len(str1) == 1:
        return str1
    return str1[0]

def isPal(seq):
    nucs = ["A", "C", "G", "T"]
    n = len(seq)
    if n%2 == 1:
        return False
    for i in range(n//2):
        if nucs.index(seq[i]) + nucs.index(seq[n - i - 1]) != 3:
            return False
    return True

def getTud(df, len1):
    nucs = ["A", "T", "C", "G"]
    for n in nucs:
        df[n+"Count"] = df["Mer"].map(lambda x: x.count(n))
    df['Expected'] = (df["A"].astype(float)/df["Seq Length"].astype(float))**df["ACount"] * (df["C"].astype(float)/df["Seq Length"].astype(float))**df["CCount"] * (df["G"].astype(float)/df["Seq Length"].astype(float))**df["GCount"] * (df["T"].astype(float)/df["Seq Length"].astype(float))**df["TCount"]
    return df['MerCount in Sequence']/(df['Expected']*(df["Seq Length"] - len1 + 1)).astype(float)

def zeroToOne(x):
    if x == 0:
        return 1
    return x

def getAdjustedTud(df, len1):
    nucs = ["A", "T", "C", "G"]
    for n in nucs:
        df[n+"Count"] = df["Mer"].map(lambda x: x.count(n))
    df['Expected'] = (df["A"].astype(float)/df["Seq Length"].astype(float))**df["ACount"] * (df["C"].astype(float)/df["Seq Length"].astype(float))**df["CCount"] * (df["G"].astype(float)/df["Seq Length"].astype(float))**df["GCount"] * (df["T"].astype(float)/df["Seq Length"].astype(float))**df["TCount"]
    return np.log(df['MerCount in Sequence'].map(zeroToOne)/(df['Expected']*(df["Seq Length"] - len1 + 1)).astype(float))
    
def mergeCountStuff(df, col):
    a = df.groupby(col).size()
    df1 = pd.DataFrame(zip(list(a.index), list(a)), columns = [col, col + ' Count'])
    df = pd.merge(df, df1, left_on = col, right_on = col, how = 'left')
    return df    
    
def clusterToInt(str1):
    if str1 == 'Unclustered':
        return 3
    return (ord(str1[0]) - ord("A"))**(1/3.0)    

with open("C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\Myco"+str(n)+".txt", "r") as f:
   data = f.read()
   bleh = pd.read_json(data)

bleh["Name"] = "M Smeg"
bleh["Seq Length"] = bleh["A"]+bleh["C"]+bleh["G"]+bleh["T"]
bleh["isPal"] = bleh["Mer"].map(isPal)
bleh["Smeg TUD"] = getTud(bleh.copy(), n)

with open("C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\Data"+str(n)+".txt", "r") as f:
   data = f.read()
   yay = pd.read_json(data)

yay = yay[["Name", "SubCluster", "Mer", "MerCount in Sequence", "A", "T", "C", "G"]]
yay['Seq Length'] = yay['A'] + yay['T'] + yay['C'] + yay['G']
yay['GC'] = (yay['C'] + yay['G'])/yay['Seq Length']
yay['isPal'] = yay["Mer"].map(isPal)
yay["TUD"] = getTud(yay.copy(), n)
yay["Adjusted TUD"] = getAdjustedTud(yay.copy(), n)
yay['Cluster'] = yay["SubCluster"].map(mainCluster)

pals = yay[yay["isPal"]].copy()
pals['Z scored Adjusted TUD'] = (yay["Adjusted TUD"] - yay["Adjusted TUD"].mean())/yay["Adjusted TUD"].std()

#==============================================================================
# zScoreCutoff = 2
# over = pals[(pals["Z scored Adjusted TUD"] > zScoreCutoff)]
# under = pals[(pals["Z scored Adjusted TUD"] < -zScoreCutoff)]
#==============================================================================

tudCutoff = 2
over = pals[(pals["TUD"] > tudCutoff)]
under = pals[(pals["TUD"] < 1.0/tudCutoff)]

overSize = over.groupby("Mer", as_index = False).size()
underSize = under.groupby("Mer", as_index = False).size()

overSize = pd.DataFrame(zip(list(overSize.index), list(overSize)), columns = ["Mer", "Over"])
underSize = pd.DataFrame(zip(list(underSize.index), list(underSize)), columns = ["Mer", "Under"])

table = pd.merge(overSize, underSize, left_on = "Mer", right_on = "Mer", how = "outer")
table = pd.merge(table, bleh[bleh["isPal"]][["Mer", "Smeg TUD"]], left_on = "Mer", right_on = "Mer", how = "outer")
table.fillna(0, inplace = True)

table.to_csv("overUnderTable"+str(n)+".csv")





