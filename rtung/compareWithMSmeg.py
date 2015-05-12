# -*- coding: utf-8 -*-
"""
Created on Tue May 05 21:21:02 2015

@author: Ray
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

n = 4

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

with open("C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\rtung\\Data"+str(n)+".txt", "r") as f:
   data = f.read()
   yay = pd.read_json(data)

yay = yay[["Name", "SubCluster", "Mer", "MerCount in Sequence", "A", "T", "C", "G"]]
yay['Seq Length'] = yay['A'] + yay['T'] + yay['C'] + yay['G']
yay['GC'] = (yay['C'] + yay['G'])/yay['Seq Length']
yay['isPal'] = yay["Mer"].map(isPal)
yay["TUD"] = getTud(yay.copy(), n)
yay["Adjusted TUD"] = getAdjustedTud(yay.copy(), n)
yay['Cluster'] = yay["SubCluster"].map(mainCluster)
#=====================================
with open("C:\\Users\\Ray\\Documents\\GitHub\\pallid_pallas\\rtung\\Myco"+str(n)+".txt", "r") as f:
   data = f.read()
   bleh = pd.read_json(data)

bleh["Name"] = "M Smeg"
bleh["Seq Length"] = bleh["A"]+bleh["C"]+bleh["G"]+bleh["T"]
bleh["isPal"] = bleh["Mer"].map(isPal)
bleh["Smeg TUD"] = getTud(bleh.copy(), n)
bleh["Smeg Adjusted TUD"] = getAdjustedTud(bleh.copy(), n)

#lol = yay[yay["Cluster"] == "M"][["Mer", "Adjusted TUD"]].copy()
lol = yay[["Mer", "Adjusted TUD"]].copy()
res1 = lol.groupby("Mer", as_index=False).mean()
a = lol.groupby("Mer").std()
res2 = pd.DataFrame(zip(list(a.index), list(a['Adjusted TUD'])), columns = ["Mer", "std"])
res = pd.merge(res1, res2, left_on = "Mer", right_on = "Mer", how = "left")
res = pd.merge(res, bleh[["Mer", "Smeg Adjusted TUD"]], left_on = "Mer", right_on = "Mer", how = "left")
res['isPal'] = res["Mer"].map(isPal)
res["aboveLower"] = (res["Smeg Adjusted TUD"] > res["Adjusted TUD"] - res["std"])
res["belowUpper"] = (res["Smeg Adjusted TUD"] < res["Adjusted TUD"] + res["std"])
res["inBounds"] = (res.aboveLower & res.belowUpper).astype(int)

#X = res[res["isPal"]==True]["Adjusted TUD"]
#Y = res[res["isPal"]==True]["Smeg Adjusted TUD"]
X = res["Adjusted TUD"]
Y = res["Smeg Adjusted TUD"]
plt.scatter(X, Y)
X = sm.add_constant(X)
model = sm.OLS(Y, X)
results = model.fit()
m = results.params[1]
b = results.params[0]
print(results.summary())
