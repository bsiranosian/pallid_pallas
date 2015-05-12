# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 19:59:26 2015

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
yay['Cluster'] = yay["SubCluster"].map(mainCluster)
yay['Div'] = (yay['A']*yay['A'] + yay['C']*yay['C'] + yay['G']*yay['G'] + yay['T']*yay['T'])/(yay['Seq Length']*yay['Seq Length'])

blah = yay[yay["isPal"]][["Name", "MerCount in Sequence"]].groupby("Name").sum()
blah1 = pd.DataFrame(zip(list(blah.index), list(blah["MerCount in Sequence"])), columns = ["Name", "PalCount in Sequence"])
lol = pd.merge(blah1, yay[["Name", "Seq Length", "SubCluster", "GC", "A", "T", "C", "G"]].drop_duplicates(), left_on = "Name", right_on = "Name", how = 'left')
lol["Pal Density"] = lol["PalCount in Sequence"]/(lol["Seq Length"] - n + 1)
#print(lol[["Name", "SubCluster", "Pal Density"]].sort("Pal Density"))


zeroPals = yay[(yay['TUD'] == 0) & yay['isPal']].copy()
zeroPals.drop(['TUD', 'MerCount in Sequence', 'isPal'], axis=1, inplace=True)
for i in ["Name", "Mer"]:
    zeroPals = mergeCountStuff(zeroPals, i)

zeroNotPals = yay[(yay['TUD'] == 0) & (yay['isPal'] == False)].copy()
zeroNotPals.drop(['TUD', 'MerCount in Sequence', 'isPal'], axis=1, inplace=True)
for i in ["Name", "Mer"]:
    zeroNotPals = mergeCountStuff(zeroNotPals, i)

palPlot = zeroPals[["Name", "Name Count"]].drop_duplicates()
notPalPlot = zeroNotPals[["Name", "Cluster", "SubCluster", "GC", "Name Count", "Seq Length", "Div"]].drop_duplicates()
palPlot = palPlot.rename(columns={'Name Count': 'Pal Name Count'})
notPalPlot = notPalPlot.rename(columns={'Name Count': 'Not Pal Name Count'})
toPlot = pd.merge(palPlot, notPalPlot, left_on = "Name", right_on = "Name", how = 'outer')
toPlot.fillna(0, inplace = True)



#==============================================================================
# xmin = .5
# xmax = .72
# Y = lol["Pal Density"]#/toPlot["Pal Name Count"]
# X = lol["Seq Length"]
# plt.scatter(X, Y, c = toPlot["Cluster"].apply(clusterToInt))
#==============================================================================

#==============================================================================
# xmin = .5
# xmax = .72
# Y = toPlot["Not Pal Name Count"]#/toPlot["Pal Name Count"]
# X = toPlot["GC"]
# plt.scatter(X, Y, c = toPlot["Cluster"].apply(clusterToInt))
#==============================================================================

xmin = 40000
xmax = 180000
Y = toPlot["Pal Name Count"]#/toPlot["Pal Name Count"]
X = toPlot["Seq Length"]
plt.scatter(X, Y, c = toPlot["Cluster"].apply(clusterToInt))

#==============================================================================
# xmin = 0
# xmax = 1
# Y = toPlot["Not Pal Name Count"]/(4**6 - 4**3)
# X = toPlot["Pal Name Count"]/4**3
# plt.scatter(X, Y, c = toPlot["Cluster"].apply(clusterToInt))
# plt.axis([0, 1, 0, .25])
#==============================================================================

X = sm.add_constant(X)
model = sm.OLS(Y, X)
results = model.fit()
m = results.params[1]
b = results.params[0]
print(results.summary())
plt.plot([xmin, xmax], [m*xmin + b, m*xmax + b])

#==============================================================================
# zeroPals['Density'] = zeroPals['Name Count']/zeroPals["Seq Length"]
# 
# xmin = .5
# xmax = .72
# toPlot = zeroPals[["Name", "Cluster", "SubCluster", "GC", "Name Count", "Seq Length", "Density"]].drop_duplicates()
# #toPlot = toPlot[toPlot.Cluster != "A"]
# plt.scatter(toPlot["GC"], toPlot["Density"], c = toPlot["Cluster"].apply(clusterToInt))
# plt.axis([xmin, xmax, 0, 0.0012])
# plt.title("GC Content vs Density of 6-pals \n with TUD 0 For Phage with Pals of TUD 0", fontsize=20)
# plt.ylabel('Number of TUD 0 over seq len', fontsize=14)
# plt.xlabel('GC Content', fontsize=14)
# 
# Y = toPlot["Density"]
# X = toPlot["GC"]
# X = sm.add_constant(X)
# model = sm.OLS(Y, X)
# results = model.fit()
# m = results.params[1]
# b = results.params[0]
# print(results.summary())
# plt.plot([xmin, xmax], [m*xmin + b, m*xmax + b])
#==============================================================================
