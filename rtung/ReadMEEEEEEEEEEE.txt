Data4.txt: pandas DataFrame containing columns {Name of phage, Mer of length 4, # occurences of 4-Mer in phage}
Data6.rar: pandas DataFrame containing columns {Name of phage, Mer of length 6, # occurences of 6-Mer in phage}
Myco4.txt: pandas DataFrame containing columns {Name "M smeg", Mer of length 4, # occurences of 4-Mer in M Smeg}
Myco6.txt: pandas DataFrame containing columns {Name "M smeg", Mer of length 4, # occurences of 4-Mer in M Smeg}
chartThing.csv: ordering of # occurence of 4-mer in phage by cluster, with top being least and bottom being most
compareWithMSmeg.py: produces scatterplot of AdjustedTUD in Msmeg vs Average Adjusted TUD for a 4-mer
getBaseData.py: generates Data4.txt and Data6.txt (which is a .rar in github)
guessProb.py: helps guess probability of 0 counts of a given mer for given conditions
lengthAndNoPalCounts.py: scatterplot of length of phage genome vs # pals no occurrence for a given phage
orderMerByCluster.py: generates chartThing.csv
overUnderTable.py: generates overUnderTable4.csv and overUnderTable6.csv
overUnderTable4.csv: shows number of times a 4-pal is over or under expressed and the TUD of smeg
overUnderTable6.csv: shows number of times a 6-pal is over or under expressed and the TUD of smeg

Results:

4pal occurence linked to msmeg 4 pal occurrence:
   4palPhageTUDvsMSmegAdjustedTUD.png
   4palPhageTUDvsMSmegAdjustedTUDStats.PNG
   chartThingColored.xlsx

6pal no occurrence per phage not linked to sequence length:
   6numNonOccurringVSSequenceLength.png
   6numNonOccurringVSSequenceLengthStats.PNG
   
General Characterizing:
   overUnderTable4.csv
   overUnderTable6.csv

