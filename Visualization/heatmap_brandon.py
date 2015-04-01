# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:10:48 2015

@author: Brandon
"""

import sys
import json
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------

# sort the data by cluster

#------------------------------------------------------------------------------

def createHeatmap(pals, zscores):
    data = np.array(json.load(open(zscores, "r")))
    x = np.array(json.load(open(pals, "r")))
    
    plt.pcolor(data,cmap=plt.cm.Reds)
    plt.xticks(np.arange(0, len(x))+0.5, x)
#   plt.yticks(np.arange(0,10)+0.5,rows)
    plt.show()
    
#------------------------------------------------------------------------------

def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python heatmap_brandon.py filename \n \
                  filename is the file path of the zscores')
    else:
        file = sys.argv[0]

    plt.hist(createHeatmap(file))
    
#------------------------------------------------------------------------------
    
plt.hist(createHeatmap("pals.txt", "zscores.txt"))
    
if __name__ == '__main__':
    main()