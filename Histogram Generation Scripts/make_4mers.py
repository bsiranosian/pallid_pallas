# -*- coding: utf-8 -*-
"""
This script takes in a 2-d array of zscores, pulls that information out of the
json file (formatted such that each row is a different phage, and each column
is a different pal), writes a file of histograms (one for each pal, all
labeled) to an output file of whatever name the user specifies. # of bins is
arbitrarily set to 100 for now, could easily be made to be variable.

make_pal_histograms(z, p*, o)
z, required: the array of zscores
p, aspirationally required: list of pals (will be required as soon as I figure
out how to get it to work, for now it's hard-coded the 4-mer pals in Ray's order)
o, required: name of desired output file

"""

import sys

def make_hist_list(z, o):
    import numpy as np
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('multipage.pdf')
    matplotlib.get_backend()
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    import json
    
    with open(z, 'r') as f:
        l1 = json.load(f)

    #with open('pals.txt', 'r') as f:
    #l2 = json.load(f)
    #Throws error, not sure why. Got around this by hard coding the 4-mer pals,
    #but we should probably figure this out.

    #Temporary hard-coded fix for generating pal names:
    l2 = {0:'AATT', 1:'CATG', 2:'GATC', 3:'TATA', 4:'ACGT', 5:'CCGG', 6:'GCGC', 7:'TCGA', 8:'AGCT', 9:'CGCG', 10:'GGCC',  11:'TGCA', 12:'ATAT', 13:'CTAG', 14:'GTAC', 15:'TTAA'}


    possible_pals = list(range(16))
    num_phage = int(len(l1))
    x = []


    for c in possible_pals:
        for r in range(num_phage): 
            x.append(l1[r][c])
            #print x
            mean = np.mean(x)
            std = np.std(x)
            num_bins = 100
        # the histogram of the data
        n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='blue', alpha=0.5)        
        # add a 'best fit' line
        y = mlab.normpdf(bins, mean, std)
        plt.plot(bins, y, 'r--')
        plt.xlabel('Z-scores')
        plt.ylabel('Frequency')
        plt.title('Phage Z-scores of Pal: ' + str(l2[c]))

        # Tweak spacing to prevent clipping of ylabel
        plt.subplots_adjust(left=0.15)
        np.sum(n * np.diff(bins))
        #plt.show() useful test, prints graphs to GUI for testing
        plt.savefig(pp, format='pdf')
        plt.clear()
        x=[]
        
make_hist_list('zscores.txt', 'output.png')

#My attempt to implement comandline functionality; full disclosure, I have no clue if it works
# a main function is reserved for what is run when the program is calle from the command line
#def main():
    # sys.argv is a list that contains the command line arguments. 
    # sys.argv[0] is the name of the python call. Everything after that is 
    # arguments to the script. Here we're checking for 3 additional arguments, 
    # so we want the list to be of length 4.
    #if len(sys.argv) != 3:
        # if the correct arguments aren't given, exit with a helpful message
        #sys.exit('USAGE: python make_hist_list.py input_zscores.txt output_file \n \
                  #input_zscores is the location of a 2d array populated with z scores, where rows are unique phages and columns are unique pals, \n \
                  #output_file is the name of a file to save the results to.')
    # assign variables from the arguments
   # else:
        #z= sys.argv[1]
        #o = sys.argv[2]

    # run the program you designed
    #make_hist_list(z, o)

# heres where you define all your functions
#def make_hist_list(z, o):
    #return 0
