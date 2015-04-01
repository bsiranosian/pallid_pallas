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

#print pals


import sys
import os

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
    #l2 = {0:'AATT', 1:'CATG', 2:'GATC', 3:'TATA', 4:'ACGT', 5:'CCGG', 6:'GCGC', 7:'TCGA', 8:'AGCT', 9:'CGCG', 10:'GGCC',  11:'TGCA', 12:'ATAT', 13:'CTAG', 14:'GTAC', 15:'TTAA'}

    keys = []
    ps = ["AAATTT", "CAATTG", "GAATTC", "TAATTA", "ACATGT", "CCATGG", "GCATGC", "TCATGA", "AGATCT", "CGATCG", "GGATCC", "TGATCA", "ATATAT", "CTATAG", "GTATAC", "TTATAA", "AACGTT", "CACGTG", "GACGTC", "TACGTA", "ACCGGT", "CCCGGG", "GCCGGC", "TCCGGA", "AGCGCT", "CGCGCG", "GGCGCC", "TGCGCA", "ATCGAT", "CTCGAG", "GTCGAC", "TTCGAA", "AAGCTT", "CAGCTG", "GAGCTC", "TAGCTA", "ACGCGT", "CCGCGG", "GCGCGC", "TCGCGA", "AGGCCT", "CGGCCG", "GGGCCC", "TGGCCA", "ATGCAT", "CTGCAG", "GTGCAC", "TTGCAA", "AATATT", "CATATG", "GATATC", "TATATA", "ACTAGT", "CCTAGG", "GCTAGC", "TCTAGA", "AGTACT", "CGTACG", "GGTACC", "TGTACA", "ATTAAT", "CTTAAG", "GTTAAC", "TTTAAA"]
    l2 = dict()

    for i in range(64):
        l2[i] = ps[i]

    possible_pals = list(range(64))
    num_phage = int(len(l1))
    x = []

    results_dir = '/Users/juliagross/Desktop/6tud-resultsdir4115'
    #All the images will be saved here
    if not os.path.isdir(results_dir):
        os.mkdir(results_dir)
    
    for c in possible_pals:
        for r in range(64): 
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

        #create file name
        file_identifier = l2[c] + ".png"
        file_name = os.path.join(results_dir, file_identifier)
        
        #plt.show() useful test, prints graphs to GUI for testing
        plt.savefig(file_name, format='png')
        plt.clf()
        x=[]
        
make_hist_list('/Users/juliagross/Desktop/Phage_Informatics_Project/kmer_analysis/tuds6.txt', 'output.png')
