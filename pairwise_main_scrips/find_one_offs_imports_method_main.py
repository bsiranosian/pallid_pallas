#Imports
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from os import listdir
from os.path import isfile, join

def find_one_offs(pairwise, base):
    counter = 0
    for p in pairwise:
        t_ilist = Pairwise.get_t_ilist(p)
        o_seq = Pairwise.get_o_seq(p)
        tem = "GATC"
        for i in t_ilist:
            oth = o_seq[i:i+4]
            if (oth != tem) & (base == "G") & (oth[1:4] == "ATC"):
                counter += 1
            elif (oth != tem) & (base == "A") & (oth[0] == "G") & (oth[2:4] == "TC"):
                counter += 1
            elif (oth != tem) & (base == "T") & (oth[0:2] == "GA") & (oth[3] == "C"):
                counter += 1
            elif (oth != tem) & (base == "C") & (oth[0:3] == "GAT"):
                counter += 1
                
    return counter
      
#Main Method (flow of control starts here from command line call):
def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
        input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        mypath = '/Users/juliagross/Desktop/phaedrus_parsed'
        file_names = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) ]
        file_names = file_names[1:]

        total_A_off = 0
        total_T_off = 0
        total_G_off = 0
        total_C_off = 0
        
        for f in file_names:
            input_fasta_name = f
            input_fasta = open(input_fasta_name, "r+")
            input_fasta_copy = input_fasta.readlines()
            pairwise = Pairwise.read_fasta(input_fasta_copy)
            total_A_off += find_one_offs(pairwise, "G")
            total_T_off += find_one_offs(pairwise, "A")
            total_G_off += find_one_offs(pairwise, "T")
            total_C_off += find_one_offs(pairwise, "C")
            
            input_fasta.close()
        
        output_x = [0, 1, 2, 3]
        output_data = [total_A_off, total_T_off, total_G_off, total_C_off]
           
        #Output File Generation                
        plt.bar(output_x, output_data, align='center')
        plt.title("Off the one offs: ")
        #plt.xlabel("Base not present")
        plt.ylabel("Total Instances")
        plt.savefig('pairwise_output.png')
         
        
if __name__ == '__main__':
    main()