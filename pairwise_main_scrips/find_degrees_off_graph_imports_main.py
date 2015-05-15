import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from os import listdir
from os.path import isfile, join


#... (see any infrasatructure script)

#Main Method (flow of control starts here from command line call):
def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
        input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        mypath = '/Users/juliagross/Desktop/phaedrus_parsed'
        file_names = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) ]
        file_names = file_names[1:]

        total_0_off = 0
        total_1_off = 0
        total_2_off = 0
        total_3_off = 0
        total_4_off = 0
        
        for f in file_names:
            input_fasta_name = f
            input_fasta = open(input_fasta_name, "r+")
            input_fasta_copy = input_fasta.readlines()
            pairwise = Pairwise.read_fasta(input_fasta_copy)
            total_0_off += find_degrees_off(pairwise, 0)
            total_1_off += find_degrees_off(pairwise, 1)
            total_2_off += find_degrees_off(pairwise, 2)
            total_3_off += find_degrees_off(pairwise, 3)
            total_4_off += find_degrees_off(pairwise, 4)
            
            input_fasta.close()
        
        output_x = [0, 1, 2, 3, 4]
        output_data = [total_0_off, total_1_off, total_2_off, total_3_off, total_4_off]
           
        #Output File Generation                
        plt.bar(output_x, output_data, align='center')
        plt.title("Where B3 phage have GATC, the other Bs have: ")
        plt.xlabel("Bases Different")
        plt.ylabel("Total Instances")
        plt.savefig('pairwise_output.png')
         
        
if __name__ == '__main__':
    main()