"""This is a program designed in OOP style, which:
a. reads in a .fasta file in the following format: 
template name, template sequence, other name, other sequence
b. constructs a Pairwise object from every such 4-line group
c. has useful getter/setter methods for the 6 attributes
of pairwise objects (temp/oth name, seq, and ilist)
d. has custom analysis methods, designed for batch analysis
of an array of Pairwise objects from a given fasta
e. has a main method which that organizes flow of control, 
determines desired output for any given program run, and 
can be run from the comand-line. 

@author Julia Gross
@date 5/2/15. """

#Imports
#import numpy as np
import sys

    #Main method: This is where the flow of control is. It first calls read_fasta on a fasta and saves the output to pairwise. It next calls either one or a series of analysis methods, with the goal of producing data in a form that's ammenible to being read in and processesed by a seperate data visualization script. 

class Pairwise:
    
    #Helper function for setting up ilists:
    @staticmethod
    def find_starts(seq):
        gatc_start_indecies = []
        for i, c in enumerate(seq):
            if i+3 <= (len(seq)-1):
                tet = seq[i] + seq[i+1] + seq[i+2] + seq[i+3]
                if tet == "GATC":
                    gatc_start_indecies.append(i)
        return gatc_start_indecies
    
    
    #Pairwise object constructor, designed to read in groups of 4 lines from correctly formatted fasta files. Outputs a pairwise object with a t_name, t_seq, t_ilist, o_name, o_seq, and o_ilist.
    def __init__(self, t_name, t_seq, o_name, o_seq):
        self.t_name = t_name
        self.t_seq = t_seq
        self.o_name = o_name
        self.o_seq = o_seq
        self.t_ilist = Pairwise.find_starts(t_seq)
        self.o_ilist = Pairwise.find_starts(o_seq)
    
        
    #Getter Methods (t/o name, seq, ilist):
    def get_t_name(self):
        return self.t_name
    
    def get_o_name(self):
        return self.o_name
    
    def get_t_seq(self):
        return self.t_seq
    
    def get_o_seq(self):
        return self.o_seq
    
    def get_t_ilist(self):
        return self.t_ilist
    
    def get_o_ilist(self):
        return self.o_ilist
    
    
    #Method to actually read in the fasta. Returns array of pairwise objects.
    @staticmethod
    def read_fasta(input_fasta_copy):
        #Initializes variables
        pairwise_list = []
        t_name = ""
        t_seq = ""
        o_name = ""
        o_seq = ""
        #Checks if the file has the correct formatting, ie # of lines 
        #evenly divisible by 4
        for j, line in enumerate(input_fasta_copy):
            pass
        if j % 4 != 0:
            #print "j = " + str(j)
            return "This file is incorrectly formatted. Please try again."
        
        #If so, itterates through the file, making a pairwise object for 
        #every 4 lines, and returning a list of those objects.     
        else:    
            for i, line in enumerate(input_fasta_copy):
                if i % 4 == 0:
                    t_name = line
                if i % 4 == 1:
                    t_seq = line
                if i % 4 == 2:
                    o_name = line
                if i % 4 == 3:
                    o_seq = line
                    new_pairwise = Pairwise(t_name, t_seq, o_name, o_seq)
                    pairwise_list.append(new_pairwise)
            return pairwise_list
    
    
    #Prints a list of pariwise objects, or a single one if given only one.
    @staticmethod
    def print_single_pairwise(one):
        print Pairwise.get_t_name(one)
        print Pairwise.get_t_seq(one)
        print Pairwise.get_o_name(one)
        print Pairwise.get_o_seq(one)
        print Pairwise.get_t_ilist(one)
        print Pairwise.get_o_ilist(one)
    
    @staticmethod
    def print_pairwise_list(pairwise):
        for p in range(len(pairwise)):
            Pairwise.print_single_pairwise(pairwise[p])
            print
    
    #Analysis Functions:
    
    #This function takes in a list of pairwise objects, iterates through it based on the t_ilist and counts what happens at all those sites: 4/4 match, 3/4 match, 2/4 match, 1/4 match, 0/4 match. The function takes all the final counts in each category for each pairwise and make them into a list of 5 #s. Each of the lists is then appended to a master characterization list, which the function outputs.
    def characterize_deviations_data(pairwise):
        characterizations = []
        for p in pairwise:
            ilist = p.get_t_ilist()
            t_seq = p.get_t_seq()
            o_seq = p.get_o_seq()
            for i in ilist:
                characterization = []
                same_count = 0.0
                one_off = 0.0
                two_off = 0.0
                three_off = 0.0
                four_off = 0.0
                tem = t_seq[i:i+4]
                oth = o_seq[i:i+4]
                for j in range(len(tem)):
                    diff_score = 0
                    if tem[j] != oth[j]:
                        diff_score += 1
            if diff_score == 0:
                same_count += 1
            if diff_score == 1:
                one_off += 1
            if diff_score == 2:
                two_off += 1
            if diff_score == 3:
                three_off += 1
            if diff_score == 4:
                four_off += 1
        characterization = [same_count, one_off, two_off, three_off, four_off]
        characterizations.append(characterization)
        return characterizations
        
    
def main():
    if len(sys.argv) != 2:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
            input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        file_id = sys.argv[1]
        input_fasta_name = '/Users/juliagross/Desktop/' + file_id
        input_fasta = open(input_fasta_name, "r+")
        input_fasta_copy = input_fasta.readlines()
        pairwise = Pairwise.read_fasta(input_fasta_copy)
        Pairwise.print_pairwise_list(pairwise)
        input_fasta.close()
        
if __name__ == '__main__':
    main()
    
    

