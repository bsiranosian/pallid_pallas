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
import numpy as np
import sys
import os

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
        
    def __str__(self):
        return "Pairwise()"
    
        
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
#END OF PAIRWISE CLASS
            
#Analysis Methods:
            
def find_same_counts(pairwise):
    same_counts = []
    for p in pairwise:
        t_ilist = Pairwise.get_t_ilist(p)
        o_ilist = Pairwise.get_o_ilist(p)
        same_count = 0
        for i in range(len(t_ilist)):
            if t_ilist[i] in o_ilist:
                same_count += 1
        same_counts.append(same_count)
    return same_counts

def find_g_to_base(pairwise, base):
    g_to_x = 0
    for p in pairwise:
        t_ilist = Pairwise.get_t_ilist(p)
        o_seq = Pairwise.get_o_seq(p)
        for i in t_ilist:
            oth = o_seq[i:i+4]
            transversion = base + "ATC"
            if (oth == transversion):
                g_to_x += 1
    return g_to_x

def find_degrees_off(pairwise, degree):
    total_off = 0
    for p in pairwise:
        t_ilist = Pairwise.get_t_ilist(p)
        o_seq = Pairwise.get_o_seq(p)
        tem = "GATC"
        for i in t_ilist:
            oth = o_seq[i:i+4]
            same_count = 0
            for t in range(len(tem)):
                if oth[t] == tem[t]:
                    same_count += 1
            if same_count == (4 - degree):
                total_off += 1
    return total_off
      
        
#Main Method (flow of control starts here from command line call):
def main():
    if len(sys.argv) != 2:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
        input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        #Generates list of pairwise objects from fasta file, assigns it to pairwise        
        file_id = sys.argv[1]
        input_fasta_name = '/Users/juliagross/Desktop/' + file_id
        input_fasta = open(input_fasta_name, "r+")
        input_fasta_copy = input_fasta.readlines()
        pairwise = Pairwise.read_fasta(input_fasta_copy)
        
        #Runs Analyis functions
        
        """Same Count test code"""
        #same_counts = find_same_counts(pairwise)
        #print np.sum(same_counts)
        
        """G to X test code"""
        #g_to_a = find_g_to_base(pairwise, "A")
        #g_to_t = find_g_to_base(pairwise, "T")
        #g_to_c = find_g_to_base(pairwise, "C")
        #print "G to A: " + str(g_to_a)
        #print "G to T: " + str(g_to_t)
        #print "G to C: " + str(g_to_c)
        
        """Degree off test code:"""
        #none_off = find_degrees_off(pairwise, 0)
        #one_off = find_degrees_off(pairwise, 1)
        #two_off = find_degrees_off(pairwise, 2)
        #three_off = find_degrees_off(pairwise, 3)
        #four_off = find_degrees_off(pairwise, 4)
        #print "0 off: " + str(none_off)
        #print "1 off: " + str(one_off)
        #print "2 off: " + str(two_off)
        #print "3 off: " + str(three_off)
        #print "4 off: " + str(four_off)
        
        
        
        #Writes analysis results to output file
        
        
        #Pairwise.print_pairwise_list(pairwise)
        
        #Closes file
        input_fasta.close()
        
if __name__ == '__main__':
    main()
    
    

