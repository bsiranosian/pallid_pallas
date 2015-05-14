#Main Method (flow of control starts here from command line call):
def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
        input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        mypath = '/Users/juliagross/Desktop/phaedrus_parsed'
        file_names = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) ]
        file_names = file_names[1:]

        for f in file_names:
            input_fasta_name = f
            file_identifier = f[42:55]
            input_fasta = open(input_fasta_name, "r+")
            input_fasta_copy = input_fasta.readlines()
            pairwise = Pairwise.read_fasta(input_fasta_copy)  
                
            pairwise_output = open("pairwise_output", "ab")
            pairwise_output.write(file_identifier + ": " + str(find_g_to_base(pairwise, "A")) + "   ") 
            pairwise_output.write(file_identifier + ": " + str(find_g_to_base(pairwise, "T")) + "   ") 
            pairwise_output.write(file_identifier + ": " + str(find_g_to_base(pairwise, "C")) + "   ") 
            pairwise_output.write('\n')
            
            input_fasta.close()
        pairwise_output.close()
         
        
if __name__ == '__main__':
    main()