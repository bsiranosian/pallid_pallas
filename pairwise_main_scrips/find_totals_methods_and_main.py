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
    return len(same_counts)

def find_GATC_counts(pairwise, reference):
    reference_counter = 0
    for p in pairwise:
        if reference == "T":    
            t_ilist = Pairwise.get_t_ilist(p)
            reference_counter += len(t_ilist)
        if reference == "O":
            o_ilist = Pairwise.get_o_ilist(p)
            reference_counter += len(o_ilist)
    return reference_counter
      
        
#Main Method (flow of control starts here from command line call):
def main():
    if len(sys.argv) != 1:
        sys.exit('USAGE: python pariwise.py input_fasta \n \
        input_fasta is the location of a fasta file with the BLASTED pairwise alignment.')
    else: 
        mypath = '/Users/juliagross/Desktop/phaedrus_parsed'
        file_names = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) ]
        file_names = file_names[1:]
        
        total_GATC_t = 0
        total_GATC_o = 0
        total_same_count = 0
        for f in file_names:
            input_fasta_name = f
            #file_identifier = f[42:55]
            input_fasta = open(input_fasta_name, "r+")
            input_fasta_copy = input_fasta.readlines()
            pairwise = Pairwise.read_fasta(input_fasta_copy)  
                
            pairwise_output = open("pairwise_output", "ab")
            total_GATC_t += find_GATC_counts(pairwise, "T")
            total_GATC_o += find_GATC_counts(pairwise, "O")
            total_same_count += find_same_counts(pairwise)
            
            input_fasta.close()
            
        pairwise_output.write(str(total_GATC_t))
        pairwise_output.write('\n')
        pairwise_output.write(str(total_GATC_o))
        pairwise_output.write('\n')
        pairwise_output.write(str(total_same_count))
        
        pairwise_output.close()
         
        
if __name__ == '__main__':
    main()

