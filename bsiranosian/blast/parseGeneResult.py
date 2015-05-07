# script to parse gene_result file downloaded from genbank and a fasta file to 
# extract the sequences of genes, given the coordinates in the file 
import sys
import os

def main():
	if len(sys.argv) != 5:
		sys.exit("USAGE: python parseGeneResult.py gene_result.txt fasta_file phage_name output_dir")

	# result file
	rf = sys.argv[1]
	# parse it and get relevant info, coordinates of the genes
	with open(rf, 'r') as rff:
		lines = [a.strip().split('\t') for a in rff.readlines()]
		coords = [[int(a[12]), int(a[13])] for a in lines[1:]]
		# sort coordinates
		def getKey(item):
			return item[0]
		coords = sorted(coords, key=getKey)
		
	# fasta file
	ff = sys.argv[2]
	with open(ff, 'r') as fff:
		lines = [a.strip().upper() for a in fff.readlines()]

		if lines[0][0] == '>':
			lines = lines[1:]
		# condense into a string
		seq = ''
		for a in lines: seq+=a 

	# loop through coordinates and make files with each of the genes 
	for c in range(len(coords)):
		out_name = os.path.join(sys.argv[4], sys.argv[3]+'_gp'+str(c+1)+'.fasta')
		seq_to_write = seq[coords[c][0]-1:coords[c][1]]
		with open(out_name, 'w') as outf:
			outf.write('>'+sys.argv[3]+'_gp'+str(c+1)+'\n')
			outf.write(seq_to_write + '\n')


if __name__ == '__main__':
	main()