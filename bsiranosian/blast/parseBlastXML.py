import sys
import xml.etree.ElementTree as ET

def main(inFile, outFile, eFilter=100):
	# load in xml
	tree = ET.parse(inFile)
	root = tree.getroot()

	with open(outFile,'w') as outf:
		# loop over HSPs in the xml
		for hit in range(len(root[8][0][4])):
			matchName = root[8][0][4][hit][2].text
			for hsp in range(len(root[8][0][4][hit][5])):
				#implementation of evalue filter here
				if float(root[8][0][4][hit][5][hsp][3].text) <= eFilter:
					outf.write('>query ' + str(hit+1) + '.' + str(hsp+1) + '\n')
					outf.write(root[8][0][4][hit][5][hsp][14].text + '\n')
					outf.write('>' + matchName + '\n')
					outf.write(root[8][0][4][hit][5][hsp][15].text + '\n')


if __name__ == '__main__':

	if len(sys.argv) == 3:
		main(sys.argv[1],sys.argv[2])
	elif len(sys.argv) == 4:
		main(sys.argv[1],sys.argv[2],float(sys.argv[3]))
	else:
		sys.exit('USAGE: python parseBlastXML.py input_XML output_fasta [e-value filter] \n \
			alignments with e-values greater than the specified value are removed if \n \
			that option is specified')
