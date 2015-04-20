# make fasta map file

cd ~/projects/phage_2015
tail -n +2 ~/Downloads/PhagesDB_Data\ \(4\).txt  | tr '\t' > fasta_map_1.txt

paste <(cat fasta_map_1.txt) <(ls -d -1 $PWD/genomes/*.*) -d ',' > fasta_map.txt
# FIX STUPID ERROR ON BAEE #

cd ~/projects/tango_final/kmer_analysis
python kmerCountGenome.py ~/projects/phage_2015/fasta_map.txt 4 ~/projects/phage_2015/count_4mers.txt
python kmerCountGenome.py ~/projects/phage_2015/fasta_map.txt 6 ~/projects/phage_2015/count_6mers.txt


python compareTUD.py --k 4 --s ~/projects/phage_2015/dev_4mers.txt ~/projects/phage_2015/fasta_map.txt ~/projects/phage_2015/crap.nex  
python compareTUD.py --k 6 --s ~/projects/phage_2015/dev_6mers.txt ~/projects/phage_2015/fasta_map.txt ~/projects/phage_2015/crap.nex  

grep '\-B' fasta_map.txt > fasta_map_B.txt
grep -v 'B3' fasta_map_B.txt > fasta_map_B_noB3.txt
for i in $(cut -f2 -d, fasta_map_B_noB3.txt); do cat $i >> all_B_noB3.fasta; done
makeblastdb -dbtype 'nucl' -in all_B_noB3.fasta -title 'mycobacteriophages cluster B, no B3' 

