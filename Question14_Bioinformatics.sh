### The first three awk commands are aiming to join the chromosome and position into a single id since that uniquely identifies a locus and sorts them in prep for subsequent join. 
awk -F'\t' '{print "chr"$0}' association.txt > association_with_chr_word.txt
awk -F'\t' 'BEGIN{OFS="\t";} {print $0,$1"_"$3}' association _with_chr_word.txt | sort -k7 > association_chr_position_combined.txt
awk -F'\t' 'BEGIN{OFS="\t";} {print $0,$1"_"$2}' fromVCF.txt | sort -k7 > fromVCF_chr_position_combined.txt


### This join command takes the column indexes for the combined chromosome and start position and makes a file 
join -1 7 -2 7 -t $'\t' association_chr_position_combined.txt fromVCF_chr_position_combined.txt | cut -f2,3,4,7,11,10 > association_final.txt