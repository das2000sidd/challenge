bcftools query -l query_file.vcf > samples_to_remove.txt
awk '{print "^"$0}' samples_to_remove.txt > samples_to_remove_with_symbol.txt
bcftools view -S samples_to_remove_with_symbol.txt query_file.vcf > query_file_no_samples.vcf