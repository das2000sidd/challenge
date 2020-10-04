## To remove duplicated positions
gunzip -c duplicates.vcf.gz | grep -v "#" | cut -f1,2 |  sort -k1,1 -k2,2n | uniq -d > duplicated_loci.txt
vcftools --gzvcf duplicates.vcf.gz --exclude-positions duplicated_loci.txt --recode --out unique_variants

### To remove duplicate positions with duplicate alleles
gunzip -c duplicates.vcf.gz | grep -v "#" | cut -f1,2,4,5 | sort -k1,1n -k2,2n -k3,3 -k4,4 | uniq -d > duplicated_loci_both_allele.txt
cut -f1-2 duplicated_loci_both_allele.txt > duplicated_loci_both_allele_only_locus.txt
vcftools --gzvcf duplicates.vcf.gz --exclude-positions duplicated_loci_both_allele_only_locus.txt --recode --out unique_no_duplicate_loci



### To remove where  loci and allele both are duplicated
gunzip -c duplicates.vcf.gz | grep -v "#" | cut -f1,2,4 | sort | uniq -d > duplicated_loci_ref_allele.txt
cut -f1-2 duplicated_loci_ref_allele.txt > duplicated_loci_ref_allele_to_remove.txt
vcftools --gzvcf duplicates.vcf.gz --exclude-positions duplicated_loci_ref_allele_to_remove.txt --recode --out unique_no_duplicate_loci_and_allele


