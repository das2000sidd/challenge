### To convery VCF file to Plink binary format file using a combination of vcftools and plink
vcftools --gzvcf duplicates.vcf.gz --plink-tped --out duplicates_plink 
plink --tfile duplicates_plink --make-bed --out duplicates_plink_binary


## To convert plink binary file to vcf file
plink --bfile duplicates_plink_binary --recode vcf --out vcf_from_bfile


### HERE ASSUMING THAT duplicates.vcf.gz is a vcf file with sample genotype info