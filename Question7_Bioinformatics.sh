### First using vcftools

### Pulling out the sites with a difference with the local one being the gzdiff sites file below
vcftools --gzvcf  /Chr21_not_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz –gzdiff /Chr21_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz --diff-site --out SItes_not_matching_chr21


### Now using bcftools

### First indexing the individual vcf files and then using bcftools isec to pull out files with overlap and unique variant
bcftools index -f /Chr21_not_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz
bcftools index -f /Chr21_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz
bcftools isec /Chr21_local/ ALL.chr21_GRCh38_sites.20170504.vcf.gz /Chr21_not_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz -p Ch21_bcftools_check


