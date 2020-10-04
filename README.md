# ITG Challenge

## Description/Instructions

Please find below a few questions that mimic some common problems we encounter at ITG. They are grouped by broad theme. You will notice there are many questions, the goal is not to answer them all but to pick a few questions to focus on (10 is a good number, but pick as many as you want). You should pick from all three categories, but there are many more bioinformatics questions so you should mainly pick from those. We encourage you to choose your questions according to your areas of expertise but also to try and answer questions that are as varied as possible.

For programmatic questions, you can use the language and libraries of your choice, but we will assess whether your choice of language was optimal. Try and aim for a minimal solution in terms of code length. If you use a shell script, you can assume that common non-core packages will be installed (e.g. `awk`, `sed`, `perl`, `python`, `sponge`, `wget` or `jq`). You can use the shell of your choice, if not otherwise specified we will assume `bash`. Assume that all common bioinformatics tools `bcftools`, `bedtools`, `vcftools`, `plink` and others are all installed.

We are primarily interested in how you would solve these problems if you encountered them in real life. Whenever the command line or programming is used, please include your code along with your answer. Not all questions are programmatic, and some are open-ended. Feel free to include details and to discuss potential issues if you don't think there is a clear-cut answer.

To submit your results, please clone this repository and make your edits. Once you're done, send us a link to your repository, or compress it and send us a link to the archive.

## Questions

### Support/resource management/Shell
1. A user has several versions of R installed in their path. Each version of R has a number of locally installed libraries. The user is confused, and would like to know which library is installed for each version of R. Can you write a command to help them out?

Following is the command I used to determine the libraries installed for the various R version in our institutional HPC cluster:

ls -ltrh /home/sdi0596/R/x86_64-pc-linux-gnu-library/*


Following is a sample output for two different version of R:

/home/sdi0596/R/x86_64-pc-linux-gnu-library/4.0:
total 54K
drwxrwxr-x  7 sdi0596 sdi0596 4.0K Jul 31 12:06 bitops
drwxrwxr-x  7 sdi0596 sdi0596 4.0K Jul 31 12:07 caTools
drwxrwxr-x  8 sdi0596 sdi0596 4.0K Jul 31 12:10 utf8

/home/sdi0596/R/x86_64-pc-linux-gnu-library/3.2:
total 5.5K
drwxrwxr-x  7 sdi0596 sdi0596 4.0K Dec 12  2018 BSgenome.Mmusculus.UCSC.mm10
drwxrwxr-x  7 sdi0596 sdi0596 4.0K Dec 12  2018 BSgenome.Mmusculus.UCSC.mm9


3. A user wants to install an `R` package and gets the following [error log](data/error.log). What is likely to cause the error and how can they solve it?

The most important error message in that error log is :

unrecognised command line option ‘-std=c++11’

This prevents compilation of the subsequent packages being installed and is due to an old gcc compiler being used which is not compatible with the version of compiler required for the tools being used. 
A possible resolution would to update to the latest version of the gcc compiler.

4. A user is running commands like this one `cat file1 <(cut -d " " -f 1-15,17,18 file2) > file3`. What does this command do? It runs fine on the command line, but then the user includes it into a file with other commands, saves it and runs `chmod +x` on it. However, that line of code throws the following error : `syntax error near unexpected token '('`. What has the user forgotten?

This command appends file2 to the end of file1.

The user gets the error because he has forgotten to escape the brackets with a backslash.



5. A collaborator has sent you [this script](data/EasyQCWrapper.sh). It is a wrapper for a bioinformatics software called `EasyQC`.  Running it, you get the following error: 

    ```bash
    ./test.EasyQC-START.R: line 6: syntax error near unexpected token 'EasyQC'
    ./test.EasyQC-START.R: line 6: 'library(EasyQC)'
    ```

     You need to run this script now, but your collaborator is unavailable for a few days. What is causing the error? (Hint: Nothing is wrong with the `.ecf` EasyQC script.)
     
The library call is a R language specific syntax. It is being attempted to run as bash script and hence this R error is being thrown. 

To avoid this error, it needs to run using the “Rscript” command

     

7. Bioinformaticians often work on a computing cluster. The cluster runs a software called a job scheduler that attributes resources to users depending on the requirements of their jobs. In this case, let's imagine the cluster is running IBM LSF. You do not need to know it to answer this question. The `bjobs` command lists all jobs submitted by the user (manpage [here](https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.2/lsf_command_ref/bjobs.1.html)). It gives this kind of output:
    ```
    JOBID   USER             STAT  QUEUE      FROM_HOST EXEC_HOST JOB_NAME SUBMIT_TIME
    9670681 current_user     RUN   basement   head_node node1     job1     Oct 24 10:24
    9740051 current_user     RUN   basement   head_node node1     job2     Oct 24 17:41
    9670681 current_user     RUN   normal     head_node node2     job3     Oct 24 10:24
    9740981 current_user     PEND  basement   head_node           job4     Oct 24 17:44

    ```
     - Given the [following output](data/farm-snapshot.txt) of `bjobs -all`, which users are the top 5 users of the cluster?
     - How many jobs does the user `pathpip` have running in all queues?
     - A user wants to know how many jobs they have pending (`PEND`) and running (`RUN`) in each queue. Write a command line to do that (You can use the log above to check your command line). How would they display this on their screen permanently, in real time?
     
     Using the following command, a list of number of jobs per used was determined in descending order:
cat farm-snapshot.txt | awk '{print $2}' | sort | uniq -c | sort -k1,1nr


Based on the output above, the top 5 users are:

7146	km18
3475	ro4
2321	igs
1655	nw17
1521	pathpip



The user papthpip has a total of 1521 jobs running. Following is the command used to determine that value:

cat farm-snapshot.txt | grep "pathpip" | wc -l

To get the number of pending jobs, following is the command:

grep "PEND" farm-snapshot.txt | wc -l

The number of pending jobs was 13336.


To get the number of running jobs, following is the command:

grep "RUN" farm-snapshot.txt | wc -l

The number of running jobs was 5580.


9. All major computational tasks in your lab are done via SSH connection to mainframe servers or HPC clusters. A user comes from a Linux (mostly command-line) background but IT only support Windows 10 for laptops. How would you advise them to configure their laptop to make their transition easier?


A suitable advice would be to download and install a SSH client for windows such as PuTTY or WinSCP. If configured correctly, they should be able to SSH into the HPC clusters with their institutional login id and password.

### Bioinformatics

2. From an existing VCF with an arbitrary number of samples, how do you produce a VCF file without any samples using `bcftools`?

The logic here was first generate the unique list of samples for a vcf file using bcftools query and then remove them using bcftools view with -S flag by appending a “^” before the sample names.

bcftools query -l query_file.vcf > samples_to_remove.txt

awk '{print "^"$0}' samples_to_remove.txt > samples_to_remove_with_symbol.txt

bcftools view -S samples_to_remove_with_symbol.txt query_file.vcf > query_file_no_samples.vcf


Here it is assumed query_file.vcf is a multi-sample VCF with at least two samples.


4. How do you convert a gzipped VCF to the `bimbam` format? (you may choose to script a solution yourself, or not)

A gzipped VCF can be converted to a bimbam file using the following code. This was run using Plink 1.9. A sample format would be as shown below:

plink –vcf vcf_file –snps-only –recode-bimbam –out vcf_bimbam


6. How would you change the chromosome numbers in the file above to chromosome names (e.g. "chr1" instead of "1")?
    - How would you change the names back to the original? Would your solution work if an additional column containing text of arbitrary length and content is appended at the left of the file?
    - These positions are extracted from a VCF. Convert this file to the BED format.
    
Following is the command use to add “chr” in front of pre existing chromosome name:
 
awk 'OFS="\t" {$1="chr"$1; print}' rand.chrpos.txt > rand.chrpos.with.chr.word.txt

To remove chr word and return back to original, following can be done:

sed 's/chr//g' rand.chrpos.with.chr.word.txt > no.chr.word.txt

 To convert VCF file to bed file, following is the command:

gunzip -c duplicates.vcf.gz | grep -v "#" | awk '{FS="\t";OFS="\t";print $1,$2,$2+length($5)-length($1),$4,$5}' > duplicates.bed

In the code above $2+length($5)-length($1) is to determine the length of deletion in case there is one. 
   
    
    
7.	Download the 1000 Genomes sites VCF file for chromosome 21 [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr21_GRCh38_sites.20170504.vcf.gz). We want to compare it to [a locally stored file](data/ALL.chr21_GRCh38_sites.20170504.vcf.gz).
    - What is the fastest way to check the integrity of, or compare, any such downloaded file?
    - If you find that the files are indeed different, how do you find their differences? Keep in mind that this kind of file can be very large (>100Gb), your solution should be as fast and memory-efficient as possible.
    - If you found no differences in the end, what could cause a false alarm?
    
    The fastest way to check the integrity of the downloaded file and locally stored file would be to compare the MD5 sum of both. It is a 32 character hash value unique to the file and will not match if the contents of the two files do not match.

There are multiple ways to find the differences using vcftools, bcftools as well as bedtools.
For vcftools, following is a sample command:

vcftools --gzvcf  ALL.chr21_GRCh38_sites.20170504.vcf.gz –gzdiff ALL.chr21_GRCh38_sites.20170504.vcf.gz --diff-site --out SItes_not_matching_chr21

If there is a difference in the files, then the common unique sites of each file will be written SItes_not_matching_chr21.recocde,.vcf

For bcftools, following is a sample set of commands:

bcftools index -f / Chr21_not_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz 
bcftools index -f / Chr21_local/ ALL.chr21_GRCh38_sites.20170504.vcf.gz
bcftools isec / Chr21_local/ ALL.chr21_GRCh38_sites.20170504.vcf.gz / Chr21_not_local/ALL.chr21_GRCh38_sites.20170504.vcf.gz -p Ch21_bcftools_check

From the above command, three files were generated named 0000.vcf, 0001.vcf and 0002.vcf. The unique variants go to 0000.vcf and 0001.vcf while the common variants go to 0002.vcf. In this case all the variants went to 0002.vcf suggesting there wasn’t a difference in the VCF files.

    
   
8.	What is the p-value corresponding to standard normal z-scores of 10.35, 29.7, 45.688 and 78.1479?

Following are the p values corresponding to the standard normal distributions:

2*pnorm(-abs(10.35))= 4.184858e-25
2*pnorm(-abs(29.7))= 7.678615e-194
2*pnorm(-abs(45.688)) = 0
2*pnorm(-abs(78.1479)) = 0

Here it is assumed that a 2 sided p value has been desired. If not, then the multiplication by 2 has to be removed.



9.	We want to round a column of numbers to `n` decimal places, with values with 5 as their rightmost significant digit rounded up. Use the language of your choice.

Assuming the language is R and a vector of values called ‘numbers’ is available, following is the command:

round(numbers,digits = 5)



10.  Is [this HRC-imputed file](https://drive.google.com/open?id=1dOYlwIlAyz9-i4gVy2sgpQv_4YX9Y5nU) missing any chromosomes? Try to find out in seconds if you can.

Following was the command using to determine total number of unique chromosomes:

gunzip -c hrc.positions.txt.bgz.gz | cut -f2 | sort | uniq | wc -l

The command is unzipping the file and  piping the column 2 from this file into the sort command and subsequently the unique command is pulling out the unique chromosome ids.
There was a total of 20 chromosomes. Since this is from HRC which is human data, we can say that chromosome X and Y are missing.


12. How would you convert a VCF file to the Plink binary format? How would you do the reverse, and what kind of problems do you anticipate?

Following is a sample code to convert vcf file to Plink binary format using a combination of vcftools and plink.

vcftools --gzvcf duplicates.vcf.gz --plink-tped --out duplicates_plink
plink --tfile duplicates_plink --make-bed --out duplicates_plink_binary


To convert plink binary file to vcf file, following would be a sample command:

plink --bfile bfile_generated --recode vcf --out vcf_from_bfile

One of the problems in the above command is that Plink version 1.9 may truncate indels that are too long and hence the number of variants in the Plink format map files and ped files may not match with that of VCF file.



13. Write a snippet to reformat a PED file so as to output a file with the following header `sample_name genotype_SNP1 genotype_SNP2 ...` where genotypes are coded VCF-style (e.g `A/C`, the order of the alleles in the output is not important).

A R code was used for this.

Following is the syntax to run it:

Rscript ./Bioinformatics_14_command_line_args.R --pedfile=ped_subset.ped --mapfile=map_subset.map


Where ped_subset.ped and map_subset.map should be plink generated Ped and Map file.
It will generate a file named “Ped_in_Vcf_format.txt



14. A genetic association pipeline starts with a VCF and produces summary statistics for every position. The VCF contains multiallelics and indels. Unfortunately, a program in the pipeline trims all alleles to their first character. Why might allele frequencies not always be equal for a given variant? Find a way to correct the alleles in the association file by using the information from the VCF. Select columns are provided for [the association file](https://github.com/hmgu-itg/challenge/raw/master/data/association.txt.gz). We also provide [a file](https://github.com/hmgu-itg/challenge/raw/master/data/fromVCF.txt.gz) that was created from the VCF using `bcftools query -f '%CHROM %POS %REF %ALT %AN %AC\n'`.

Allele frequencies may not be equal for a given variant due to several reasons. 
1.	If one of the alleles happens to be deleterious in a homozygous genotype, then natural selection will work against the selection of that allele.
2.	Another reason could be a founder effect where due to inbreeding, a large section of the population is homozygous in genotype for a particular allele and hence one allele has a disproportionately greater allele frequency than the other.
3.	A third possible reason could be that the variant has rose recently in the population and has not had a history long enough for it to have a frequency close enough to the major allele at that locus.

awk -F'\t' '{print "chr"$0}' association.txt > association_with_chr_word.txt
awk -F'\t' 'BEGIN{OFS="\t";} {print $0,$1"_"$3}' association _with_chr_word.txt | sort -k7 > association _chr_position_combined.txt
awk -F'\t' 'BEGIN{OFS="\t";} {print $0,$1"_"$2}' fromVCF.txt | sort -k7 > fromVCF_chr_position_combined.txt


join -1 7 -2 7 -t $'\t' association _chr_position_combined.txt fromVCF_chr_position_combined.txt | cut -f2,3,4,7,11,10 > association_final.txt

One of the caveats of this approach is that with the VCF file reports multiple alleles for indels and thus the joined file 'association_ final.txt' has more lines than original association file

Eg: for position 10629861 on chr10, the association file has the following positions:
10	chr10:10629861[b38]	10629861	C	C	0.015

But the VCF file has the following:
chr10	10629861	C	CAAA	2616	46
chr10	10629861	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	2616	20
chr10	10629861	C	CAAAAA	2616	24
chr10	10629861	C	CAAAAAAAAAAAAAAAAAAAAA	2616	18
chr10	10629861	C	CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA	2616	13
chr10	10629861	C	CAA	2616	19



16. A researcher wants to conduct a disease association study. However, colleagues warn him that the dataset contains related individuals. He would like to remove relatedness in his dataset, but given his disease is rare, he would also like to maximise the number of cases kept in. Using [a list of samples with disease status](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.pheno.tsv) and [a file containing pairs of individuals above a relatedness threshold](https://github.com/hmgu-itg/challenge/raw/master/data/relateds.tsv), create an exclusion list of samples to remove to help the researcher achieve their goal.





### Statistical genetics

3. A common practice when performing genetic association studies is to perform an ethnicity check as a quality control step. Can you explain how this is done?
    - You are dealing with whole-genome sequencing data in 2,326 Bulgarian samples. How would you perform such an ethnicity check, and which projection dataset would you use? 
    
    Ethnicity check is done to make sure cases and controls are not from different populations whereby a spurious association may arise simply because one population has certain variants and the other does not due to underlying differences. 
Ethnicity check is done using principal component analysis. Usually representative control population data from the thousand genomes data is combined together in plink format file with the project analysis data. Subsequently the data is projected on to the first two principal components since they account for the highest and second highest amount of variation in the data. Here it would be expected if there is not an issue of ethnicity, then the study data (especially the unaffected controls) will cluster with the control data from it’s own ethnicity. If there is a population stratification issue, then the controls from the study will not cluster with the control population data.

For the Bulgarian samples, since it is an Eastern European country, I would use the EUR dataset of the thousand genomes project. 

In order to this analysis, I would use plink. I would generate a binary PED format file for the Bulgarian samples and also for the thousand genomes EUR data. Assuming these would VCF files, I would use vcftools to generate Plink compatible map and PED file. Then using Plink I would first remove SNPs with more than 10% of samples having no data and only use autosomal SNPs. Subsequently I would use a maf filter of 0.05 and Hardy Weinberg equilibrium violating SNPs to remove certain SNPs.  Then I would combine the plink binary files for both my data and the EUR data and generate a N*N matrix of genomewide pairwise distances matrix using plink –mds-plot flag with 4 principal components. Subsequently I would plot values of the positions on the first and second dimension determined by the MDS analysis and see if all the individuals of interest are clustering with the EUR population or not. If it is, we will be able to conclude that the data is not showing population stratification under control conditions.

    
    
    
4. You perform a single-point association with blood lipids and find a variant with MAF=0.7% associated at p=1e-14. Do you expect the effect size to be large or small? What would be your next steps investigating this signal?

Based on the observation that common variants tend to have smaller effect sizes and rare variants larger, I would expect this rare variant (since minor allele frequency is < 1%) to have a large effect size. 
Once the variant loci is found, I would annotate this variant to see where in the genome it lies. If it causes a change in the coding sequence, then initial conservation check could be done to see how conserved that particular position is across species and that would give us some direction to how strong the effect of that mutation could be.
If it is not in a coding sequence, then I would see if it is in other important regulatory regions of a gene such as a splice site, promoter or an enhance region. If so, then cis-eQTL analysis can be done to see if this SNP effects gene expression.



6. An analyst studies a population of remote villages in Eastern Europe. They are interested in a particular variant, and compare the frequency in their villages (3.5%) to the EUR population frequency in the 1000 Genomes (0.03%). They conclude that the variant has increased in frequency in their villages. Do you agree, and if not, what would your advice be?

The first step in this study would be to do a principal component analysis of the 1000 genomes EUR data and the samples from the remote villages in Europe to see if the remote villages cluster with the EUR population. If they do not, then the right population group from the thousand genomes need to be investigated to see which is the population group these isolated villages are close to the most.
My assumption would be these samples would not cluster well with most of the thousand genomes. This is probably a founder effect whereby an isolated population arose from a small number of individuals who settled in those remote villages a long time ago and have not outbred leading to a large frequency of individuals carrying the rare variant.


7.  The same analyst sends you association summary statistics for random glucose.
    - Which checks would you perform on such a dataset?
    - You wish to include this dataset in a meta-analysis. Which additional information should you ask for in your next email to your colleague?
    - In this dataset, you observe  &#955;=1.25. The analyst has adjusted for age, age-squared, sex, and fasting status. What would you suggest they do?
    Following are the checks than can be performed:

Following are the checks to perform:

1.	The p values for the individual SNP have been adjusted for multiple testing correction or not
2.	SNPs with MAF less than 0.05 have been removed or not
3.	Individuals with genotyping type less than 95% for all the SNPs on the array have been excluded
4.	The genomic inflation factor which should be around 1 suggesting the statistics reported for the SNPs will not be inflated.



Based on that genomic inflation factor value, I would advice them to divide the chi square test statistic for each variant by that value which would make it less likely to observe extreme p values.


    
8. You are a co-author on a manuscript. You receive a draft, in which the main author has used the traditional &#945;=5e-8 as the significance threshold. The paper describes an analysis of 10 related blood phenotypes (platelet count, platelet volume, immature platelet fraction ...) using the fixed markers of the Infinium ImmunoArray on 897 individuals. What do you think about the chosen threshold, and what would you suggest to the first author? What would be your comments on the design of the study? 


The Infinium Immunoarray from Illumina has about 250,000 fixed markers as per Illumina documentation. Assuming only the fixed markers are used for the analysis, they are all tag SNPs and the p value threshold is set at 0.05 for an individual SNP, a Bonferroni correction p value would be 0.05/250,000 which equals to 2e-7. Thus, the threshold of 5e-8 used the main author would be quite strict and could miss out on some true associations. Thus, my advice would be to lower the significance threshold to 2e-7.


