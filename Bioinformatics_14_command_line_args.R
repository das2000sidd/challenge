setwd("/Users/sdi0596/Documents/GitHub/challenge")


args <- commandArgs()


pedfile <-sub('--pedfile=','',args[grep('--pedfile=',args)])
mapfile <-sub('--mapfile=','',args[grep('--mapfile=',args)])

pedfile
mapfile
ped=read.delim(file=pedfile,header = F,sep="\t",stringsAsFactors = F)
map=read.delim(file=mapfile,header = F,sep="\t",stringsAsFactors = F)
row=1
first_ped_allele_column=7
ped_rows=nrow(ped)
ped_cols=ncol(ped)
map_rows=nrow(map)
map_cols=ncol(map)
table_vcf_format=matrix(0,nrow=ped_rows,ncol=map_rows+1) ## table for desired output
table_vcf_format=as.data.frame(table_vcf_format)

map$snp="genotype_SNP"
map$id=1:nrow(map)
map$snp_id=paste(map$snp,map$id,sep="") ## giving unique names to snps
colnames(table_vcf_format)[1]="sample_name"
colnames(table_vcf_format)[2:ncol(table_vcf_format)]=map$snp_id

ped_matrix=as.matrix(ped)
col_desired_output=2
no_of_times_run=(ncol(ped)-6)/2
run=1
for(index_row in 1:ped_rows){
  table_vcf_format[index_row,1]=ped_matrix[index_row,1]
  first_ped_allele_column=7 ## needs to be set to 7 for a new sample since that is where the first genotype starts
  for(run in 1:no_of_times_run){
    first_allele=ped_matrix[index_row,first_ped_allele_column] ## first_ped_allele_column is for extracting the allele from ped file
    second_allele=ped_matrix[index_row,first_ped_allele_column+1]
    genotype=paste(first_allele,second_allele,sep="/")
    print(genotype)
    table_vcf_format[index_row,run+1]=genotype
    col_desired_output=col_desired_output+1 ## col_desired_output is for controlling in which column of final file the genotype goes
    first_ped_allele_column=first_ped_allele_column+2 ## first_ped_allele_column incremented by two since the the first allele of next genotype in 2 columns away from first allele of current genotype
    
  }
}






write.table(table_vcf_format,file="Ped_in_Vcf_format.txt",col.names = T,row.names = F,sep="\t",quote = F)



## ran using 

Rscript ./Bioinformatics_14_command_line_args.R --pedfile=ped_subset.ped --mapfile=map_subset.map



