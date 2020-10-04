setwd("~/Documents/GitHub/challenge/data/Bioinfo_Que16")


related_samples=read.table(file="relateds.tsv",header = T,sep="\t",stringsAsFactors = F)
pheno_samples=read.table(file="relateds.pheno.tsv",header = T,sep="\t",stringsAsFactors = F)


all_related_samples=c(related_samples$IID1,related_samples$IID2)


all_related_samples=unique(all_related_samples)


all_related_samples=as.data.frame(all_related_samples)

library(dplyr)

pheno_samples_those_related=left_join(all_related_samples,pheno_samples,by=c("all_related_samples"="IID"))


to_exclude=subset(pheno_samples_those_related,pheno_samples_those_related$pheno==0)
to_not_exclude=subset(pheno_samples_those_related,pheno_samples_those_related$pheno==1)


length(unique(to_not_exclude$all_related_samples,related_samples$IID1))
length(unique(to_not_exclude$all_related_samples,related_samples$IID2))


## to not exclude has related individuals with the disease
pheno_cases=subset(pheno_samples,pheno_samples)


pheno_samples$status=ifelse(pheno_samples$pheno==0,"control","case")

pheno_samples=pheno_samples[,c(1,3)]
related_samples_stats1=left_join(related_samples,pheno_samples,by=c("IID1"="IID"))
colnames(related_samples_stats1)[3]="IID1_status"
related_samples_stats2=left_join(related_samples_stats1,pheno_samples,by=c("IID2"="IID"))  
colnames(related_samples_stats2)[4]="IID2_status"


controls_to_remove=subset(related_samples_stats2,related_samples_stats2$IID1_status=="control" & related_samples_stats2$IID2_status=="control")


control_samples=unique(c(controls_to_remove$IID1,controls_to_remove$IID2))

IID1_unique=unique(controls_to_remove$IID1)
IID2_unique=unique(controls_to_remove$IID2)


length(intersect(IID1_unique,IID2_unique))


IID2_ids=unique(controls_to_remove$IID2)
vector_of_removals=character(1)

for(row in 1:nrow(controls_to_remove)){
  IID1_samples=controls_to_remove[row,1]
  if(IID1_samples %in% IID2_ids==TRUE) {
    #print("found a match")
    print(IID1_samples)
    vector_of_removals = c(vector_of_removals,IID1_samples)
  }
}

final_removal_list=c(vector_of_removals,IID2_unique)

write.table(final_removal_list,file="Control_samples_to_remove.txt",col.names = F,row.names = F,sep="\t",quote = F)

