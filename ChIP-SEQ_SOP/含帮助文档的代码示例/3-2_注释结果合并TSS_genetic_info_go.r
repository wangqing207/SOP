args<-commandArgs(T)
file_to_ann_dir<-args[1]
refgene_fname<-args[2]
tax_id<-args[3]
type<-args[4]
NCBI_database_dir<-args[5]
print(file_to_ann_dir)
print(NCBI_database_dir)


if (! length(args==5)){
 print("usage:   /data4/sniu/programs/R-3.0.1/bin/Rscript  ../Rscript/TSS_genetic_info_go_zjl.R  /data5/zhangjl/Project/BC150047-1-chip/diffBind/myc_vs_IgG/ /data4/zhangjl/refData/hg19/hg19_refGene_uscs.bed  9606 genetic /data4/sniu/database/NCBI_gene/  ")
q()
}

macs_genetic_fnames <- list.files(file_to_ann_dir,pattern="^annotated_genetic_.+$")
f_names <- macs_genetic_fnames
print(f_names)





for(f_name in f_names){
print(paste("start:",f_name))
macs_genetic_fname <- f_name
macs_tss_fname <- gsub("genetic","tss",f_name)

## gene_info_fname <- "gene_info_10090"
## gene2refseq_fname <- "gene2refseq_10090"
## gene2go_fname <- "gene2go_10090"

## ncbi_fname    <- "GeneID_refseq_info_go_9606.txt"

###################################################################################
###############
###################################################################################
macs_genetic_tb <- read.table(macs_genetic_fname, sep="\t",header=T,stringsAsFactors=F,check.names=F)
macs_tss_tb <- read.table(macs_tss_fname, sep="\t",header=T,stringsAsFactors=F,check.names=F)

macs_tb <- merge(macs_tss_tb,macs_genetic_tb,by=colnames(macs_genetic_tb)[1:(ncol(macs_genetic_tb)-1)],all=T,sort=F )
macs_tb <- cbind(macs_tb,refseqID=sapply(macs_tb[,"genetic"],function(x){strsplit(x,"-")[[1]][1]}))
macs_tb <- unique(macs_tb)

refgene_tb1 <- read.table(refgene_fname, sep="\t",header=F,stringsAsFactors=F,check.names=F)
refgene_tb <- refgene_tb1[,c(1:4,6)] 
colnames(refgene_tb) <- c("refseq_chr","refseq_start","refseq_end","refseqID","refseq_strand")


setwd(NCBI_database_dir)

## ncbi_tb <- read.table(ncbi_fname, sep="\t",header=T,stringsAsFactors=F,check.names=F)
####################################
######### gene2refseq
gene2refseq_tb_1 <- read.csv(paste("gene2refseq",tax_id,sep="_"),sep="\t",header=F,stringsAsFactors=F)
gene2refseq_tb <- gene2refseq_tb_1[,c(2,4)] 
gene2refseq_tb[,2] <- sapply(gene2refseq_tb[,2],function(x){strsplit(x,"\\.")[[1]][1]})
gene2refseq_tb <- unique(gene2refseq_tb)
gene2refseq_tb <- gene2refseq_tb[gene2refseq_tb[,2]!="-",]
colnames(gene2refseq_tb) <- c("GeneID","refseqID")
######  gene_info
gene_info_tb_1   <- read.csv(paste("gene_info",tax_id,sep="_"),sep="\t",header=F,stringsAsFactors=F) 
gene_info_tb <- unique(gene_info_tb_1[,c(2,3,9)]) 
colnames(gene_info_tb) <- c("GeneID","symbol","description") 
###### gene2go
gene2go_tb_1     <- read.csv(paste("gene2go",tax_id,sep="_"),sep="\t",header=F,stringsAsFactors=F)  
gene2go_tb <-  unique(gene2go_tb_1[,c(2,3,6,8)])
colnames(gene2go_tb) <- c("GeneID","GO_ID","GO_term","GO_category")

######################
gene_refseq_info_tb <- merge(gene2refseq_tb,gene_info_tb,by.x="GeneID",by.y="GeneID",all.x=T,all.y=F,sort=F)

###################
macs_refgene_tb1 <- cbind(macs_tb,refgene_tb[ match(macs_tb[,"refseqID"],refgene_tb[,"refseqID"])  ,]) 
macs_refgene_tb1 <- macs_refgene_tb1[,-grep("refseqID",colnames(macs_refgene_tb1))[2]]
macs_refgene_tb1[,"refseqID"] <- as.character(macs_refgene_tb1[,"refseqID"])
###################
macs_refgene_info_tb1 <- merge(macs_refgene_tb1,gene_refseq_info_tb,by="refseqID",all.x=T,all.y=F,sort=F) 
macs_refgene_info_tb2 <- macs_refgene_info_tb1[,c(colnames(macs_refgene_tb1), colnames(gene_info_tb))]
##################
gene2go_tb_exist <- gene2go_tb[!is.na( match(gene2go_tb[,"GeneID"], macs_refgene_info_tb1[,"GeneID"])),]

macs_refgene_info_go_tb1 <- merge(macs_refgene_info_tb1,gene2go_tb_exist,by="GeneID",all.x=T,all.y=F,sort=F)
macs_refgene_info_go_tb2 <- macs_refgene_info_go_tb1[,c(colnames(macs_refgene_tb1), colnames(gene_info_tb),colnames(gene2go_tb)[-1])]

setwd(file_to_ann_dir)
############################################################################################
## with GO annotation
write.table(macs_refgene_info_go_tb2 ,file=paste("annotated_tss_genetic_GO_",  gsub("annotated_genetic_","",macs_genetic_fname),sep="")
        ,sep="\t",col.names=T,row.names=F)
## without GO annotation
write.table(macs_refgene_info_tb2,file=paste("annotated_tss_genetic_",  gsub("annotated_genetic_","",macs_genetic_fname),sep="")
        ,sep="\t",col.names=T,row.names=F)
## if no gene annotation
#write.table(macs_tb,file=paste("annotated_tss_genetic_",  gsub("annotated_genetic_","",macs_genetic_fname),sep="")
#        ,sep="\t",col.names=T,row.names=F)        
        
########
print(paste("finished:",f_name))

}
