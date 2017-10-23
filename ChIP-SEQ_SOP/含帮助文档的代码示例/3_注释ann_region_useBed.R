args=commandArgs(T)
annotate_file=args[1] #bed file  :/data4/zhangjl/
gene_type=args[2] ###annotated type: CGI genetic tss 
file_to_ann=args[3] ###regrex subfix  : peaks.xls 
output_dir<-args[4] ### output dir 
annotate_dir<-args[5] ### the annotated existed dir 

if (!length(args)==5){
 print ("usage:
         Rscript ann_region_useBed.R ucsc.bed genetic/CGI peaks.xls/diffBind.xls output_dir dir_of_peaks_or_diffBind")
q()
}


library(GenomicRanges)
library(GenomicFeatures)



#gene_tb <- read.csv("mm10_ucsc_genetic.bed",sep="\t",header=F)
gene_tb <- read.csv(annotate_file,sep="\t",header=F)
## gene_tb <- read.csv("mm10_CGI_shore_shelf.bed",sep="\t",header=F)
## gene_tb <- read.csv("mm10_TSS_updown_2k5k.bed",sep="\t",header=F)



#gene_type <- "tss"  ## tss or genetic
## gene_type <- "genetic"
## bed must be 4 columns
setwd(annotate_dir)
temp<-list.files();
f_names<-temp[grep(paste("^[^an].+",file_to_ann,"$",sep=""),temp)]
#f_names <- list.files(pattern="^[^ag].+peaks.xls$")

#f_names<-temp[grep(paste(file_to_ann,"$",sep=""),temp)]

for(f_name in f_names){
print(paste("start",f_name,sep=":"))
######  read in differentially methylated sites or tiles
setwd(annotate_dir)
diff_df <- read.table(f_name,sep="\t",header=T, stringsAsFactors=F,check.names=F)
headDF<-names(diff_df)
headDF<-sub(pattern="stop",replace="end",headDF)
names(diff_df)<-headDF

gr_diff <-GRanges(seqnames=diff_df[,"chr"]
            ,ranges=IRanges(start=diff_df[,"start"],end=diff_df[,"end"])
		   # ,strand=diff_df[,"strand"]
           # ,name =diff_df[,"id"]
            )
######  gtf file to GRanges
gr_gene <-GRanges(seqnames=gene_tb[,"V1"]
            ,ranges=IRanges(start=gene_tb[,"V2"],end=gene_tb[,"V3"])
            ,name =gene_tb[,"V4"]
            )
######  find overlap
gr_diff_gene_overlap <- findOverlaps(gr_diff, gr_gene,type="any",select="all")
gr_diff_gene_tb1 <- cbind(diff_df[gr_diff_gene_overlap@queryHits,]
							,type=as.character(gr_gene@elementMetadata[gr_diff_gene_overlap@subjectHits,"name"]))
######  add not overlaped sites or tiles
gr_diff_gene_tb2 <- merge(diff_df,gr_diff_gene_tb1,by=colnames(diff_df),all=T,sort=F)
colnames(gr_diff_gene_tb2)[ncol(gr_diff_gene_tb2)] <- gene_type
######  output annotated table file
setwd(output_dir)
write.table(gr_diff_gene_tb2,file=paste("annotated",gene_type,f_name,sep="_"), sep="\t",col.names=T,row.names=F)


print(paste("finished",f_name,sep=":"))
############################################################################################
############################################################################################
#######################  annotation statistics
# f1_out <- file(paste("stats_",paste("annotated_",f_name,sep=""),sep=""),"w")
# ann_unique  <- unique(gr_diff_gene_tb2[,c("id","type")])
# ann_count <- table(ann_unique[,"type"],useNA="always")
# writeLines(paste(names(ann_count),collapse="\t"),f1_out)
# writeLines(paste(as.character(ann_count),collapse="\t"),f1_out)
# close(f1_out)
}

