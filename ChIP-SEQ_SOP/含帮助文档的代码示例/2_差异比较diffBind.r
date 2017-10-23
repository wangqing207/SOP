args <- commandArgs(T)
f_name<- args[1] ####csvFile of the information 
if(!length(args==1)){
q()
}


library(DiffBind) 
#group1 <- 1
#group2 <- 2
######################################################################

tamoxifen = dba(sampleSheet=f_name,peakCaller="macs")
# tamoxifen = dba.count(tamoxifen,peaks="group_comparison_macs_peaks.xls",bCalledMasks=TRUE,minOverlap=1)
tamoxifen = dba.count(tamoxifen,minOverlap=0.4)


# tamoxifen = dba.contrast(tamoxifen, categories=3)
#tamoxifen = dba.contrast(tamoxifen, group1=group1, group2=group2)
tamoxifen = dba.contrast(tamoxifen, tamoxifen$masks$treatment, tamoxifen$masks$control)
tamoxifen = dba.analyze(tamoxifen,method=DBA_DESEQ)
tamoxifen.DB = dba.report(tamoxifen,method=DBA_DESEQ
                          ,th=1, bUsePval=FALSE, fold=0,bNormalized=TRUE
                          ,bCalled=T, bCounts=T,bCalledDetail=F
                          ,DataType=DBA_DATA_FRAME )
#######
#group_comp_tb <- read.table("group_comparison_macs_peaks.xls",sep="\t",header=T,stringsAsFactors=F)
#
#tamoxifen.DB <- merge(tamoxifen.DB,group_comp_tb[,c("chr","start","end","embDNA1","endDNA1")]
#                      ,by.x=c("chr","start","end")
#                      ,by.y=c("chr","start","end"),all=T,sort=F)
#
                          
tamoxifen.DB <- tamoxifen.DB[order(tamoxifen.DB[,"Chr"]),]                        
tamoxifen.DB_diff  <- tamoxifen.DB[tamoxifen.DB[,"p-value"]<=0.05,] 
tamoxifen.DB_hyper <- tamoxifen.DB[tamoxifen.DB[,"p-value"]<=0.05 & tamoxifen.DB[,"Fold"]>0,]  
tamoxifen.DB_hypo  <- tamoxifen.DB[tamoxifen.DB[,"p-value"]<=0.05 & tamoxifen.DB[,"Fold"]<0,]

colnames(tamoxifen.DB)[1:3] <- c("chr","start","end")
colnames(tamoxifen.DB_diff)[1:3] <- c("chr","start","end")
colnames(tamoxifen.DB_hyper)[1:3] <- c("chr","start","end")
colnames(tamoxifen.DB_hypo)[1:3] <- c("chr","start","end")
 

 
 
  

write.table(tamoxifen.DB,file=paste("total_DiffBind_",gsub(".csv","",f_name),".xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(tamoxifen.DB_diff,file=paste("diff_DiffBind_",gsub(".csv","",f_name),".xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(tamoxifen.DB_hyper,file=paste("hyper_DiffBind_",gsub(".csv","",f_name),".xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
write.table(tamoxifen.DB_hypo,file=paste("hypo_DiffBind_",gsub(".csv","",f_name),".xls",sep=""),sep="\t",col.names=T,row.names=F,quote=F)



