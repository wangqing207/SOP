#!usr/bin/Rscript
###usage step1_methylation_preprocess.R 450k_annotation.txt
### input files include : 
##1: Methylation_EPIC (.idat) files 
##2: sample_sheet.csv 
##3: Annotation_EPIC.txt 
### output files : total.txt(include the annotation and betavalue),pheno(need group information) #######
args = commandArgs(T)
library(mime)
library(markdown)
library(beanplot)
library(illuminaio)
library(nor1mix)
library(siggenes)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(multtest)
library(splines)
library(stats)
library(base)
library(Biobase)
library(lattice)
library(reshape)
library(GenomicRanges)
library(IRanges)
library(XVector)
library(Biostrings)
library(foreach)
library(iterators)
library(locfit)
library(minfi)
library(bumphunter)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) 
#source("http://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
#### get the dir of filename ####
baseDir = getwd()
dir.create(paste(getwd(),"1_Data_QC",sep="/"))
dir.create(paste(getwd(),"2_Preprocess_data",sep="/"))
dir.create(paste(getwd(),"3_Dif_methylation",sep="/"))
######input data 

baseDir = getwd()
######RGset <- read.metharray.exp(file.path(baseDir, "200483200014"))
file.table = read.metharray.sheet(baseDir)  #读取csv工作表
file.table
exp<-read.metharray.exp(target = file.table,force=TRUE)
exp
pd <- pData(exp)
pd
#####################
###step 2:QC report 
## samples
a = qcReport(exp, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group, pdf = "1_Data_QC/QCReport.pdf",maxSamplesPerPage = 68)
###########
pdf("1_Data_QC/Densityplot.pdf")
par(cex=0.75)
c = densityPlot(exp,sampGroups = pd$Sample_Name,xlab = "Beta",pal = colorRampPalette(c("blue","green", "orange", "red"))(16))
print(c)
dev.off()
############
## QC-probesets
detP <- detectionP(exp)
failed <- detP>0.01
colMeans(failed) # Fraction of failed positions per sample
failed_probesets = sum(rowMeans(failed)>0.5) 
percent = failed_probesets/nrow(detP)
callrate = list(failed_probesets,(1-percent))
names(callrate) =c("failed_probesets","success_percent")
write.table(callrate,file = "1_Data_QC/probesetsQC.txt",sep = "\t",row.names = F,col.names = T)
###step 3 :preprocess
###reprocess1 These functions implements preprocessing for Illumina methylation microarrays as used in Genome Studio, the standard software provided by Illumina.
### qc_sample
MSet.raw <- preprocessRaw(exp)
qc_sample <- getQC(MSet.raw)
pdf("1_Data_QC/sampleQC.pdf")
k = plotQC(qc_sample)
print(k)
dev.off()
##########################
###preprocess2 Subset-quantile Within Array Normalisation (SWAN) for nfinium I and II type probes on a single array to be normalized together.
MSet.swan <- preprocessSWAN(exp, MSet.raw)
pdf("2_Preprocess_data/MDS.pdf")
c= mdsPlot(MSet.swan, numPositions = 866836,pch =2,sampGroups = pd$Sample_Group, sampNames = pd$Sample_Name)
print(c)
dev.off()
pdf("1_Data_QC/preprocess_compare.pdf")
par(mfrow=c(1,2))
d =  plotBetasByType(MSet.raw[,1], main = "Raw")
e =  plotBetasByType(MSet.swan[,1], main = "SWAN")
print(d)
print(e)
dev.off()
##########################
### get beta of each sites 
#methy=getMeth(MSet.raw)
#unmethy=getUnmeth(MSet.raw)
beta_value = getBeta(MSet.swan, type = "Illumina", offset=100, betaThreshold=0)

colnames(beta_value) = pd$Sample_Name
## get pheno and group informations 


##get pheno and group informations
pheno =as.data.frame(colnames(beta_value))
colnames(pheno) = "Samplename"
write.table(pheno,file = "pheno.txt",row.names = F,col.names = T)
### beta boxplot 
pdf("2_Preprocess_data/boxplot.pdf")#,width=10,height=10)
f = boxplot(beta_value,col = "red",boxwex = 0.1,main = "The boxplot of samples",xlab = "sample name",ylab = "beta value")
print(f)
dev.off()
#####sample cluster ######

d <- as.dist(dist(t(beta_value), method="euclidean"))
fit <- hclust(d, method="ward.D")
pdf("2_Preprocess_data/samplecluster.pdf")
a = plot(fit, main = "sample cluster",xlab = "sample",cex = 1)
print(a)
dev.off()
####
colnames(beta_value) =  paste(pd$Sample_Name,"AVG_Beta",sep =".")
detP <- detectionP(exp)
colnames(detP) = paste(pd$Sample_Name,"Detection.Pval",sep =".")
########### conbind the beta_value and p_values ##
beta_value = as.data.frame(beta_value)
detP = as.data.frame(detP)
betamatrix = matrix()

   for(i in c(1:ncol(beta_value))){
  betamatrix = cbind(betamatrix,beta_value[i],detP[i])
}
betamatrix = betamatrix[,-1]
###########

###########注释文件 MethylationEPIC_v-1-0_B2.csv :Annotation_EPIC.txt
annotation = read.delim(file = "Annotation_EPIC.txt",header = T)
#annotation = read.delim(file = "Annotation_EPIC.txt",header = T)
# dim(annotation)
#[1] 867531     47
Betamatrix = betamatrix[with(betamatrix,order(rownames(betamatrix))),]
#head(annotation[,c(1:3)])
Name<-rownames(Betamatrix)
Betamatrix_02<-cbind(Name,Betamatrix)
#rm(Betamatrix)
alldata<-merge(Betamatrix_02,annotation,by.x="Name")

#title0<-read.delim(file="total_titel.txt",sep="\t",header=T,quote=TRUE)  #列名改成全部大写，在excel表完成修改了，包括第5列"INFINIUM_DESIGN_TYPE"
#colnames(alldata)<-title0$title
# colnames(alldata)[28]<-"UCSC_REFGENE_NAME"
# colnames(alldata)[30]<-"UCSC_REFGENE_GROUP"
# colnames(alldata)[32]<-"RELATION_TO_UCSC_CPG_ISLAND"   #改850k
# colnames(alldata)[31]<-"UCSC_CPG_ISLANDS_NAME" 

write.table(alldata,file = "total.txt",sep = "\t",row.names = F,col.names = T, quote=FALSE)
###
## sed 's/\"//g' total1.txt > total.txt ## delete the ""  ###
###rm -rf total1.txt 
##获取甲基化和非甲基化值

#signal=cbind(Name,unmethy,methy,Betamatrix)
#write.table(signal,file = "Signal_unmethy_methy.txt",sep = "\t",row.names = F,col.names = T, quote=FALSE)
#methy_02=cbind(Name,methy)
#unmethy_02=cbind(Name,unmethy)
#write.table(methy_02,file = "methy.txt",sep = "\t",row.names = F,col.names = T, quote=FALSE)
#write.table(unmethy_02,file = "unmethy.txt",sep = "\t",row.names = F,col.names = T, quote=FALSE)
