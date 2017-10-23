#!/usr/bin/Rscript 
###usage  step3_IMA_sig.R  case_groupname   control_groupname 
rm(list=ls())
options(stringAsfactors = FALSE)
library(bioDist)
library(parallel)
library(KernSmooth) 
library(WriteXLS)
library(MASS)
library(base)
library(gmp)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(stats)
library(limma)
library(preprocessCore)
library(WriteXLS)
library(dplR)
library(IMA)
############################################################################################################
############################################################################################################
########################Prepare all the parameters here#####################################################
############################################################################################################
############################################################################################################

#######################options in IMA.methy450R#############################################################
######################Load data#############################################################################
libPaths = getwd()                            ### Specify the location of your R library
MethyFileName = "total.txt" ### Specfiy the original methylation data produced by the GenomeStudio
PhenoFileName = "pheno.txt"                 ### Specify the phenotype for each sample
############################################################################################################
args = commandArgs(T)
casegroup = args[1]
controlgroup = args[2]

###########output file######################################################################################
siteleveltest = paste(casegroup,"vs",controlgroup,"sigsite.txt",sep = "") ### Specify the path and name for the site-level testing result
list11excel = paste(casegroup,"vs",controlgroup,"sigregion.xls",sep = "")     ### Specify the path and name for region-level analysis results
list11Rdata = paste(casegroup,"vs",controlgroup,"sig.Rdata",sep = "")        ### Specify the path of Rdata file which stores the region-level analysis results
############################################################################################################

#################Preprocessing:IMA.methy450PP ##############################################################
samplefilterdetectP = 1e-5   ### The cutoff for sample-level detection Pvalue
samplefilterperc = 0.95      ### The percent of loci with detection Pvalue less than "samplefilterdetectP" in each sample
sitefilterdetectP = 0.05     ### The cutoff for site-level detection Pvalue
sitefilterperc = 0.05         ### The percent of samples with detection Pvalue less than "sitefilterdetectP" for each site
na.omit = TRUE               ### Remove the sites containing missing beta value
XYchrom = FALSE                ### Remove the sites on chromosome X
peakcorrection = FALSE       ### If TRUE, peak correction is performed
normalization = FALSE        ### If TRUE, quantile normalization performed
transfm = FALSE              ### If FALSE, no transform is performed; if "arcsinsqr", arcsin square root transformation is performed; if "logit", logit transformation is performed;
locidiff = FALSE             ### If FALSE, don't filter sites by the difference of group beta value. Otherwise, remove the sites with beta value difference smaller than the specified value
locidiffgroup = c(casegroup,controlgroup) ### Specify which two groups are considered to check the loci difference (if "locidiff" is not true)
snpfilter = FALSE            ### If FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site; otherwise, filter out the list of SNP containing loci by specifying the snp file name and location
### A list of SNP-containing probes (based on dbSNP v132) could be accessed by the command: snpfilter = system.file("extdata/snpsites.txt",package ="IMA")
##############################################################################################################

############sitetest/regionwrapper############################################################################
testmethod = "pooled"       ### Other options of differential testing methods: "wilcox"/"pooled"/"satterthwaite" for the comparison between two group
concov = "OFF"             ### If "ON", covariates is continuous variable
gcase = casegroup              ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = controlgroup           ### Specify the control group index in the sample.txt file (if "concov" is "ON")
Padj = "BH"                ### Options for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="median"        ### Options for deriving an index of overall methylation value of each region. mean/median/tbrm: "tbrm" is Tukey's Biweight robust average 
paired = FALSE             ### If ture, the differential test methods would change to the corresponding paired-test methods
####################################output the differential sites#############################################
rawpcut = 0.05             ### cut off for raw pvalue 
adjustpcut = NULL          ### cut off for adjusted pvalue
betadiffcut =0.14       ### cut off for beta value difference
##############################################################################################################

##############################################################################################################
###############################End of the parameter specification#############################################
##############################################################################################################


################################ Analysis Routes #############################################################
##############################################################################################################

.libPaths(libPaths) ## Specify your R library
library(IMA)        ## load the IMA package
data =IMA.methy450R(fileName = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName) ## load the data
dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,peakcorrection = peakcorrection,transfm = transfm,samplefilterdetectP = samplefilterdetectP,samplefilterperc = samplefilterperc,sitefilterdetectP = sitefilterdetectP,locidiff = locidiff, locidiffgroup = locidiffgroup,XYchrom = XYchrom,snpfilter = snpfilter) ## QC filtering

sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut,paired = paired) ## site-level testing with the "BH" adjustment
write.table(sitetest,file=siteleveltest,row.names=TRUE,sep = "\t") ## saving the results (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)

regionswrapper(dataf,indexmethod =indexmethod,gcase = gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,concov = concov,list11excel=list11excel,list11Rdata=list11Rdata,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut,paired = paired)   ## region-level testing for all 11 categories of annotated regions
############################### End of Analysis Routes #########################################################
################################################################################################################
