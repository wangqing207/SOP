#!usr/bin/Rscript 
#### add annotation in site file ##
### usage :Rscript step4_null_annotation.R   Annotation_EPIC.txt
args = commandArgs(T)
annotation = args[1]
filedir = getwd()
b = list.files(path =filedir,pattern = "*site.txt$")
b
#########
## site1 for annotation and site2 for annotaion and plot #########
anno = read.delim(file = annotation,header = T)
anno=anno[order(anno$ILMNID),]

CpGsite = function(x){
input = x
data = read.delim(file = input,header = T)
class(data)
data1 = data[rownames(data) %in% anno$ILMNID, ]
anno1 =  anno[anno$ILMNID %in% rownames(data1),]
data2 = cbind(data1,anno1)
Target_ID = rownames(data2)
data2 = data.frame(Target_ID,data2)
attach(data2)
write.table(data2,file = paste(input,"null",sep ="."),sep = "\t",row.names = F,col.names = T)

 }

for(i in 1:length(b)){
  CpGsite(b[i])
}
