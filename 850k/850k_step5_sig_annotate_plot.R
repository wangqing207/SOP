#!usr/bin/Rscript  
######## program for format the CpDsite.txt ###################
##usage Rscript 850k_step5_sig_annotate_plot.R
### input files: chromInfo_human_UCSC.bed  Annotation_EPIC.txt ##
### output files : xxsigsite.txt   xxsigsite_plot.pdf in Dif_methylation ######
filedir = getwd()
b = list.files(path =filedir,pattern = "*sigsite.txt$")
#b
#########
#### site1 for annotation and site2 for annotaion and plot #########
anno = read.delim(file = "Annotation_EPIC.txt",header = T)
anno=anno[order(anno$ILMNID),]
CpGsite = function(x){
input = x
data = read.delim(file = input,header = T)
data1 = data[rownames(data) %in% anno$ILMNID, ]
anno1 =  anno[anno$ILMNID %in% rownames(data1),]
data2 = cbind(data1,anno1)
Target_ID = rownames(data2)
data_all = data.frame(Target_ID,data2)
attach(data2)
write.table(data_all,file = input,sep = "\t",row.names = F,col.names = T)
}

for(i in 1:length(b)){
  CpGsite(b[i])
}

b = list.files(path =filedir,pattern = "*sigsite.txt$")
b
CpGsite = function(x){
input = x
data2 = read.delim(file = input,header = T)#################
attach(data2)
data3 = data.frame(CHR,MAPINFO,Beta.Difference)
end = data3[,2] +1 
data4 = data.frame(data3[,c(1:2)],end,data3[,3])
colnames(data4) = c("chr","start","end","dataponit")
a = paste("chr",data4[,1],sep ="")
data4[,1] = a 
data4_red = data4[data4$dataponit > 0,]
dim(data4_red)
value1 = c(rep("red",nrow(data4_red)))
data5_1 = data.frame(data4_red,value1)
colnames(data5_1) = c("chr","start","end","dataponit","value")
data4_green =data4[data4$dataponit < 0,]
dim(data4_green)
value2 = c(rep("green",nrow(data4_green)))
data5_2 = data.frame(data4_green,value2)
colnames(data5_2) = colnames(data5_1)
data5 = rbind(data5_1,data5_2)
result_name = unlist(strsplit(input,"\\."))[1]
write.table(data5,file = paste(result_name,"plot",".txt",sep=""),sep = "\t")
}
for(i in 1:length(b)){
  CpGsite(b[i])
}
################
################################################################
################################################################
drawPics=function (x) {
  chr_fname <- "chromInfo_human_UCSC.bed"
  region_fname <- x
  chrs <- paste("chr",c("X","Y",22:1),sep="")

################################################################
################################################################
  chr_tb1 <- read.table(chr_fname,sep="\t",header=F)
  chr_tb2 <- chr_tb1[grep("^chr",chr_tb1[,1]),]
  max_chr_len <- max(chr_tb2[,2]) 
  max_plot_len <- 10e6*ceiling(max_chr_len/10e6)
  pdf_name = unlist(strsplit(x,"\\."))[1]
  pdf(file=paste(pdf_name,".pdf",sep=""));
  
  
  plot(c(1, max_plot_len), c(1, nrow(chr_tb2)), type= "n",xlab="Position (M)", ylab="Chromosome",axes=F)
  par(mai=c(1,0,0.9,1))
  for( i in 1:nrow(chr_tb2)){
    rect(xleft=1,xright=chr_tb2[i,2]+1,ytop=match(chr_tb2[i,1],chrs)-0.1,ybottom=match(chr_tb2[i,1],chrs) +0.1,col="gray")
  }

  text(rep(-5e6,nrow(chr_tb2)),1:nrow(chr_tb2),chrs,cex=0.8,adj=c(1,0.5))
  
##axis(side = 2,las = 1,lwd=1,at =1:nrow(chr_tb2),labels= chrs,cex.axis = 1)
  axis(side = 1,las = 1,lwd=1,at =seq(0,max_plot_len,10e6),  labels=paste(0:(as.integer(max_plot_len/10e6))*10,"",sep="")  ,cex.axis = 1)

##axis(side = 2,las = 1,lwd=1,at =c(-ymin,ymin),labels=c(signif(-ymin,2),signif(ymin,2)),cex.axis = 0.5)
######################################################################################
#################   plot methylation state

  region_tb1 <- read.delim(region_fname,sep="\t",header=T)
## region_tb1[,"chr"] <- paste("chr",region_tb1[,"chr"],sep="")

  for(i in 1:nrow(region_tb1)){
    
    if( region_tb1[i,"value"]=="red" ){
      rect(xleft=region_tb1[i,"start"],xright=region_tb1[i,"end"]
           ,ytop=as.numeric(match(region_tb1[i,1],chrs)+0.4)
           ,ybottom=as.numeric(match(region_tb1[i,1],chrs))
           ,col= "red"  ,border="red" )
    }
    if( region_tb1[i,"value"]=="green" ){
      rect(xleft=region_tb1[i,"start"],xright=region_tb1[i,"end"]
           ,ytop=as.numeric(match(region_tb1[i,1],chrs))
           ,ybottom=as.numeric(match(region_tb1[i,1],chrs)-0.4)
           ,col= "green"  ,border="green" )
    }
  }
  dev.off();
##file.remove(region_name)
}#End of function drawPics

files=list.files(pattern="*.plot.txt$")
for(k in 1:length(files)){
  drawPics(files[k]);
}


