#!/usr/bin/Rscript 
##usage Rscript step7_sigsite_sum_total.R 
#setwd("&&project/Dif_methylation")
#a = list.files(pattern = "*sigsite.txt$") ## com_difbeta  （com_difbeta*.txt)
#a = list.files(pattern = "com_difbeta.*txt$")
a = list.files(pattern = "*sigsite.txt$")									#批量读入文件
for (i in c(1:length(a))){
  data = read.delim(a[i],header =T)
  sum_num = nrow(data)														#行数
  hypermethylated_num = nrow(data[data$Beta.Difference > 0,]) 				#甲基化数量为Beta.Difference > 0
  hypomethylated_num = sum_num - hypermethylated_num						#低甲基化数等于行数-甲基化数
  out = data.frame(a[i],sum_num,hypermethylated_num,hypomethylated_num)
  write.table(out,file = "compare.txt",sep = "\t",append = TRUE,col.names =T,row.names = F )#生成可相加文件
}



a = list.files(pattern = "com_difbeta.*txt$")
#a = list.files(pattern = "*sigsite.txt$")
for (i in c(1:length(a))){
  data = read.delim(a[i],header =T)
  sum_num = nrow(data)
  hypermethylated_num = nrow(data[data$Beta.Difference > 0,])
  hypomethylated_num = sum_num - hypermethylated_num
  out = data.frame(a[i],sum_num,hypermethylated_num,hypomethylated_num)
  write.table(out,file = "compare.txt",sep = "\t",append = TRUE,col.names =T,row.names = F )
}