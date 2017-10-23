#!/usr/bin/bash 

source input.txt
 
if [ ! -d ${out}/pbs ];
then 
   mkdir -p ${out}/pbs
fi
if [ ! -d ${out}/bam ]
then 
mkdir -p ${out}/bam 
fi 
if [ ! -d ${out}/fq/qc ]
then 
mkdir -p ${out}/fq/qc
fi 


if [ ! -d ${out}/peaks ];
then
    mkdir -p ${out}/peaks
fi

if [ ! -d ${out}/peaks/annotation ];
then
    mkdir -p ${out}/peaks/annotation
fi

for i in $f_names
    do 
            cat > ${out}/pbs/"run_FqToSam_CHIP_"$i".pbs" <<stop
#PBS -N $project.$i
#PBS -j oe
#PBS -o ${out}
#PBS -e ${out}
#PBS -l nodes=1:ppn=4
#PBS -r y
#PBS -q blast
stop
     done


if [[ $seqtype = "single" ]];
then
    for i in $f_names
       do 
                  cat >> ${out}/pbs/"run_FqToSam_CHIP_"$i".pbs" <<stop
#### remove 1 adapter
gzip -dc  ${out}/fq/${i}.fq.gz | $filter_fq -Q 33 -q 20 -p 50 -v 2>${out}/fq/${i}.out  | \\
$trimer_fq -Q 33 -t 20 -l 20 -v 2>>${out}/fq/${i}.out | \\
$clipper_fq -Q 33 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 20 -M 10 -i - -o ${out}/fq/qc/clip_trim_q20_${i}.fq.gz -v 1>>${out}/fq/${i}.out

##### bowtie mapping
   ${bowtie2} -p 8     -I 0  -X 500     -x ${bowtie2_ref}      ${out}/fq/qc/clip_trim_q20_${i}.fq.gz -S ${out}/bam/${i}.sam   1>${out}/bam/${i}_bowtie1.out  2>${out}/bam/${i}_bowtie1.error
stop
     done
fi

if [[ $seqtype = "pair" ]];
then 
   for i in $f_names
       do
                 cat >> ${out}/pbs/"run_FqToSam_CHIP_"$i".pbs" <<stop
#### remove 1 adapter
gzip -dc  ${out}/fq/${i}_1.fq.gz | $filter_fq -Q 33 -q 20 -p 50 -v 2>${out}/fq/${i}.out  | \\
$trimer_fq -Q 33 -t 20 -l 20 -v 2>>${out}/fq/${i}.out | \\
$clipper_fq -Q 33 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 20 -M 10 -i - -o ${out}/fq/qc/clip_trim_q20_${i}_1.fq -v 1>>${out}/fq/${i}.out  && \\
gzip -dc  ${out}/fq/${i}_2.fq.gz | $filter_fq -Q 33 -q 20 -p 50 -v 2>>${out}/fq/${i}.out  | \\
$trimer_fq -Q 33 -t 20 -l 20 -v 2>>${out}/fq/${i}.out | \\
$clipper_fq -Q 33 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -l 20 -M 10 -i - -o ${out}/fq/qc/clip_trim_q20_${i}_2.fq -v 1>>${out}/fq/${i}.out  && \\

gzip ${out}/fq/qc/clip_trim_q20_${i}_1.fq
gzip ${out}/fq/qc/clip_trim_q20_${i}_2.fq

perl /home/huangf/hf_scripts/clean_no_pair.pl ${out}/fq/qc/clip_trim_q20_${i}_1.fq.gz  ${out}/fq/qc/clip_trim_q20_${i}_2.fq.gz ${out}/fq/qc/trim_paired_${i}_1.fq.gz ${out}/fq/qc/trim_paired_${i}_2.fq.gz ${out}/fq/qc/trim_single_${i}.fq.gz >>${out}/fq/${i}.out

##### bowtie mapping
   ${bowtie2} -p 8     -I 0  -X 500     -x ${bowtie2_ref}    -1 ${out}/fq/qc/trim_paired_${i}_1.fq.gz -2 ${out}/fq/qc/trim_paired_${i}_2.fq.gz -S ${out}/bam/${i}.sam   1>${out}/bam/${i}_bowtie1.out  2>${out}/bam/${i}_bowtie1.error
stop
     done
fi

#if [ -e  ${out}/bam/${i}_success ]; then 
    #echo "${i} bowtie mapping has been performed";
#else

for i in $f_names
    do
          cat >> ${out}/pbs/"run_FqToSam_CHIP_"$i".pbs" <<stop
##### sam to bam
   ${samtools} view -bS -o  ${out}/bam/${i}.bam  ${out}/bam/${i}.sam   && \\ 
   rm -f ${out}/bam/${i}.sam
###### bam sort
    ${samtools} sort ${out}/bam/${i}.bam   ${out}/bam/${i}_sorted      && \\ 
    rm -f ${out}/bam/${i}.bam
###### bam to sam
####  bam remove informal chr. 
     ${samtools}  view  -h -L ${reference_length}   ${out}/bam/${i}_sorted.bam | ${samtools} view -bS -   >${out}/bam/${i}.bam   
     ${samtools} index ${out}/bam/${i}.bam 
    
####bam to cov 

perl  /home/huangf/hf_scripts/genome_cov/split_genome_by_100k.pl ${genomeFai} | ${bedtools}  coverage -split  -counts -abam ${out}/bam/${i}.bam -b - >${out}/bam/${i}_genome.cov 
  perl -lane '\$t=log(\$F[3]+1)/log(10);print join("\t",@F[0..2],\$t)' ${out}/bam/${i}_genome.cov >${out}/bam/${i}_genome.log 

####macs call peaks

cd ${out}/peaks
/home/huangf/bin/python ${macs} -f BAM --pvalue ${macsPvalue} -g ${genomeSize} -t ${out}/bam/${i}.bam -n ${i}_macs 1>${i}_macs.error 2>${i}_macs.out  
$bedtools getfasta -fi $genomefa -bed ${i}_macs_peaks.bed -fo ${i}_macs.fasta
/home/huangf/bin/Rscript /home/huangf/hf_scripts/chrPeaks_hf.r $cytoband

stop
    done 

echo "#!/bin/bash" >${out}/pbs/run_FqToSam.CHIP.sh
for i in $f_names
do 
  echo "qsub run_FqToSam_CHIP_"$i".pbs" >> ${out}/pbs/run_FqToSam.CHIP.sh
done

chmod 777 ${out}/pbs/run_FqToSam.CHIP.sh
