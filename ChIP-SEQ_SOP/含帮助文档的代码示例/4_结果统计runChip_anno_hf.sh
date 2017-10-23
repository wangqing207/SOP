#!/bin/bash

source Rscript/Rvariable.txt
source Rscript/Rvariable_anno.txt

if [ ! -d ${output}/pbs ];
then
  mkdir -p ${output}/pbs
fi

if [ ! -d ${output}/annotate ];
then 
   mkdir -p ${output}/annotate
fi 

if [ ! -d ${output}/fasta ];
then
   mkdir -p ${output}/fasta
fi

source ../input.txt

#  echo "$path2R medip-dmr.r $treatment $control $BSgenome_species $species $output $outfile_prefix $extend $ws $shift $uniq $FC"
##########
##########DIFFBIND_ANNOTATION###########
##########

samples=`ls ../peaks/*peaks.bed |cut -d "/" -f 3|cut -d "_" -f 1`
for i in $samples
do
/home/huangf/programs/BEDTools-Version-2.16.1/bin/bedtools getfasta -fi $genomefa -bed  ../peaks/${i}_macs_peaks.bed -fo fasta/${i}_macs_peaks.fasta
done


cat >$output/pbs/"run_ChIP_anno.pbs" <<stop
#PBS -N CHIPS_anno
#PBS -j oe
#PBS -o ${output}
#PBS -e ${output}
#PBS -l nodes=1:ppn=4
#PBS -r y
#PBS -q blast
stop

cat >>$output/pbs/"run_ChIP_anno.pbs" <<stop

$path2R /home/huangf/hf_scripts/Medips/Rscript/ann_region_useBed.R  $annotated_file_tss $gene_type_tss $file_to_ann_diffBind $output_dir $annotate_dir
$path2R /home/huangf/hf_scripts/Medips/Rscript/ann_region_useBed.R  $annotated_file_genetic $gene_type_genetic $file_to_ann_diffBind $output_dir $annotate_dir
$path2R /home/huangf/hf_scripts/chip-seq/Rscript/TSS_genetic_info_go_zjl.R $output_dir $refgene_fname $tax_id genetic $NCBI_database
$path2R /home/huangf/hf_scripts/Medips/Rscript/ann_region_useBed.R  $annotated_file_tss $gene_type_tss $file_to_ann_peaks $peaks_dir $peaks_dir
$path2R /home/huangf/hf_scripts/Medips/Rscript/ann_region_useBed.R  $annotated_file_genetic $gene_type_genetic $file_to_ann_peaks $peaks_dir $peaks_dir
$path2R /home/huangf/hf_scripts/chip-seq/Rscript/peaks_tss_genetic_zjl.r $peaks_dir $refgene_fname $tax_id genetic $NCBI_database
$path2R /home/huangf/hf_scripts/Medips/Rscript/peak_chromosome_length_distribution.r $BSgenome_species $species $peaks_dir
$path2R /home/huangf/hf_scripts/Medips/Rscript/genetic_distribution.r  $peaks_dir


stop


