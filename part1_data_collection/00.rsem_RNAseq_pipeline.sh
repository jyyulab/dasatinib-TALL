#!urs/bin/bash
#Run RSEM bowtie2 RNAseq pipeline
#Coded by Jingjing Liu (Jingjing.liu@stjude.org) 

#Reference genome: hg38 (GENCODE, V30)
#fastqc v0.11.5
#salmon v0.13.1
#picard v2.9.4
#java v1.8.0_60
#RSEM v1.2.28
#Bowtie2 v2.3.5 

rf=../../BOWTIE2_RSEM_index/hg38_bowtie2

#####0.convert bam to fastq#####
#use sample1 as example
java -java ../../path_to_picard/picard.jar SamToFastq I= ../../path_to_bams/sample1.bam F=../../path_to_fastq/sample1_R1.fastq.gz F2=../../path_to_fastq/sample1_R2.fastq.gz VALIDATION_STRINGENCY=SILENT

#####1.qc#####
for sample in sample1 sample2 sample3 ....
do
	fastqc ../../path_to_fastq/$sample_R1.fastq.gz -o ./$sample
done

#####2.quantification#####
for sample in sample1 sample2 sample3 ....
do
	rsem-calculate-expression --bowtie2 --bowtie2-path ~/bin/bowtie2 -p 8 --no-bam-output --paired-end  ../../path_to_fastq/$sample_R1.fastq.gz ../../path_to_fastq/$sample_R2.fastq.gz $rf ../../rsem_bowtie2_output/$sample_
done
