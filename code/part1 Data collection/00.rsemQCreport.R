#!/usr/bin/env Rscript
##prepare mouse developement gene expression eset
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#RSEM QC report
#####0.load packages#####
rm(list = ls())
require(dplyr)
require(openxlsx)
setwd("./QC/RSEM_report/")
source("./utility.R")

#####1.Discovery#####
res_dir<-("../../resem_result_path")
phe_info<-read.xlsx("./phe_info.xlsx")
phe_info<-phe_info[phe_info$original_group=="Discovery",]
sample_names<-as.character(phe_info$sampleID)
out_qc<-rsemQCreport(res_dir = res_dir,sample_names = sample_names)
write.xlsx(out_qc,file="Discovery_rsemQC.xlsx")

#####2.Validation#####
res_dir<-("../../resem_result_path")
phe_info<-read.xlsx("./phe_info.xlsx")
phe_info<-phe_info[phe_info$original_group=="Validation",]
sample_names<-as.character(phe_info$sampleID)
out_qc<-rsemQCreport(res_dir = res_dir,sample_names = sample_names)
write.xlsx(out_qc,file="Validation_rsemQC.xlsx")

#####3.TARGET#####
res_dir<-("../../resem_result_path")
phe_info<-read.xlsx("./phe_info.xlsx")
phe_info<-phe_info[phe_info$original_group=="TARGET",]
sample_names<-as.character(phe_info$sampleID)
out_qc<-rsemQCreport(res_dir = res_dir,sample_names = sample_names)
write.xlsx(out_qc,file="TARGET_rsemQC.xlsx")
