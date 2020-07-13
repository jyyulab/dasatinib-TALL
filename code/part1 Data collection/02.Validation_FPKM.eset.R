#!/usr/bin/env Rscript
##prepare Validation expression eset
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
setwd("./DATA/Validation")
source("./utility.R")

#####1.load matrix#####
res_dir<-("../../resem_result_path")
phe_info<-read.xlsx("./phe_info.xlsx")
phe_info<-phe_info[phe_info$original_group=="Validation",]

choose_format<-"FPKM"
choose_level<-"gene"
sample_names<-phe_info$sampleID
out_mat<-rsem2matrix(res_dir = res_dir,choose_format = choose_format,choose_level = choose_level,sample_names = sample_names)

#####2.generate eset#####
exp_mat<-out_mat[,-1]
exp_mat<-apply(exp_mat,2,as.numeric)
rownames(exp_mat)<-out_mat$gene
rownames(phe_info)<-phe_info$sampleID
eset<-NetBID2::generate.eset(exp_mat = exp_mat,phenotype_info=phe_info)
Validation_FPKM.eset<-eset

exprs(eset)<-log2(exprs(eset)+0.1)
Validation_log2FPKM.eset<-eset
#####3.QC#####
draw.eset.QC(Validation_FPKM.eset, outdir = "./QC", do.logtransform = FALSE, prefix = 'Validation_FPKM_',
             intgroup = "condition", choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.text' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)
draw.eset.QC(Validation_log2FPKM.eset, outdir = "./QC", do.logtransform = FALSE, prefix = 'Validation_log2FPKM_',
             intgroup = "condition", choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.text' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)
#save data
save(Validation_FPKM.eset,Validation_log2FPKM.eset,file = "./Validation_esets")
