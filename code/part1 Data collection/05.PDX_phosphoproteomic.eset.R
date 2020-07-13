#!/usr/bin/env Rscript
##prepare phosphoproteomic eset
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#create phosphoproteomic eset

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
setwd("./DATA/phosphoproteomic")

#####1.load data#####
data.dir<-"../../path_to_phosphoproteomic_rawdata"
d<-read.xlsx(sprintf("%s/deep_phospho_v01.3.xlsx",data.dir),sheet = 1)
pd<-read.xlsx(sprintf("%s/pheno_info_phosphoproteomic.xlsx",data.dir),sheet=1)

#####2.create eset#####
p_mat<-d[,9:19]
rownames(p_mat)<-d$Peptides

fd<-d[,1:8]
rownames(fd)<-fd$Peptides

pd<-pd[match(colnames(p_mat),pd$sampleID),]
rownames(pd)<-pd$sampleID
phospho.eset<-generate.eset(exp_mat = p_mat,phenotype_info = pd,feature_info = fd,annotation_info = "phosphoproteomic")
exprs(phospho.eset)<-log2(exprs(phospho.eset)) 

#####3.QC#####
draw.eset.QC(phospho.eset, outdir = "./QC", do.logtransform = FALSE, prefix = 'phospho_log2intensity_',
             intgroup = NULL, choose_plot = c("heatmap", "pca", "density", "correlation", "meansd"), generate_html = TRUE, correlation_strategy = "pearson", plot_all_point = FALSE,
             emb_plot_type='2D.text' # "2D", "2D.interactive", "2D.ellipse", "2D.text" or "3D" 
)

#save data
save(phospho.eset,file = "./phospho.eset")