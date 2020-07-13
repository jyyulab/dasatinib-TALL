#!/usr/bin/env Rscript

##NetBID analysis of TARGET cohort
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
setwd("./Analysis/TARGET")

##### 1.load eset and network #####
load("./DATA/TARGET/TARGET_esets")
eset<-TARGET_FPKM.eset
load("./SJAR/TARGET_TALL_merge.network.RData")
network<-merge.network

#####2.select features #####
exp_mat<-exprs(eset)
choose1 <- apply(exp_mat<= 0, 1, sum)<= ncol(exp_mat) * 0.50
exp_mat<-exp_mat[choose1,]
exp_mat<-log2(exp_mat+0.1)
choose2<-IQR.filter(exp_mat = exp_mat,use_genes = rownames(exp_mat),thre = 0.4)
exp_mat<-exp_mat[choose2,]
eset<-generate.eset(exp_mat = exp_mat,phenotype_info = pData(eset)[colnames(exp_mat),],feature_info = fData(eset)[rownames(exp_mat),])
draw.eset.QC(eset = eset,outdir="./QC/TARGET",intgroup='condition',do.logtransform=FALSE,prefix='after_QC_')

#####3.convert ID #####
NetBID2::db.preload(use_level = "gene",use_spe = "human")
use_genes<-fData(eset)$gene
transfer_tab<-NetBID2::get_IDtransfer(from_type ="ensembl_gene_id_version",to_type = "external_gene_name",dataset = "hsapiens_gene_ensembl",use_genes = use_genes,ignore_version = TRUE)
fData(eset)$ensembl_gene_id<-gsub("\\..*$","",fData(eset)$gene)
transfer_tab<-left_join(fData(eset),transfer_tab,by="ensembl_gene_id")
eset<-NetBID2::update_eset.feature(use_eset = eset, use_feature_info = transfer_tab,from_feature = "gene",to_feature = "external_gene_name")
TARGET_cal.eset<-eset
#####4.calculate activity #####
ac_mat <- NetBID2::cal.Activity(target_list=merge.network$target_list,cal_mat=exprs(eset),es.method='weightedmean')
TARGET_ac.eset <- NetBID2::generate.eset(exp_mat=ac_mat,phenotype_info=pData(eset)[colnames(ac_mat),],feature_info=NULL,annotation_info='TARGETnet')
NetBID2::draw.eset.QC(eset=TARGET_ac.eset,outdir="./QC/TARGET",intgroup='condition',do.logtransform=FALSE,prefix='AC_')

save(TARGET_ac.eset,file = "./DATA/TARGET/TARGET_ac.eset")
save(TARGET_cal.eset,file = "./DATA/TARGET/TARGET_cal.eset")