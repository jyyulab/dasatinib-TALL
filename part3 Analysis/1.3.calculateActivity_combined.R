#!/usr/bin/env Rscript

##NetBID combined dataset activity
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
setwd("./Analysis/combined")
##### 1.load eset and network #####
load("./DATA/Discovery/Discovery_esets") #Discovery_FPKM.eset
load("./DATA/Validation/Validation_esets") #Validation_FPKM.eset
load("./DATA/TARGET/TARGET_esets") #TARGET_FPKM.eset
load("./DATA/transfer_tab.RData") #transfer_tab
load("./DATA/common_feature.RData") #common_f
load("./SJAR/TARGET_TALL_merge.network.RData")
network<-merge.network #network list
phe_info<-read.xlsx("./DATA/pheno_info.xlsx")
#####2.combine eset#####
eset<-NetBID2::merge_eset(eset1=Discovery_FPKM.eset,eset2 = Validation_FPKM.eset,group1 = "Discovery",group2 = "Validation",use_col = "sampleID")
eset<-NetBID2::merge_eset(eset1=eset,eset2 =TARGET_FPKM.eset,group1 = "Discovery_Validation",group2 = "TARGET",use_col = "sampleID")
eset<-NetBID2::update_eset.feature(use_eset = eset,use_feature_info = transfer_tab,from_feature = "gene",to_feature = "external_gene_name")
pd<-data.frame(sampleID=pData(eset)$sampleID)
pd<-left_join(pd,phe_info,by="sampleID")
rownames(pd)<-pd$sampleID
pData(eset)<-pd
exprs(eset)<-log2(exprs(eset)+0.1)
#####3.rmbatch effect#####
NetBID2::draw.eset.QC(eset=eset,outdir = "./QC/combined",do.logtransform = F,intgroup = "batch1",prefix = "3batch_")
NetBID2::draw.eset.QC(eset=eset,outdir = "./QC/combined",do.logtransform = F,intgroup = "batch2",prefix = "2batch_")
eset<-eset[common_f,]
NetBID2::draw.eset.QC(eset=eset,outdir = "./QC/combined",do.logtransform = F,intgroup = "batch1",prefix = "3batch_comfeature_")
NetBID2::draw.eset.QC(eset=eset,outdir = "./QC/combined",do.logtransform = F,intgroup = "batch2",prefix = "2batch_comfeature_")
exprs(eset)<-limma::removeBatchEffect(exprs(eset),batch = pd$batch2)
NetBID2::draw.eset.QC(eset=eset,outdir = "./QC/combined",do.logtransform = F,intgroup = "batch2",prefix = "rm2batch_comfeature_")
Combined_cal.eset<-eset

ac_mat <- NetBID2::cal.Activity(target_list=merge.network$target_list,cal_mat=exprs(eset),es.method='weightedmean')

Combined_ac.eset <- NetBID2::generate.eset(exp_mat=ac_mat,phenotype_info=pData(eset)[colnames(ac_mat),],feature_info=NULL,annotation_info='TARGETnet')
NetBID2::draw.eset.QC(eset=Combined_ac.eset,outdir="./QC/Combined",intgroup='batch2',do.logtransform=FALSE,prefix='AC_')

save(Combined_ac.eset,file = "./DATA/combined/Combined_ac.eset")
save(Combined_cal.eset,file = "./DATA/combined/Combined_cal.eset")