#!/usr/bin/env Rscript

##analysis master regulators in GSE15907
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
#####1.load data#####
load("./DATA/Mouse_develop/mouse.log2hgSymbol.eset")
load("./SJAR/TARGET_TALL_merge.network.RData") #merge.network
eset<-mouse.log2hgSymbol.eset
exp_mat <- exprs(eset)
choose1 <- apply(exp_mat<= quantile(exp_mat, probs = 0.05), 1, sum)<= ncol(exp_mat) * 0.90
exp_mat <- exp_mat[choose1,]
cal_eset <- generate.eset(exp_mat=exp_mat, phenotype_info=pData(eset)[colnames(exp_mat),],feature_info = fData(eset)[rownames(exp_mat),])

mouse_cal.eset<-cal_eset
ac_mat <- cal.Activity(target_list =merge.network$target_list,cal_mat=exprs(cal_eset),es.method='weightedmean')
mouse_ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(cal_eset)[colnames(ac_mat),],feature_info=NULL,annotation_info='activity in TARGETnet')

save(mouse_ac.eset,file = "./DATA/Mouse_develop/mouse_ac.eset")
save(mouse_cal.eset,file = "./DATA/Mouse_develop/mouse_cal.eset")