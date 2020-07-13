#!/usr/bin/env Rscript

##calculate biomarker score
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
#####1.load data#####
load("./DATA/combined/Combined_ac.eset")
load("./DATA/combined/preTCR_dasatinib_Targets")
ms_tab<-("./DATA/Discovery/Discovery_msTable.xlsx")
#####2.select drivers#####
ms_sel<-ms_tab[ms_tab$Size>50&ms_tab$Size<500&ms_tab$P.Value.sensitive.Vs.resistant_DA<1e-5,]
ms_sel<-arrange(ms_sel,ms_sel$P.Value.sensitive.Vs.resistant_DA)
ms_sel<-ms_sel[!duplicated(ms_sel$geneSymbol),]
targets<-unique(unlist(preTCR_dasatinib_Targets))
ms_sel<-ms_sel[ms_sel$geneSymbol%in%targets,]
ms_sel<-arrange(ms_sel,desc(ms_sel$Z.sensitive.Vs.resistant_DA))
tmp<-ms_sel[,c("originalID_label","Z.sensitive.Vs.resistant_DA")]
tmp$target<-tmp$originalID_label
tmp$MI<-abs(tmp$Z.sensitive.Vs.resistant_DA)
tmp$spearman<-sign(tmp$Z.sensitive.Vs.resistant_DA)
tmp<-tmp[,c("target","MI","spearman")]
rownames(tmp)<-tmp$target
marker_list<-list()
marker_list[["30biomarker"]]<-tmp
#####3.calculate biomarker score
eset<-Combined_ac.eset
ac_mat<-cal.Activity(target_list = marker_list,cal_mat = exprs(eset),es.method = "weightedmean")
ac_mat<-data.frame(t(ac_mat),check.names = F)
ac_mat$sampleID<-rownames(ac_mat)
out_df<-left_join(ac_mat,pData(eset),by="sampleID")
write.xlsx(out_df,file = "./DATA/combined/biomarker_score.xlsx")
