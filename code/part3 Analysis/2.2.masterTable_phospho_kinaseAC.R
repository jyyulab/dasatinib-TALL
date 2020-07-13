#!/usr/bin/env Rscript

##analysis infered kinase activity in phosphoproteomic data
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
#####1.load data#####
load("./DATA/phosphoproteomic/phospho_ac.eset")
load("./DATA/phosphoproteomic/human_kinase_substrate.gsc")
#####2.differential kinase activity#####
eset<-phospho_ac.eset
phe_info<-pData(eset)

comps<-rbind(c("Sensitive_no_drug","Sensitive_Dasatinib_10nM"),c("Resistant_no_drug","Resistant_Dasatinib_10nM"),c("Sensitive_no_drug","Resistant_no_drug"),c("Resistant_Dasatinib_10nM","Sensitive_Dasatinib_10nM"))
DA<-list()
for(i in 1:dim(comps)[1]){
  comp_name <- sprintf('%s.Vs.%s',comps[i,2],comps[i,1]) 
  G0  <- rownames(phe_info)[which(phe_info$condition==comps[i,1])] # 
  G1  <- rownames(phe_info)[which(phe_info$condition==comps[i,2])] # 
  DA_gene_limma <- getDE.limma.2G(eset=eset,G1=G1,G0=G0,G1_name=comps[i,1],G0_name=comps[i,2])
  DA[[comp_name]] <- DA_gene_limma
}
use_comp<-names(DA)
combine_info<- lapply(use_comp,function(DA,x,z_col='Z-statistics',display_col=c('logFC','P.Value'),type="phospho"){
  ID<-rownames(DA[[1]])
  avg_col<-colnames(DA[[x]])[grep('^Ave',colnames(DA[[x]]))]
  DA_info <- DA[[x]][,c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col)))]
  colnames(DA_info) <- paste0(colnames(DA_info),'.',x,paste0("_",type))
  colnames(DA_info)[1] <- paste0('Z.',x,paste0("_",type))
  out <- data.frame(DA_info[ID,])
  out
},DA=DA)
ms_tab<-do.call(cbind.data.frame, combine_info)
ms_tab$gene<-rownames(ms_tab)
ms_tab<-left_join(fData(eset),ms_tab,by="gene")
tmp<-unlist(lapply(human_kinase_substrate.gsc, length))
tmp<-data.frame(tmp)
names(tmp)<-"Size"
tmp$gene<-rownames(tmp)
ms_tab<-left_join(ms_tab,tmp,by="gene")
file_name<-"./DATA/phosphoproteomic/kinaseAC_ms-tab.xlsx"
out2excel(ms_tab,out.xlsx = file_name)