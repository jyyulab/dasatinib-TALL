#!/usr/bin/env Rscript

##analysis master regulators in discovery cohort
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
#####1.load data#####
load("./DATA/Discovery/Discovery_ac.eset")
load("./DATA/transfer_tab.RData") 
load("./DATA/Discovery/Discovery_esets")
load("./SJAR/TARGET_TALL_merge.network.RData")
exp.eset<-Discovery_FPKM.eset;ac.eset<-Discovery_ac.eset
exp.eset<-NetBID2::update_eset.feature(use_eset = exp.eset,use_feature_info = transfer_tab,from_feature = "gene",to_feature = "external_gene_name") #convert ID
exprs(exp.eset)<-log2(exprs(exp.eset)+0.1)
#####2.differential expression and activity#####
DA<-list()
DE<-list()
phe_info<-pData(ac.eset)
comps<-c("sensitive","resistant")
comp_name <- sprintf('%s.Vs.%s',comps[1],comps[2])
)
G0  <- rownames(phe_info)[which(phe_info$condition==comps[2])] 
G1  <- rownames(phe_info)[which(phe_info$condition==comps[1])]

DE_gene_limma <- getDE.limma.2G(eset=exp.eset,G1=G1,G0=G0,G1_name=comps[1],G0_name=comps[2])
DA_driver_limma <- getDE.limma.2G(eset=ac.eset,G1=G1,G0=G0,G1_name=comps[1],G0_name=comps[2])
DE[[comp_name]]<-DE_gene_limma
DA[[comp_name]]<-DA_driver_limma
#####3.master table#####
db.preload(use_level='gene',use_spe='human',update=FALSE)
all_comp <- names(analysis.par$DE) 
ms_tab<-generate.masterTable(use_comp=all_comp,DE=DE,DA=DA,
                     target_list =merge.network$target_list,
                     tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                     main_id_type='external_gene_name')
write.xlsx(ms_tab,file = "./DATA/Discovery/Discovery_msTable.xlsx")