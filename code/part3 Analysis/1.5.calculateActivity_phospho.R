#!/usr/bin/env Rscript

##infer Kinase activity by phosphoproteomic data
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
#####1.load data#####
load("./DATA/phosphoproteomic/phospho.eset")
load("./DATA/phosphoproteomic/human_Kinase_substrateNet")
#####2.clean dataset
eset<-phospho.eset
fd<-fData(eset)
phe_info<-pData(eset)
phe_info$condition<-paste(phe_info$Dasatinib_sensitivity,phe_info$Treatment,sep = "_")
comps<-rbind(c("Sensitive_no_drug","Sensitive_Dasatinib_10nM"),c("Resistant_no_drug","Resistant_Dasatinib_10nM"),c("Sensitive_no_drug","Resistant_no_drug"),c("Resistant_Dasatinib_10nM","Sensitive_Dasatinib_10nM"))
DE<-list()
for(i in 1:dim(comps)[1]){
  comp_name <- sprintf('%s.Vs.%s',comps[i,2],comps[i,1]) 
  G0  <- rownames(phe_info)[which(phe_info$condition==comps[i,1])] # 
  G1  <- rownames(phe_info)[which(phe_info$condition==comps[i,2])] # 
  DE_gene_limma <- getDE.limma.2G(eset=eset,G1=G1,G0=G0,G1_name=comps[i,1],G0_name=comps[i,2])
  DE[[comp_name]] <- DE_gene_limma
}
use_comp<-names(DE)
combine_info<- lapply(use_comp,function(DA,x,z_col='Z-statistics',display_col=c('logFC','P.Value'),type="D-phospho"){
  ID<-rownames(DA[[1]])
  avg_col<-colnames(DA[[x]])[grep('^Ave',colnames(DA[[x]]))]
  DA_info <- DA[[x]][,c(z_col,avg_col,setdiff(display_col,c(z_col,avg_col)))]
  colnames(DA_info) <- paste0(colnames(DA_info),'.',x,paste0("_",type))
  colnames(DA_info)[1] <- paste0('Z.',x,paste0("_",type))
  out <- data.frame(DA_info[ID,])
  out
},DA=DE)

ms_tab<-do.call(cbind.data.frame, combine_info)
ms_tab$gene<-rownames(ms_tab)
ms_tab<-left_join(fData(eset),ms_tab,by=c("Peptides"="gene"))
df<-ms_tab[,c("Peptides","GN","Mod.sites","P.Value.Sensitive_Dasatinib_10nM.Vs.Sensitive_no_drug_D.phospho")]
names(df)[4]<-"Pval"
splitSite<-function(x,df){
  tmp<-df[df$Peptides==x,]
  if(grepl(",",tmp$Mod.sites)){
    sites<-unlist(strsplit(tmp$Mod.sites,","))
    n<-length(sites)
    out<-data.frame(Peptides=rep(tmp$Peptides,n),GN=rep(tmp$GN,n),Mod.sites=sites,Pval=rep(tmp$Pval,n))
  }else{out<-tmp}
  return(out)
}
tmp<-lapply(as.character(df$Peptides), splitSite,df=df)
tmp<-do.call("rbind",tmp)
tmp$GENE_SITE<-paste(tmp$GN,tmp$Mod.sites,sep = "_")
tmp<-arrange(tmp, tmp$Pval)
tmp<-tmp[!duplicated(tmp$GENE_SITE),]

#####3.create kinase-substrate list
net<-human_Kinase_substrateNet
net$KINASE<-toupper(net$KINASE)
net$SUB_GENE_SITE<-toupper(net$SUB_GENE_SITE)
net.sel<-net[net$SUB_GENE_SITE%in%tmp$GENE_SITE,]
net2list<-function(net,k_col,sub_site_col){
  df<-net[,c(k_col,sub_site_col)]
  names(df)<-c("kinase","sub_site")
  out_list<-list()
  k_names<-unique(df$kinase)
  for(k_name in k_names){
    tmp<-df[df$kinase==k_name,"sub_site"]
    tmp<-tmp[!is.na(tmp)]
    out_list[[k_name]]<-unique(as.character(tmp))
  }
  return(out_list)
}
gsc<-net2list(net = net.sel,k_col="GENE",sub_site_col = "SUB_GENE_SITE")

#########4.infer kinase activity
eset<-update_eset.feature(use_eset = eset,use_feature_info = tmp,from_feature = "Peptides",to_feature = "GENE_SITE")
ac_mat<-cal.Activity.GS(use_gs2gene = gsc,cal_mat = exprs(eset),es.method = "mean")
ac_eset<-generate.eset(exp_mat = ac_mat,phenotype_info = phe_info[colnames(ac_mat),])
phospho_cal.eset<-eset
phospho_ac.eset<-ac_eset
save(phospho_cal.eset,file = "./DATA/phosphoproteomic/phospho_cal.eset")
save(phospho_ac.eset,file = "./DATA/phosphoproteomic/phospho_ac.eset")