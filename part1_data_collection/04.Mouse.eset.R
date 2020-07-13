#!/usr/bin/env Rscript
##prepare mouse developement gene expression eset
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(NetBID2)
require(dplyr)
require(arrayQualityMetrics)
setwd("./DATA/Mouse_develop")

#####1.load GEO#####
gset <- NetBID2::load.exp.GEO(out.dir='./',
                         GSE='GSE15907',
                         GPL='GPL6246',
                         getGPL=TRUE,
                         update=FALSE)



#####2.QC#####
exprs(gset)<-limma::normalizeQuantiles(exprs(gset))
exprs(gset)<-log2(exprs(gset))

pData(gset)$condition<-gsub("#.*","",pData(gset)$title)
table(pData(gset)$condition)

#2.1 ID convert
mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
a<-getBM(attributes=c('affy_mogene_1_0_st_v1', 'external_gene_name',"ensembl_gene_id"), 
         filters = 'affy_mogene_1_0_st_v1',
         values = as.character(unique(fData(gset)$ID)), 
         mart = mouse)
mm.gene<-as.character(unique(a$ensembl_gene_id))
m2h.gene <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = mm.gene , mart = mouse, attributesL = c("external_gene_name","ensembl_gene_id"), martL = human, uniqueRows=T)

names(m2h.gene)<-c("mm.ensembl_gene_id","geneSymbol","hg.ensembl_gene_id")
m2h.gene<-m2h.gene[!is.na(m2h.gene$hg.ensembl_gene_id),]

eset<-gset
fData(eset)$IQR<-apply(exprs(eset),1,IQR)
fd<-fData(eset)[,c("ID","IQR")]
fd<-left_join(fd,a,by=c("ID"="affy_mogene_1_0_st_v1"))
fd<-arrange(fd,desc(IQR))
r<-which(duplicated(fd$ID)|is.na(fd$ensembl_gene_id)|duplicated(fd$ensembl_gene_id))
fd<-fd[-r,]
eset<-NetBID2::generate.eset(exp_mat = exp,phenotype_info = pData(eset),feature_info = fd)

transfer_tab<-left_join(fd,genesV2,by=c("ensembl_gene_id"="mm.ensembl_gene_id"))
transfer_tab<-transfer_tab[!is.na(transfer_tab$hg.ensembl_gene_id)|!duplicated(transfer_tab$hg.ensembl_gene_id),]
r<-which(is.na(transfer_tab$hg.ensembl_gene_id)|duplicated(transfer_tab$hg.ensembl_gene_id))
transfer_tab<-transfer_tab[-r,]

eset<-NetBID2::update_eset.feature(use_eset = eset,from_feature = "ID",to_feature = "geneSymbol",use_feature_info = transfer_tab)

#2.2 qc
exp_mat<-exprs(eset)
choose1 <- IQR.filter(exp_mat = exp_mat,thre=0.5)
exp_mat<-exp_mat[choose1,]
qc.eset<-NetBID2::generate.eset(exp_mat = exp_mat,phenotype_info = pData(eset))

arrayQualityMetrics(expressionset = qc.eset,
                    outdir = out.dir,
                    force = TRUE,
                    intgroup = c('phenotype markers:ch1'),reporttitle = "IQR0.5_phenomarker")

arrayQualityMetrics(expressionset = qc.eset,
                    outdir = out.dir,
                    force = TRUE,
                    intgroup = c('source_name_ch1'),reporttitle = "IQR0.5_source")

rm_samples<-("GSM777042","GSM777067","GSM854306","GSM854307","GSM854308","GSM854309","GSM854310","GSM854311","GSM854312","GSM854313","GSM854314","GSM854326","GSM854338","GSM854339","GSM854340")
eset<-eset[,!rm_samples]



#####3.sample selection#####
pd<-read.xlsx("./pheno_info_modified.xlsx") #S2 table1 from Mingueneau, M. et al. The transcriptional landscape of alphabeta T cell differentiation. Nat Immunol 14, 619-632, doi:10.1038/ni.2590 (2013)
rownames(pd)<-pd$geo_accession
pd<-pd[eset$geo_accession,]
pData(eset)<-pd
sel_samples<-pd$geo_accession[!is.na(pd$type_annotation)]
mouse.log2hgSymbol.eset<-eset[,sel_samples]
save(mouse.log2hgSymbol.eset,file = "./mouse.log2hgSymbol.eset")