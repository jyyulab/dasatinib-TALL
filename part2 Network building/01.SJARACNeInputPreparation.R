#!/usr/bin/env Rscript
##prepare TARGET TALL Network input 
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(dplyr)
require(openxlsx)
require(NetBID2)
require(arrayQualityMetrics)
setwd("./DATA/TARGET/")
sjar.dir<-("./SJAR")
#####1.hub gene list#####
sig_tf <- read.xlsx("SIG_TF_list.xlsx")
hub_sig<-sig_tf[grepl('SIG',sig_tf$funcType),"geneSymbol"]
hub_tf<-sig_tf[grepl('TF',sig_tf$funcType),"geneSymbol"]

#####2.QC and select samples#####
load("TARGET_TALL.RNASeq.rawcount.265.eset")
eset<-TARGET_TALL.RNASeq.rawcount.265.eset
exp_mat<-exprs(eset)
exp_mat<- NetBID2::RNASeqCount.normalize.scale(mat = as.matrix(exp_mat), total = 50000000) 
exp_mat <- log2(exp_mat+ 1)
eset<-NetBID2::generate.eset(exp_mat = exp_mat,feature_info = fData(TARGET_TALL.RNASeq.rawcount.265.eset),phenotype_info = pData(TARGET_TALL.RNASeq.rawcount.265.eset))

exp_mat<-exprs(eset)
choose1 <- NetBID2::IQR.filter(exp_mat = exp_mat,thre=0.5)
exp_mat<-exp_mat[choose1,]
qc.eset<-NetBID2::generate.eset(exp_mat = exp_mat,phenotype_info = pData(eset))

arrayQualityMetrics(expressionset = qc.eset,
                    outdir = "./QC/TARGET/",
                    force = TRUE,
                    intgroup = c('group'),reporttitle = "IQR0.5_group")

eset<-eset[apply(exprs(eset)<=2,1,sum)<=ncol(eset)*0.90,]
rm_samples<-c("SJALL016444_D1.TARGET.10.PATKYI.09A.01","SJTALL002078_D1.TARGET.10.PATMXN.09A.01D","SJTALL002050_D1","SJTALL022093_D2")
eset.sel<-eset[,!pData(eset)$sampleName%in%rm_samples]

arrayQualityMetrics(expressionset = eset.sel,
                    outdir = "./QC/TARGET/",
                    force = TRUE,
                    intgroup = c('group'),reporttitle = "rmbad_group")

#####3.select features#####
eset.sel<-eset.sel[apply(exprs(eset.sel)<=2,1,sum)<=ncol(eset.sel)*0.90,]
exp_mat <- exprs(eset.sel)
fd <- fData(eset.sel);fd$IQR <- apply(exp_mat, 1, IQR)
th<-quantile(as.numeric(exp_mat),0.1);th
eset.sel<-eset.sel[apply(exprs(eset.sel)>=th,1,sum)/ncol(eset.sel) >=0.1,]
IQR.cutoff <- 0.1; IQR.rescue <- 0 
fd.cutoff <- fd[fd$IQR >= quantile(fd$IQR, IQR.cutoff),]
fd.cutoff.sig <- intersect(fd.cutoff$geneSymbol, hub_sig)
fd.cutoff.tf <- intersect(fd.cutoff$geneSymbol, hub_tf)

fd.rescue <- fd[fd$IQR >= quantile(fd$IQR, IQR.rescue) &fd$geneSymbol%in%union(hub_sig,hub_tf),]
fd.rescue.sig <- intersect(fd.rescue$geneSymbol, hub_sig)
fd.rescue.tf <- intersect(fd.rescue$geneSymbol, hub_tf)
  
fd.sel.sig <- union(fd.cutoff.sig, fd.rescue.sig)
fd.sel.tf <- union(fd.cutoff.tf, fd.rescue.tf)
n.sig.total <- length(hub_sig); n.sig.sel <- length(fd.sel.sig); r.sig = n.sig.sel/n.sig.total; n.sig.cutoff <- length(fd.cutoff.sig); n.sig.rescue <- n.sig.sel - n.sig.cutoff
n.tf.total <- length(hub_tf); n.tf.sel <- length(fd.sel.tf); r.tf = n.tf.sel/n.tf.total; n.tf.cutoff <- length(fd.cutoff.tf); n.tf.rescue <- n.tf.sel - n.tf.cutoff
cat("Hub\t#Total\t#Filtered\t%Filtered\t#Cutoff\t#Rescue\nSIG\t", n.sig.total, "\t", n.sig.sel, "\t", r.sig, "\t", n.sig.cutoff, "\t", n.sig.rescue, "\nTF\t", n.tf.total, "\t", n.tf.sel, "\t", r.tf, "\t", n.tf.cutoff, "\t", n.tf.rescue, "\n")
  
#####4.write outputs#####
fd.filtered.index <- union(fd.cutoff$geneSymbol, fd.rescue$geneSymbol);
fd.filtered <- fd[fd$geneSymbol%in%fd.filtered.index,];
exp_mat.filtered <- exp_mat[row.names(exp_mat)%in%fd.filtered.index,]
exp_mat.out <- data.frame(geneSymbol = fd.filtered$geneSymbol, exp_mat.filtered)
SJAR.hub_genes.tf <-file.path(sjar.dir, 'tf.txt')
SJAR.hub_genes.sig <-file.path(sjar.dir, 'sig.txt')
SJAR.expression_matrix <-file.path(sjar.dir, 'input.exp')

write.table(
	exp_mat.out,
    file = SJAR.expression_matrix,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE)

cat(intersect(fd.filtered$geneSymbol, hub_tf),file = SJAR.hub_genes.tf,sep = '\n')
cat(intersect(fd.filtered$geneSymbol, hub_sig),file = SJAR.hub_genes.sig,sep = '\n')