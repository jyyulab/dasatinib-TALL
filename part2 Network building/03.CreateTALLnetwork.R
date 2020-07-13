#!/usr/bin/env Rscript
##create TALL_network
##Coded by Jingjing Liu(Jingjing.liu@stjude.org)
##R version 3.6.3 (2020-02-29)

#####0.load packages#####
rm(list = ls())
require(dplyr)
require(openxlsx)
require(NetBID2)
setwd("./SJAR")
out.dir<-"./QC/network"
#####1.create TALL_network#####
tf_file<-"../../TF/SJARACNE_out.final/consensus_network_ncol_.txt"
sig_file<-"../../SIG/SJARACNE_out.final/consensus_network_ncol_.txt"

tf_network<-NetBID2::get.SJAracne.network(network_file = tf_file)
sig_network<-NetBID2::get.SJAracne.network(network_file = sig_file)
NetBID2::draw.network.QC(tf.network$igraph_obj,outdir=out.dir,prefix='TF_net_')
NetBID2::draw.network.QC(sig.network$igraph_obj,outdir=out.dir,prefix='SIG_net_')
merge.network<-NetBID2::merge_TF_SIG.network(TF_network = tf_network,SIG_network = sig_network)
save(merge.network,file ="TARGET_TALL_merge.network.RData")
