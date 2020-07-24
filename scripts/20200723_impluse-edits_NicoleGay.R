#!/bin/R

TISSUE = 'Kidney'
numCores = as.numeric(1)

BiocManager::install("darch")
install.packages("darch")
install.packages("darch")

library(data.table)
library(edgeR)
library(ggplot2)
library(clusterProfiler)
library(org.Rn.eg.db)
library(DESeq2)
library(darch)
library(pracma)
library(svMisc)

library(foreach)
library(doParallel)
registerDoParallel(numCores)


WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
source(paste0(WD,'/scripts/20200723_impulse-fx-edits_NicoleGay.R')) # load ancillary functions 

TISSUE_LAB = gsub(' ','_',TISSUE)
TIMEPOINTS = c('0h','0.5h','1h','4h','7h','24h','48h')

# setwd('/projects/motrpac/PASS_ANALYSIS/rna/impulse_de/')
# outdir = sprintf('rdata/optim/%s_%s',TISSUE_LAB,SEX)
# system(sprintf('mkdir -p %s',outdir))

# load DESeq log2FCs or TMM-normalized data along with metadata 
load(sprintf('/projects/motrpac/PASS_ANALYSIS/rna/STARTING_DATA/unadj_%s_%s_control_processed_data.RData', TISSUE_LAB, SEX))
# logfcres, padjres, de_dt, counts_round, meta_data, tmm 

padj = padjres # DESeq adjusted p-values for comparison between each tp and combined controls (data.table)
logfc = logfcres # DESeq log2FCs for comparison between each tp and combined controls (data.table)
colnames(padj) = gsub('.*_|h','',colnames(padj))
colnames(logfc) = gsub('.*_|h','',colnames(logfc))

# select genes however you want 
keep = c()
# some gene selection process to define the "keep" list... 
# ...
# ...
logfc_filt = logfcres[gene %in% keep]



# parallelize process by gene instead, which should always be maxing out a CPU
# saves one RData object per gene
impulse_optimization = foreach (i = 1:length(keep), .verbose=T, .combine = rbind) %dopar% {
        maximize_single_gene(keep[i], logfc_filt, outdir, .N=100, .maxiter=1000)
        #.N: number of different initializations
        #.maxiter: max number of iterations to optimize SSE
}

q()

###################################################################################################
# RUN THIS AFTER THE OPTIMIZATION IS DONE
# merge impulse results

setwd('/projects/motrpac/PASS_ANALYSIS/rna/impulse_de/')
indir = sprintf('rdata/optim/%s_%s',TISSUE_LAB,SEX)

# read in gene-level results
# merge dfs
dflist = list()
# load saved dfs to get a sense of the data 
for (df_file in list.files(path=indir,
                           pattern='RData',
                           full.names = T)){
        load(df_file)
        dflist[[df_file]] = df
}
impulse_optimization = data.table(rbindlist(dflist))