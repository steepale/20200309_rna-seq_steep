#'---
#' title: "PASS1A (Rat) RNA-Seq: Outliers and Duplicate Samples"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200309"
#' output:
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'     
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Goals of Analysis
#' * Investigate samples for a correlation of 1 (duplicate samples) and annotate them
#' * Examine for any obvious outliers
#'     
#' ## Setup the Environment

#+ Setup Environment, message=FALSE, results='hide', warning = FALSE
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("org.Rn.eg.db")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2","GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr","org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather","PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","reshape2","xtable")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

############################################################
##### Functions ############################################
############################################################

# Set select
select <- dplyr::select

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and AUthor
################################################################################
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "steep"
################################################################################

## explicit gc, then execute `expr` `n` times w/o explicit gc, return timings
################################################################################
benchmark <- function(n = 1, expr, envir = parent.frame()) {
        expr <- substitute(expr)
        gc()
        map(seq_len(n), ~ system.time(eval(expr, envir), gcFirst = FALSE))
}
################################################################################

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

#' ## Load & Clean Data
#' ##### Data files to load:
#' * Count Matrix and Metadata Table from:
#'     * RNA-Seq from Mt. Sinai
#'         * 3 sequencing batches & metadata
#'     * RNA-Seq from Stanford
#'         * 2 sequencing batches & metadata

#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Files last saved in: 20200309_exploration-rna-seq-phase1_steep.R

# Count matrix
in_file <- paste0(WD,'/data/20200309_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv')
count_data <- read.table(in_file,sep = ',', header = TRUE,row.names = 1,check.names = FALSE)

# Meatdata table
in_file <- paste0(WD,'/data/20200309_rnaseq-meta-pass1a-stanford-sinai_steep.txt')
col_data <- read.table(in_file, header = TRUE, check.names = FALSE, sep = '\t')
row.names(col_data) <- col_data$sample_key

# Adjust columns
col_data$animal.key.exlt4 <- as.factor(col_data$animal.key.exlt4)



#+ Examination of Duplicate Samples
################################################################################
######### Duplicate Samples ####################################################
################################################################################

# Collect vial_labels that are duplicated
dup_vials <- col_data$vial_label[duplicated(col_data$vial_label)] %>% unique()

# Let's examine duplicates within sequencing batches

dups_seq <- vector()
for(dup in dup_vials){
        # Collects the number of sequencing batches per sample
        n_seq <- col_data %>% filter(vial_label %in% dup) %>%
                select(Seq_batch) %>% unique() %>% unlist() %>% length()
        if(n_seq > 1){
                dups_seq <- c(dups_seq, dup)
        }
}
# Collect dataframe of only duplicate samples
dup_data <- col_data %>%
        filter(vial_label %in% dup_vials)
dup_data$vial_label <- factor(dup_data$vial_label)

#' #### Duplicate samples were included across 2 batches: There are 40 pairs of duplicates
#+results='asis'
print(xtable(table(dup_data$vial_label, dup_data$Seq_batch)), type = 'html', include.rownames = T)
#table(dup_data$vial_label, dup_data$Seq_batch)

#+ Visualize Variance Between Duplicates & Sequencing Batch
################################################################################
######### PCA & Correlation Matrix #############################################
################################################################################


#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(status) == colnames(all.data.m))

#' ## PCA Visualization of Sequencing Batches (Unsupervised)
#' TODO: Generate a summary of major inferences
#' 

#+ PCA Sequencing Batches (Unsupervised)
##########################################################################
############ PCA Sequencing Batches (Unsupervised) #######################
##########################################################################

# Perform unsupervised clustering
# Build model
count_data <- all.data.m
col_data <- status # To look at all samples 
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet
dds

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
#counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
#counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
#counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
#counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")

#' #### Let's examine the stucture in the data to see if there is significant sample correlation
#' It looks like much of the sample correlation is driven by sequencing batch
s <- svd(norm_counts)
dim(s$v)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
image ( cor(norm_counts) ,col=cols,zlim=c(-1,1))

#' #### Variance-explained plot: explained variance for PCs
#' ##### This is what independent data would look like:
norm_counts0 <- matrix( rnorm( nrow(norm_counts)*ncol(norm_counts) ) , nrow(norm_counts), ncol(norm_counts) )
d0 <- svd(norm_counts0)$d
plot(d0^2/sum(d0^2),ylim=c(0,.25))

#' ##### This is what these data look like: Which shows that a majority of variance is explained by just 2 principle components
plot(s$d^2/sum(s$d^2))





# Outlier Detection
################################################################################
#DESeq2::plotPCA(vstd, intgroup ="animal.registration.sex") +
#        guides(color=guide_legend(title="Sex"))

#par(mar=c(8,5,2,2))
#boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#assays(dds)[["cooks"]]

#https://support.bioconductor.org/p/35918/
################################################################################



#' #### End of Script Notes

#+ End of Script Notes
################################################################################
################ End of Script Notes ###########################################
################################################################################


#' #### Session Information
################################################################################
####################### Session Info ###########################################
################################################################################
session_info()

# Stops evaluation of code
#knitr::opts_chunk$set(eval = FALSE)

