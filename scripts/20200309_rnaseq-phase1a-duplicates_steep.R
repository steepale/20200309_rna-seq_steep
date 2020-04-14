#'---
#' title: "PASS1A (Rat) RNA-Seq: Duplicates"
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

#' ## Goal of Analysis:
#' Eighty pairs of lungs were sequenced across 2 batches (1 member of each pair for each batch); determine how well correlated these duplicate samples are and if they are appropriate to act as pseudo-technical replicates
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
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2","GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr","org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather","PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","reshape2","xtable","Hmisc")
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

# flattenCorrMatrix
################################################################################
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
        ut <- upper.tri(cormat)
        data.frame(
                row = rownames(cormat)[row(cormat)[ut]],
                column = rownames(cormat)[col(cormat)[ut]],
                cor  =(cormat)[ut],
                p = pmat[ut]
        )
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

#' ## Duplicate Samples

#+ Examination of Duplicate Samples
################################################################################
######### Duplicate Samples ####################################################
################################################################################

# Collect vial_labels that are duplicated
dup_vials <- col_data$vial_label[duplicated(col_data$vial_label)] %>% unique()
dup_samples <- col_data %>%
        filter(vial_label %in% dup_vials) %>%
        select(sample_key) %>% unlist() %>% as.character() %>% unique()

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
row.names(dup_data) <- dup_data$sample_key

#' #### Duplicate samples were included across 2 batches: There are 80 duplicate pairs (160 samples)
#+results='asis'
print(xtable(table(dup_data$vial_label, dup_data$Seq_batch)), type = 'html', include.rownames = T)
#table(dup_data$vial_label, dup_data$Seq_batch)

#' ## Subset, Normalize, and Transform Duplicates 

#+ Subset (Samples and Features), Normalize, and Transform Duplicates
################################################################################
############### Subset, Normalize, and Transform Duplicates  ###################
################################################################################

# Subset raw counts to only include duplicate values
dup_samples <- dup_data$sample_key %>% as.character() %>% unique() 
dup_counts <- count_data[, dup_samples]

#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
dup_counts <- dup_counts[, rownames(dup_data)]
all(rownames(dup_data) == colnames(dup_counts))

# Perform unsupervised clustering
# Build model
count_data <- dup_counts
col_data <- dup_data # To look at all duplicate samples
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 5 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowSums(counts(dds))/ncol(dds) >= 5
dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet
dds

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
for( n in 1){
        start_time <- Sys.time()
        #rld <- rlog(dds)
        rld <- vst(dds)
        end_time <- Sys.time()
        #print(end_time - start_time)
}

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")

#' ### Duplicate Sample Correlation: Heatmaps, Histograms, and Tables

#+ Duplicate Sample Correlation: Heatmaps, Histograms, and Tables
################################################################################
##### Duplicate Sample Correlation: Heatmaps, Histograms, and Tables ###########
################################################################################

#' #### Generate a heatmap of major annotations of interest
#' ##### Annotations of interest:
#' * Samples: sample_key
#' * Sequencing Batch: Seq_batch
#' * Sex: animal.registration.sex
#' * Exercise/Control Group: animal.key.anirandgroup
#' * Time of Death (binned): specimen.collection.t_death_bins.type
#' * Tissue: Tissue

# Select annotations of interest
ann_df <- col_data %>%
        select(sample_key,
               Tissue,
               Seq_batch,
               animal.registration.sex,
               animal.key.anirandgroup,
               specimen.collection.t_death_bins.type)
annotations <- c("sample_key",
                 "Tissue",
                 "Seq_batch",
                 "animal.registration.sex",
                 "animal.key.anirandgroup",
                 "specimen.collection.t_death_bins.type")
for(anno in annotations){
        ann_df[[anno]] <- as.character(ann_df[[anno]])
}
# Adjust names
names(ann_df) <- c('Sample','Tissue','Sequencing Batch',
                   'Sex','Exercise/Control Group','Time of Death (binned)')
ann_df <- ann_df %>%
        select('Sample','Exercise/Control Group','Time of Death (binned)',
               'Sex','Sequencing Batch','Tissue')
# Melt dataframe
ann_df_melt <- melt(ann_df, id.var = 'Sample')
# Adjust colors
n <- length(unique(ann_df_melt$value))
colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(n)

#' #### The distribution of Duplicate Samples by Annotation
# Plot heatmap
ggplot(ann_df_melt, 
       aes(variable, Sample)) + 
        geom_tile(aes(fill = value),
                  colour = "white") +
        scale_fill_manual(values=colors) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        theme(legend.position="none") +
        theme(axis.title.x=element_blank()) +
        theme(text = element_text(size=20))

#' #### All samples are lung, from either Mt Sinia 2 or 3, and sampled from different groups at different times of the day.
#+results='asis'
print(xtable(table(ann_df$`Time of Death (binned)`, ann_df$`Exercise/Control Group`)), type = 'html', include.rownames = T)

#' #### Let's examine the stucture in the data to see if there is significant sample correlation
#' Correlation matrix not shown, correlations too close to depict by color

# Collect the distribution of duplicate sample correlation values
cor_dups <- vector() %>% as.numeric()
for(dv in dup_vials){
        # Collect the sample keys
        ex_dups <- dup_data %>%
                filter(vial_label %in% dv) %>%
                select(sample_key) %>% unlist() %>% as.character()
        #collect the correlations
        cor_dups <- c(cor_dups, (cor(norm_counts[, ex_dups]))[1,2]) %>% as.numeric()
}
# Collect the distribution of non-duplicate sample correlation values
cor_nondups <- cor(norm_counts)[cor(norm_counts) %!in% cor_dups]
cor_nondups <- cor_nondups[cor_nondups != 1]

#' #### A histogram of correlation values demonstrates a bimodal distribution between duplicate and non-duplicate samples
mypar(1,2)
ggplot() + 
        geom_histogram(aes(x=cor_dups), 
                       position="identity", alpha=0.5,color="blue", fill="blue") +
        geom_histogram(aes(x=cor_nondups), 
                       position="identity", alpha=0.5,color="red", fill="red") +
        xlim(0.97,1.00) +
        labs(title="Correlation Values: \nDuplicate Samples across Batches (Blue) \nNon-duplicate Samples across and within Batches (Red)",
             x="Correlation Values", y = "Frequency")
#' #### Zoom in on the right distribution (no overlap)
ggplot() + 
        geom_histogram(aes(x=cor_dups), 
                       position="identity", alpha=0.5,color="blue", fill="blue") +
        geom_histogram(aes(x=cor_nondups), 
                       position="identity", alpha=0.5,color="red", fill="red") +
        xlim(0.995,0.999) +
        labs(title="Correlation Values--Bimodal: \nDuplicate Samples across Batches (Blue) \nNon-duplicate Samples across and within Batches (Red)",
             x="Correlation Values", y = "Frequency")

#' #### A boxplot and summary stats help explain the distribution
# Boxplot
ggplot() +
        geom_boxplot(aes(y=cor_dups, x="Duplicates"), 
                       color="blue", show.legend = FALSE) +
        geom_boxplot(aes(y=cor_nondups, x="Non-Duplicates"), 
                     color="red", show.legend = FALSE) +
        xlab("") +
        ylim(0.99,1)

#' #### Duplicate Sample Correlation Values:
summary(cor_dups)
#' #### Non-Duplicate Sample Correlation Values:
summary(cor_nondups)

#' ### Naive PCA of Duplicates

#+ Naive PCA of Duplicates
################################################################################
######################## Naive PCA of Duplicates ###############################
################################################################################

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, intgroup=c("Seq_batch"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Duplicate Lung Samples (80 Pairs)") +
        guides(color=guide_legend(title="Sequencing Batch"))

#' ## A Search for Hidden Duplicates
#' #### We now can see that duplicate samples sequenced across batches (atleast MSSM 2 & 3) show correlation values between 0.9973 - 0.9982 and that correlations outside of duplicatively sequenced samples do not demonstrate values within this range. When we search the entire dataset for samples with high correaltion values, we don't find any sample pair demonstrating correlation values as high as lungs from MSSM 2 & 3.

#+ A Search for Hidden Duplicates
################################################################################
#####     A Search for Hidden Duplicates      ##################################
################################################################################
# Files last saved in: 20200309_exploration-rna-seq-phase1_steep.R

# Count matrix
in_file <- paste0(WD,'/data/20200309_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv')
count_data <- read.table(in_file,sep = ',', header = TRUE,row.names = 1,check.names = FALSE)

# Meatdata table
in_file <- paste0(WD,'/data/20200309_rnaseq-meta-pass1a-stanford-sinai_steep.txt')
col_data <- read.table(in_file, header = TRUE, check.names = FALSE, sep = '\t')
row.names(col_data) <- col_data$sample_key

# Perform unsupervised clustering
# Build model
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 5 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowSums(counts(dds))/ncol(dds) >= 5
dds <- dds[keep,]
# #### Summary of counts and annotation data in a DESeqDataSet
#dds

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
for( n in 1){
        start_time <- Sys.time()
        #rld <- rlog(dds)
        rld <- vst(dds)
        end_time <- Sys.time()
        #print(end_time - start_time)
}

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)

#' ## Conclusion:
#' #### When we search for any sample pair that exhibits correlation > 0.997, we don't find any samples other than the duplicate lung samples. This exercise adds evidence that the only duplicated samples were lungs from MSSM 2 & 3.
#' 
# Examine the correlation values
res2 <- rcorr(norm_counts)
flat_df <- flattenCorrMatrix(res2$r, res2$P)
flat_df$row <- as.character(flat_df$row)
flat_df$column <- as.character(flat_df$column)
flat_df %>%
        filter(row != column) %>%
        filter(cor >= 0.997) %>% 
        filter(row %!in% dup_samples | column %!in% dup_samples)

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

