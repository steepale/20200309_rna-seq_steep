'---
#' title: 'PASS1A Rat Tissue: -- RNASeq QC Steps -- EDA'
#' author: 'Alec Steep & Jiayu Zhang' 
#' date: '20200707'
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

#' ## Goals of Analysis:
#' * 
#'
#' 
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
#source('https://bioconductor.org/biocLite.R')
#BiocManager::install('ImpulseDE2')
#install.packages('reldist')

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr','ImpulseDE2','reldist')
lapply(pacs...man, FUN = function(X) {
        do.call('library', list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
counts <- DESeq2::counts
map <- purrr::map

# Global options
options(dplyr.print_max = 100)

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
# Declare Tissue

# TISSUE: Hypothalamus, Liver, Kidney, Aorta, Adrenal, Brown Adipose, Cortex, Gastrocnemius, Heart, Hippocampus,Lung,Ovaries,PaxGene,Spleen,Testes, White Adipose

# SCN: Hypothalamus (Suprachiasmatic nucleus)
# LIV: Liver
# KID: Kidney
# AOR: Aorta
# SKM: Gastrocnemius
# HAT: Heart
# ADG: Adrenal gland
# BAT: Brown adipose tissue
# WAT: White adipose tissue
# COR: Cortex
# HIP: Hippocampus
# LUNG: Lung
# OVR: Ovaries
# SPL: Spleen
# TES: Testes

# Load the decision table
table_file <- paste0(WD,'/data/20200603_rnaseq-tissue-data-assambly-table_steep.txt')
df_tbl <- read.table(file = table_file,sep = '\t', header = T, check.names = F)

TISSUE <- "Kidney"
# for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
print(TISSUE)
TISSUE1 <- TISSUE
#TISSUE <- 'Lung'
# # Collect the formula
# design <- df_tbl %>%
# filter(Tissue == TISSUE) %>%
# select(Formula) %>% unique() %>% 
# unlist() %>% as.character() %>% as.formula()
# Collect the Outliers
OUTLIERS <- df_tbl %>%
        filter(Tissue == TISSUE) %>%
        select(Outliers) %>% unique() %>% 
        unlist() %>% as.character()
# Collect Adjusted_Variance
ADJ_VAR <- df_tbl %>%
        filter(Tissue == TISSUE) %>%
        select(Adjusted_Variance) %>% unique() %>% 
        unlist() %>% as.character()
# Collect the TIS Symbol
TIS <- df_tbl %>%
        filter(Tissue == TISSUE) %>%
        select(Tis) %>% unique() %>% 
        unlist() %>% as.character()

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

# Set a vector for Exercise/Control Levels and Colors
ec_levels <- c('Exercise - IPE',
               'Exercise - 0.5 hr',
               'Exercise - 1 hr',
               'Exercise - 4 hr',
               'Exercise - 7 hr',
               'Exercise - 24 hr',
               'Exercise - 48 hr',
               'Control - IPE',
               'Control - 7 hr')
ec_colors <- c('gold',
               'darkgoldenrod1',
               'orange',
               'darkorange',
               'darkorange2',
               'darkorange3',
               'darkorange4',
               'steelblue1',
               'steelblue4')

# Load Metadata and count data as R objects
################################################################################
# Restore the metadata object
meta_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
col_data <- readRDS(file = meta_file)
# Restore the count object
count_file <- paste0(WD, '/data/20200603_rnaseq-counts-pass1a-stanford-sinai-processed_steep.rds')
count_data <- readRDS(file = count_file)

#' ## Collect Samples of Interest and Normalize

#+ Collect Samples of Interest and Normalize
################################################################################
#####     Collect Samples of Interest and Normalize      #######################
################################################################################
gene_plot <- data.frame()
sample_plot <- data.frame()
for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
print(TISSUE)
TISSUE1 <- TISSUE
for(sex in c('Male','Female')){
# TODO Code for lung -- use only one batch
if(TISSUE == c('Gastrocnemius')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex))
}
# Adjust for sex for reproductive tissues
if(TISSUE == 'Testes'){
        sex = 'Male'
}
if(TISSUE == 'Ovaries'){
        sex = 'Female'
}
rownames(tod_cols) <- tod_cols$sample_key
# Collect samples without NA values in TOD
nona_sams <- tod_cols %>%
        select(sample_key) %>% unlist() %>% as.character()
# Collect tissue specific counts
tod_counts <- count_data[,nona_sams]
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))
# For each gene, collect the proportion of non-zero counts across samples within a tissue-sex cohort
# Samples in cohort
cohort_samples <- tod_cols %>%
        filter(Tissue == TISSUE) %>%
        filter(animal.registration.sex == sex) %>%
        select(sample_key) %>% unlist() %>% as.character()
# the subset counts
gene_counts <- count_data[,cohort_samples]

# Non-zero counts (per gene)
######################################################
# The fraction of non-zero counts per gene (function)
nonzero <- function(x) sum(x != 0)
# Non zero counts (per gene)
nzc <- apply(gene_counts,1,nonzero)
# Non zero fraction 
nzf <- nzc/length(cohort_samples)
# Gini index
ginic <- apply(gene_counts,1,reldist::gini)
# Mean gene expresison
muc <- rowMeans(gene_counts)
# SD gene expression
sd <- apply(gene_counts,1,sd)
# Coefficient of variation
cv <- sd/muc
# Save the dataframe
gene_df <- data.frame('ENSEMBL_RAT' = row.names(gene_counts),
                    'NZF' = nzf, 
                    'GINI' = ginic,
                    'CV' = cv,
                    'E_MU' = muc,
                    'SEX' = sex,
                    'TISSUE' = TISSUE1)

# Concatenate the dataframes
gene_plot <- rbind(gene_plot, gene_df)

# Non-zero counts (per sample)
######################################################
# The fraction of non-zero counts (function)
nonzero <- function(x) sum(x != 0)
# Non zero counts
nzc <- apply(gene_counts,2,nonzero)
ginic <- apply(gene_counts,2,reldist::gini)
# Non zero fraction 
nzf <- nzc/nrow(gene_counts)
# Mean gene expresison
m <- apply(gene_counts,2,mean)
# SD gene expression
sd <- apply(gene_counts,2,sd)
# Coefficient of variation
cv <- sd/m
# Save the dataframe
nz_tis_df <- data.frame('sample_key' = colnames(gene_counts),
                    'NZF' = nzf,
                    'GINI' = ginic,
                    'CV' = cv,
                    'SEX' = sex,
                    'TISSUE' = TISSUE1)

# Concatenate the dataframes
sample_plot <- rbind(sample_plot, nz_tis_df)
}
}

# Tidy up dataframes
###############################################
# Remove duplicates
gene_plot <- unique(gene_plot)
sample_plot <- unique(sample_plot)

# Adjust factors
sample_plot$TISSUE <- as.character(sample_plot$TISSUE)
sample_plot$TISSUE <- as.factor(sample_plot$TISSUE)

#' ## Plot QC Statistics

#+ Plot QC Statistics
################################################################################
################     Plot QC Statistics      ###################################
################################################################################

# Plot a 2D Histogram of nonzero genes across tissues
#################################################
# Color Housekeeping
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

# Plot the distributions
# TODO: Add a poisson regression line
p <- gene_plot %>%
        ggplot(aes(x = log(E_MU), y= NZF)) +
        stat_bin2d(bins = 50) +
        scale_fill_gradientn(colours=r, trans = "log") +
        xlab("Mean Expression (log10)") +
        ylab("Nonzero Fraction") +
        labs(fill = 'Log(Gene Counts)') +
        facet_wrap( ~ TISSUE + SEX) +
        theme_bw() +
        theme(strip.text = element_text(size=18),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-nonzero-fraction-tissues-faceted-2Dhistogram_steep.pdf'),width=26,height=14)
plot(p)
dev.off()

# Plot Violin plots of the fraction of Nonzero genes per tissue and sex
####################################################
p <- sample_plot %>%
        ggplot(aes(x = TISSUE, y= NZF, fill = SEX)) +
        geom_violin() +
        xlab("") +
        ylab("Nonzero Gene Fraction") +
        theme_light() +
        theme(axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-nonzero-fraction-tissues-faceted-violin_steep.pdf'),width=26,height=14)
plot(p)
dev.off()

# Plot a 2D Histogram of gini coefficients per gene across tissues
#################################################
# Color Housekeeping
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

# Plot the distributions
# TODO: Add a poisson regression line
p <- gene_plot %>%
        ggplot(aes(x = log(E_MU), y= GINI)) +
        stat_bin2d(bins = 50) +
        scale_fill_gradientn(colours=r, trans = "log") +
        xlab("Mean Expression (log10)") +
        ylab("Gini Coefficient") +
        labs(fill = 'Log(Gene Counts)') +
        facet_wrap( ~ TISSUE + SEX) +
        theme_bw() +
        theme(strip.text = element_text(size=18),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-gini-fraction-tissues-faceted-2Dhistogram_steep.pdf'),width=26,height=14)
plot(p)
dev.off()

# Plot Violin plots of the gini index per tissue and sex
####################################################
p <- sample_plot %>%
        ggplot(aes(x = TISSUE, y= GINI, fill = SEX)) +
        geom_violin() +
        xlab("") +
        ylab("Gini Index") +
        theme_light() +
        theme(axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-gini-tissues-faceted-violin_steep.pdf'),width=26,height=14)
plot(p)
dev.off()

# Plot a 2D Histogram of coefficient of variation per gene across tissues
#################################################
# Color Housekeeping
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

# Plot the distributions
# TODO: Add a poisson regression line
p <- gene_plot %>%
        ggplot(aes(x = log(E_MU), y= CV)) +
        stat_bin2d(bins = 50) +
        scale_fill_gradientn(colours=r, trans = "log") +
        xlab("Mean Expression (log10)") +
        ylab("Coefficient of Variation") +
        labs(fill = 'Log(Gene Counts)') +
        facet_wrap( ~ TISSUE + SEX) +
        theme_bw() +
        theme(strip.text = element_text(size=18),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-cv-tissues-faceted-2Dhistogram_steep.pdf'),width=26,height=14)
plot(p)
dev.off()

# Plot Violin plots of the coefficient of variation per tissue and sex
####################################################
p <- sample_plot %>%
        ggplot(aes(x = TISSUE, y= CV, fill = SEX)) +
        geom_violin() +
        xlab("") +
        ylab("Coefficient of Variation") +
        theme_light() +
        theme(axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 14),
              axis.title.y = element_text(size = 20),
              legend.title = element_blank(),
              legend.text = element_text(size = 16))
pdf(paste0(WD,'/plots/20200707_rnaseq-cv-tissues-faceted-violin_steep.pdf'),width=26,height=14)
plot(p)
dev.off()


















