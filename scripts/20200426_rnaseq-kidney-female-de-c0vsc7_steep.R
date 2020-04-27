#'---
#' title: "PASS1A Rat (Female) Kidney: -- Control 0 hr vs Control 7 hr Differential Gene Expression"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200426"
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
#' * Perform a differential gene expression analysis between Control 0 and Control 7 to determine DE genes
#' * Determine if DE gene list is enriched for circadian genes in the kidney
#' * Determine if DE gene list is enriched for circadian rhythm pathways
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
#BiocManager::install("EBImage")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and Author
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

# Function to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
mouse2rat_ortho <- function(x) {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
        # Create ortholog table
        ortho_df <- getLDS(attributes=c("ensembl_gene_id","rnorvegicus_homolog_orthology_confidence"),
                           filters="ensembl_gene_id", 
                           values = x$ENSEMBL_MOUSE, 
                           mart=mart_mm_ens,
                           attributesL=c("ensembl_gene_id"), 
                           martL=mart_rn_ens) # Use biomart to get orthologs
        # Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
        ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
        ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
        ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
        names(ortho_df) <- c('ENSEMBL_MOUSE','CONFIDENCE','ENSEMBL_RAT') 
        ortho_df <- ortho_df %>%
                select(-CONFIDENCE)
        
        # Assumes that 'x' has ensembl chicken gene ID's as rownames
        # Ensure that only chicken genes appear in the matrix
        #x <- x[startsWith(rownames(x), "ENSGALG"),]
        # Assign the HUGO symbols to a new column
        x <- left_join(x, ortho_df, by = "ENSEMBL_MOUSE") %>%
                filter(!is.na(ENSEMBL_RAT)) %>%
                mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))
        x
}
#########################################################################################

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

# Adjust column objects
########################
# To factors
factor_cols <- c("labelid",
                 "vial_label",
                 "animal.registration.sex",
                 "animal.key.exlt4")
for(fc in factor_cols){
        col_data[[fc]] <- as.factor(col_data[[fc]])
}

# To Dates: 03JUL2018
date_cols <- c("acute.test.d_visit",
               "acute.test.d_start",
               "animal.familiarization.d_visit",
               "animal.familiarization.d_treadmillbegin",
               "animal.familiarization.d_treadmillcomplete",
               "animal.registration.d_visit",
               "animal.registration.d_arrive",
               "animal.registration.d_reverselight",
               "specimen.collection.d_visit",
               "animal.registration.d_birth")
for(dc in date_cols){
        col_data[[dc]] <- ymd(col_data[[dc]])
}

# To Times: 10:30:00
time_cols <- c("acute.test.t_complete",
               "specimen.collection.t_anesthesia",
               "specimen.collection.t_bloodstart",
               "specimen.collection.t_bloodstop",
               "specimen.collection.t_edtafill",
               "specimen.collection.uteruscomplete",
               "specimen.collection.t_uterusstart",
               "specimen.collection.t_uterusstop",
               "specimen.collection.t_death",
               "specimen.processing.t_collection",
               "specimen.processing.t_edtaspin",
               "specimen.processing.t_freeze",
               "acute.test.howlongshock")
for(tc in time_cols){
        col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% parse_time()
}

# Releveling factors
col_data$animal.key.anirandgroup <- as.character(col_data$animal.key.anirandgroup)
col_data$animal.key.anirandgroup <- factor(col_data$animal.key.anirandgroup, 
                                           levels = c("Control - IPE",
                                                      "Control - 7 hr",
                                                      "Exercise - IPE",
                                                      "Exercise - 0.5 hr",
                                                      "Exercise - 1 hr",
                                                      "Exercise - 4 hr",
                                                      "Exercise - 7 hr",
                                                      "Exercise - 24 hr",
                                                      "Exercise - 48 hr"))


#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script 20200413_rnaseq-kidney-female-temporal-clusters_steep.R

# Circadian Genes (Kidney)
in_file <- paste0(WD,'/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt')
circ_kid <- read.table(file=in_file, sep = '\t', header = TRUE)

#' ## Normalization for Sequencing Depth (Female Kidneys)

#+ Normalization for Sequencing Depth (Female Kidneys)
################################################################################
#####     Normalization for Sequencing Depth  (Female Kidneys)   ###############
################################################################################

#' ##### Female Kidney Samples (39 unique) from 1 batch were filtered prior to normalization
fkid <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(animal.registration.sex == 'Female') %>%
        filter(!is.na(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
fkid_counts <- count_data[,fkid] %>% as.matrix()

# SUbset kidney metadata
fkid_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(animal.registration.sex == 'Female') %>%
        filter(!is.na(animal.registration.sex))
#filter(sample_key != '90042016803_SF1')
row.names(fkid_cols) <- fkid_cols$sample_key
names(fkid_cols)
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))
colnames(fkid_cols)
#design = ~ sex + animal.key.exlt4 # Primary variable needs to be last
#design = ~ 1
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object

fkid_cols <- fkid_cols %>%
        mutate(sex = ifelse(animal.registration.sex == 'Female', 1, 0))

dds1 <- DESeqDataSetFromMatrix(countData = fkid_counts,
                               colData = fkid_cols,
                               design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
# dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 10
keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
dds2 <- dds1[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
#' TODO: Critic from Jun: Here we are removing features that have a low average expression. This may be removing important features that might have zero counts in some samples and higher counts in specific groups. Consider developing an algorithm that will account for features with expression in n or more samples.
dds2
filter_n <- nrow(dds1) - nrow(dds2) - as.numeric(zero_n)
filter_p <- filter_n/(nrow(dds1) - as.numeric(zero_n))
total_n <- nrow(dds1) - nrow(dds2)

#' ##### Note: Number of genes with average counts between zero and 1 is `r zero_n` but removing reads with less than or equal to `r reads_n` removes an additional `r filter_n` features or removes `r filter_p*100`% of the non-zero reads (total of `r total_n` features removed).
dds <- dds2
#' To see the reads per million for each sample
sort(colSums(assay(dds)))/1e6

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)
#' Size facotrs are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
summary(sizeFactors(dds))



#' ### Transformation (Female Kidneys)

#+ Transformation (Female Kidneys)
################################################################################
#####     Transformation (Female Kidneys)      #################################
################################################################################

#' ##### Log2 normalization (with size factor and pseudo count)"
log_norm_counts <- log2(counts(dds, normalized=TRUE) + 1)
#rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for( n in 1){
        start_time <- Sys.time()
        rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
#' #### Here we visualize the counts comparing log2 transform and rlog transform
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] +1), cex=0.2, main = "log2 Transform", ylim=c(0,22)) # Log2 Transform
boxplot(assay(rld)[rs > 0,], cex=0.2, main = "rlog Transform", ylim=c(0,22)) # Regularized Log Transform

#' ##### Returning to our prior plot of comparing expression of 2 samples, we can see the the variance is now more stable in the rlog transform.
mypar(1,1)
plot(assay(rld)[,c(1,2)],
     xlab = "Sample 1 rlog Gene Counts", 
     ylab = "Sample 2 rlog Gene Counts", 
     cex = 0.3, main = "rlog Transform")
