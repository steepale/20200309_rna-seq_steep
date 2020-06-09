#'---
#' title: "PASS1A Rat All Tissue: -- PCA Analysis"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200503"
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
#' * Perform PCA on Tissues, collectively and Individually
#' * Determine which metadata variables correlate with PCs
#' * Perform a circadian gene subset of each PCA plot for each tissue
#' * Group tissues influenced by circadian rhythm and those that are not
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
#BiocManager::install("enrichR")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots")
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

# Functions to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
rat_mouse_ortho <- function(x, column = "ENSEMBL_RAT", direction = 'rat2mouse') {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        if ( column %!in% c('ENSEMBL_RAT','ENSEMBL_MOUSE') ){
                stop("'column' must be either 'ENSEMBL_RAT' or 'ENSEMBL_MOUSE'", class.= FALSE)
        }
        if ( direction %!in% c('rat2mouse', 'mouse2rat')){
                stop("'direction' must be either 'rat2mouse' or 'mouse2rat'")
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
        # Create ortholog table
        if (direction == 'rat2mouse'){
                ortho_df <- getLDS(attributes=c("ensembl_gene_id",
                                                "mmusculus_homolog_orthology_confidence"),
                                   filters="ensembl_gene_id", 
                                   values = x[[column]], 
                                   mart=mart_rn_ens,
                                   attributesL=c("ensembl_gene_id"), 
                                   martL=mart_mm_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_RAT') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_MOUSE = mapIds(org.Mm.eg.db, ENSEMBL_MOUSE, "SYMBOL", "ENSEMBL"))
                
        }else{
                ortho_df <- getLDS(attributes=c("ensembl_gene_id",
                                                "rnorvegicus_homolog_orthology_confidence"),
                                   filters="ensembl_gene_id", 
                                   values = x[[column]], 
                                   mart=mart_mm_ens,
                                   attributesL=c("ensembl_gene_id"), 
                                   martL=mart_rn_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_MOUSE','CONFIDENCE','ENSEMBL_RAT') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the HUGO symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_MOUSE') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, 
                                                   ENSEMBL_RAT, 
                                                   "SYMBOL", "ENSEMBL"))
        }
        # Return the output
        x
}
#########################################################################################

# Function to take lowercase strings and convert the first letter to uppercase
################################################################################
firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
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

# Adjust column objects
########################
# To factors
factor_cols <- c("labelid",
                 "vial_label",
                 "animal.registration.sex",
                 "animal.key.exlt4",
                 "X2D_barcode",
                 "BID",
                 "Seq_flowcell_lane",
                 "Seq_flowcell_run",
                 "Seq_end_type",
                 "Lib_UMI_cycle_num",
                 "pid",
                 "acute.test.staffid",
                 "acute.test.siteid",
                 "acute.test.versionnbr",
                 "acute.test.contactshock",
                 "animal.familiarization.staffid",
                 "animal.familiarization.siteid",
                 "animal.familiarization.versionnbr",
                 "animal.familiarization.compliant",
                 "animal.key.protocol",
                 "animal.key.agegroup",
                 "animal.key.batch",
                 "animal.key.intervention",
                 "animal.key.sitename",
                 "animal.registration.staffid",
                 "animal.registration.siteid",
                 "animal.registration.versionnbr",
                 "animal.registration.ratid",
                 "animal.registration.batchnumber",
                 "specimen.collection.bloodcomplete",
                 "specimen.collection.bloodtechid",
                 "specimen.collection.uterustype",
                 "specimen.collection.uterustechid",
                 "specimen.collection.deathtype",
                 "specimen.processing.versionnbr",
                 "specimen.processing.siteid",
                 "bid",
                 "specimen.processing.samplenumber",
                 "specimen.processing.techid",
                 "barcode",
                 "shiptositeid",
                 "receivedcas",
                 "receivestatuscas")
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
               "animal.registration.d_birth",
               "Seq_date")
for(dc in date_cols){
        col_data[[dc]] <- ymd(col_data[[dc]])
}

# From Dates: 2/14/2019
date_cols <- c("RNA_extr_date",
               "Lib_prep_date")
for(dc in date_cols){
        col_data[[dc]] <- mdy(col_data[[dc]])
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
               "acute.test.howlongshock",
               "acute.test.t_start")
for(tc in time_cols){
        col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% parse_time()
}

# Set a vector for Exercise/Control Levels and Colors
ec_levels <- c("Exercise - IPE",
               "Exercise - 0.5 hr",
               "Exercise - 1 hr",
               "Exercise - 4 hr",
               "Exercise - 7 hr",
               "Exercise - 24 hr",
               "Exercise - 48 hr",
               "Control - IPE",
               "Control - 7 hr")
ec_colors <- c("gold",
               "darkgoldenrod1",
               "orange",
               "darkorange",
               "darkorange2",
               "darkorange3",
               "darkorange4",
               "steelblue1",
               "steelblue4")
# Releveling factors
col_data$animal.key.anirandgroup <- as.character(col_data$animal.key.anirandgroup)
col_data$animal.key.anirandgroup <- factor(col_data$animal.key.anirandgroup, 
                                           levels = ec_levels)

#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script 20200503_rnaseq-all-tissues-circadian-gene-sets_steep.R

# Load Circadian Genes (All tissues in gene set)
yan_tissues <- c('SCN','LIV','AOR','SKM','HAT','ADG','BAT',
                'WAT','BON','PFR','WB','ATR','VEN')
circ_tis <- list()
for(tis in yan_tissues){
        in_file <- paste0(WD,
                          '/data/20200503_rnaseq-circadian-',
                          tis,'-mouse-rat-ortho_steep-yan.txt')
        circ_tis[[tis]] <- read.table(file=in_file, sep = '\t', header = TRUE)
}


#' ## PCA (All Tissues)

#+ PCA (All Tissues)
################################################################################
###########     PCA (All Tissues)   ############################################
################################################################################

#' ##### Female Kidney Samples (39 unique) from 1 batch were filtered prior to normalization
#' # TODO: Create a seperate script that explains why sample 90109015902_SN1 is considered in outlier
sams <- col_data %>%
        #filter(Tissue == 'Kidney') %>%
        #filter(sample_key != '90109015902_SN1') %>%
        #filter(!is.na(animal.registration.sex)) %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
sams_counts <- count_data[,sams] %>% as.matrix()

# SUbset kidney metadata
sams_cols <- col_data
        #filter(Tissue == 'Kidney') %>%
        #filter(sample_key != '90109015902_SN1') %>%
        #filter(animal.registration.sex == 'Female') %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        #filter(!is.na(animal.registration.sex))

row.names(sams_cols) <- sams_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(sams_cols) == colnames(sams_counts))

design = ~1 # Primary variable needs to be last.
title = paste0('Design: ',as.character(design))
dds1 <- DESeqDataSetFromMatrix(countData = sams_counts,
                               colData = sams_cols,
                               design = design)
# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
# dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 1
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

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)

rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for(n in 1){
        start_time <- Sys.time()
        #rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

# Determine which variables are associated with PC's
pca_pcs <-prcomp(t(assay(rld)), center = FALSE, scale. = FALSE)

# Systematically go through metadata variables and examine correlations (Adjusted R Sq)
pc_cor <- function(x,y){
        if(length(unique(y[!is.na(y)])) >=2){
                fit <- lm(x~y)
                summary(fit)$adj.r.squared
        }else{
                NA
        }
}

# Create a dataframe of adjusted R squared values
pc_cor_df <- data.frame(matrix(ncol = length(names(sams_cols)), nrow = ncol(pca_pcs$x[,1:4])))
names(pc_cor_df) <- names(sams_cols)
i <- 1
for(c in names(col_data)){
        print(i)
        pc_cor_df[[i]] <- apply(pca_pcs$x[,1:4],2, pc_cor, sams_cols[[c]])
        i = i + 1
}

# Select only those variables with a relatively high correlation
keeps <- vector()
removes <- vector()
for(C in names(pc_cor_df)){
        if(is.na(max(pc_cor_df[[C]]))){
                removes <- c(removes,C)
        }
        else if(max(pc_cor_df[[C]]) >= 0.2){
                keeps <- c(keeps,C)
        }
}

pc_cor_df <- pc_cor_df %>% select(all_of(keeps))
pc_cor_df$PC <- 1:ncol(pca_pcs$x[,1:4])

pc_cor_df2 <- pivot_longer(pc_cor_df, names(pc_cor_df)[names(pc_cor_df) != 'PC'],
                          names_to = "Condition",
                          values_to = "Adjusted_R_Sq")
pc_cor_df2$Condition <- factor(pc_cor_df2$Condition)
# Select only the first PCs of interest
pc_cor_plot <- pc_cor_df2 %>%
        filter(PC %in% c(1:2))
pc_cor_plot <- pc_cor_plot %>%
        filter(Adjusted_R_Sq > 0.2)
pc_cor_plot$Condition <- as.character(pc_cor_plot$Condition) %>% as.factor()

#display First PCs
ggplot(pc_cor_plot,aes(x = PC, y = Adjusted_R_Sq)) +
        geom_line(aes(col = Condition), show.legend = T) +
        #geom_text_repel(aes(label = Condition)) +
        scale_x_continuous(breaks = 1:4)

# Collect PCs of Interest
pcoi <- pc_cor_plot %>%
        filter(PC == 1) %>%
        filter(Adjusted_R_Sq > 0.2)

# PCA of tissues
data.frame(pcoi)
mypar()
# Most revealing:
voi <- c("Seq_date", 
         "Seq_batch", 
         "Lib_prep_date", 
         "RNA_extr_date", 
         "Seq_flowcell_run", 
         "Lib_RNA_conc")

# Create a table of interesting correlations (manually investigated)
pcoi <- pcoi %>% filter(Condition %in% voi) %>%
        arrange(desc(Adjusted_R_Sq))

# Plot each of the variables of interest (manually)
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("Tissue",
                                      "Seq_batch",
                                      "RNA_extr_plate_ID",
                                      "RNA_extr_date",
                                      "Lib_prep_date",
                                      "Lib_batch_ID",
                                      "Lib_frag_size..bp.",
                                      "Seq_date",
                                      "Seq_machine_ID",
                                      "Seq_flowcell_ID",
                                      "Seq_flowcell_run",
                                      "Lib_RNA_conc",
                                      "Lib_DNA_conc"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, 
                    color=Seq_date)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: All Tissues") +
        guides(color=guide_legend(title="Seq_date"))

#' ## PCA (Kidney)

#+ PCA (Kidney)
################################################################################
###########     PCA (Kidney)   ############################################
################################################################################

#' # TODO: Create a seperate script that explains why sample 90109015902_SN1 is considered in outlier
sams <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(sample_key != '90109015902_SN1') %>%
        filter(!is.na(animal.registration.sex)) %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
sams_counts <- count_data[,sams] %>% as.matrix()

# SUbset kidney metadata
sams_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(sample_key != '90109015902_SN1') %>%
        #filter(animal.registration.sex == 'Female') %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>
        filter(!is.na(animal.registration.sex))

row.names(sams_cols) <- sams_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(sams_cols) == colnames(sams_counts))

design = ~1 # Primary variable needs to be last.
title = paste0('Design: ',as.character(design))
dds1 <- DESeqDataSetFromMatrix(countData = sams_counts,
                               colData = sams_cols,
                               design = design)
# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
# dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 1
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

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)

rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for(n in 1){
        start_time <- Sys.time()
        #rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

# Adjust for Between Sex Variance
# "To adjust for batch effects, we median- centered the expression levels of each transcript within each batch and confirmed, using the correlation matrices, that the batch effects were removed after the adjustment." 
#~ Li, J. Z. et al. Circadian patterns of gene expression in the human brain and disruption in major depressive disorder. Proc. Natl. Acad. Sci. U. S. A. 110, 9950–9955 (2013).

# Here we have 2 Groups: Control - IPE and Control 7 hr; we'll median center these groups to combine the sexes.

M_samples <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Male') %>%
        filter(sample_key != '90109015902_SN1') %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
F_samples <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Female') %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Select the counts
M_counts <- assay(rld[, M_samples])
F_counts <- assay(rld[, F_samples])

# Median Center data
# Collects median of each row, then subtracts by row medians
M_medians <- apply(M_counts,1,median)
M_centered <- M_counts - M_medians
F_medians <- apply(F_counts,1,median)
F_centered <- F_counts - F_medians
counts_centered <- cbind(M_centered, F_centered)
counts_centered <- counts_centered[, colnames(assay(rld))]
assay(rld) <- counts_centered

# Determine which variables are associated with PC's
pca_pcs <-prcomp(t(assay(rld)), center = FALSE, scale. = FALSE)

# Systematically go through metadata variables and examine correlations (Adjusted R Sq)
pc_cor <- function(x,y){
        if(length(unique(y[!is.na(y)])) >=2){
                fit <- lm(x~y)
                summary(fit)$adj.r.squared
        }else{
                NA
        }
}

# Create a dataframe of adjusted R squared values
pc_cor_df <- data.frame(matrix(ncol = length(names(sams_cols)), nrow = ncol(pca_pcs$x[,1:4])))
names(pc_cor_df) <- names(sams_cols)
i <- 1
for(c in names(col_data)){
        print(i)
        pc_cor_df[[i]] <- apply(pca_pcs$x[,1:4],2, pc_cor, sams_cols[[c]])
        i = i + 1
}

# Select only those variables with a relatively high correlation
keeps <- vector()
removes <- vector()
for(C in names(pc_cor_df)){
        if(is.na(max(pc_cor_df[[C]]))){
                removes <- c(removes,C)
        }
        else if(max(pc_cor_df[[C]]) >= 0.2){
                keeps <- c(keeps,C)
        }
}

pc_cor_df <- pc_cor_df %>% select(all_of(keeps))
pc_cor_df$PC <- 1:ncol(pca_pcs$x[,1:4])

pc_cor_df2 <- pivot_longer(pc_cor_df, names(pc_cor_df)[names(pc_cor_df) != 'PC'],
                           names_to = "Condition",
                           values_to = "Adjusted_R_Sq")
pc_cor_df2$Condition <- factor(pc_cor_df2$Condition)
# Select only the first PCs of interest
pc_cor_plot <- pc_cor_df2 %>%
        filter(PC %in% c(1:2))
pc_cor_plot <- pc_cor_plot %>%
        filter(Adjusted_R_Sq > 0.2)
pc_cor_plot$Condition <- as.character(pc_cor_plot$Condition) %>% as.factor()

#display First PCs
ggplot(pc_cor_plot,aes(x = PC, y = Adjusted_R_Sq)) +
        geom_line(aes(col = Condition), show.legend = T) +
        #geom_text_repel(aes(label = Condition)) +
        scale_x_continuous(breaks = 1:4)

# Collect PCs of Interest
pcoi <- pc_cor_plot %>%
        filter(PC == 2) %>%
        filter(Adjusted_R_Sq > 0.2) %>%
        arrange(desc(Adjusted_R_Sq))

# PCA of tissues
data.frame(pcoi)
mypar()
# Most revealing:
voi <- c("animal.key.anirandgroup", 
         "specimen.collection.t_death_bins.type",
         "specimen.collection.t_death_hour")

# Create a table of interesting correlations (manually investigated)
pcoi <- pcoi %>% filter(Condition %in% voi) %>%
        arrange(desc(Adjusted_R_Sq))

# Plot each of the variables of interest (manually)
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "specimen.collection.t_death_bins.type",
                                      "specimen.collection.t_death_hour"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, 
                    color=specimen.collection.t_death_hour)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Kindey (M & F)") +
        guides(color=guide_legend(title="specimen.collection.t_death_hour"))
        #scale_color_manual(values=ec_colors)

#' ## PCA (Liver)

#+ PCA (Liver)
################################################################################
###########     PCA (Liver)   ############################################
################################################################################

#' # TODO: Create a seperate script that explains why sample 90042016803_SF1 is considered in outlier
sams <- col_data %>%
        filter(Tissue == 'Liver') %>%
        filter(sample_key != '90042016803_SF1') %>%
        filter(!is.na(animal.registration.sex)) %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
sams_counts <- count_data[,sams] %>% as.matrix()

# SUbset kidney metadata
sams_cols <- col_data %>%
        filter(Tissue == 'Liver') %>%
        filter(sample_key != '90042016803_SF1') %>%
        #filter(animal.registration.sex == 'Female') %>%
        #filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>
        filter(!is.na(animal.registration.sex))

row.names(sams_cols) <- sams_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(sams_cols) == colnames(sams_counts))

design = ~1 # Primary variable needs to be last.
title = paste0('Design: ',as.character(design))
dds1 <- DESeqDataSetFromMatrix(countData = sams_counts,
                               colData = sams_cols,
                               design = design)
# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
# dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 1
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

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)

rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for(n in 1){
        start_time <- Sys.time()
        #rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

# Adjust for Between Sex Variance
# "To adjust for batch effects, we median- centered the expression levels of each transcript within each batch and confirmed, using the correlation matrices, that the batch effects were removed after the adjustment." 
#~ Li, J. Z. et al. Circadian patterns of gene expression in the human brain and disruption in major depressive disorder. Proc. Natl. Acad. Sci. U. S. A. 110, 9950–9955 (2013).

# Here we have 2 Groups: Control - IPE and Control 7 hr; we'll median center these groups to combine the sexes.

M_samples <- col_data %>%
        filter(Tissue == 'Liver') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Male') %>%
        filter(sample_key != '90042016803_SF1') %>%
        select(sample_key) %>% unlist() %>% as.character()
F_samples <- col_data %>%
        filter(Tissue == 'Liver') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(sample_key != '90042016803_SF1') %>%
        filter(animal.registration.sex == 'Female') %>%
        select(sample_key) %>% unlist() %>% as.character()
# Select the counts
M_counts <- assay(rld[, M_samples])
F_counts <- assay(rld[, F_samples])

# Median Center data
# Collects median of each row, then subtracts by row medians
M_medians <- apply(M_counts,1,median)
M_centered <- M_counts - M_medians
F_medians <- apply(F_counts,1,median)
F_centered <- F_counts - F_medians
counts_centered <- cbind(M_centered, F_centered)
counts_centered <- counts_centered[, colnames(assay(rld))]
assay(rld) <- counts_centered

# Determine which variables are associated with PC's
pca_pcs <-prcomp(t(assay(rld)), center = FALSE, scale. = FALSE)

# Systematically go through metadata variables and examine correlations (Adjusted R Sq)
pc_cor <- function(x,y){
        if(length(unique(y[!is.na(y)])) >=2){
                fit <- lm(x~y)
                summary(fit)$adj.r.squared
        }else{
                NA
        }
}

# Create a dataframe of adjusted R squared values
pc_cor_df <- data.frame(matrix(ncol = length(names(sams_cols)), nrow = ncol(pca_pcs$x[,1:4])))
names(pc_cor_df) <- names(sams_cols)
i <- 1
for(c in names(col_data)){
        print(i)
        pc_cor_df[[i]] <- apply(pca_pcs$x[,1:4],2, pc_cor, sams_cols[[c]])
        i = i + 1
}

# Select only those variables with a relatively high correlation
keeps <- vector()
removes <- vector()
for(C in names(pc_cor_df)){
        if(is.na(max(pc_cor_df[[C]]))){
                removes <- c(removes,C)
        }
        else if(max(pc_cor_df[[C]]) >= 0.2){
                keeps <- c(keeps,C)
        }
}

pc_cor_df <- pc_cor_df %>% select(all_of(keeps))
pc_cor_df$PC <- 1:ncol(pca_pcs$x[,1:4])

pc_cor_df2 <- pivot_longer(pc_cor_df, names(pc_cor_df)[names(pc_cor_df) != 'PC'],
                           names_to = "Condition",
                           values_to = "Adjusted_R_Sq")
pc_cor_df2$Condition <- factor(pc_cor_df2$Condition)

# Select only the first PCs of interest
pc_cor_plot <- pc_cor_df2 %>%
        filter(PC %in% c(1:2))
pc_cor_plot <- pc_cor_plot %>%
        filter(Adjusted_R_Sq > 0.2)
pc_cor_plot$Condition <- as.character(pc_cor_plot$Condition) %>% as.factor()

#display First PCs
ggplot(pc_cor_plot,aes(x = PC, y = Adjusted_R_Sq)) +
        geom_line(aes(col = Condition), show.legend = T) +
        #geom_text_repel(aes(label = Condition)) +
        scale_x_continuous(breaks = 1:4)

# Collect PCs of Interest
pcoi <- pc_cor_plot %>%
        filter(PC == 1) %>%
        filter(Adjusted_R_Sq > 0.2) %>%
        arrange(desc(Adjusted_R_Sq))

# PCA of tissues
data.frame(pcoi)
mypar()
# Most revealing:
voi <- c("animal.key.anirandgroup", 
         "specimen.collection.t_death_bins.type",
         "specimen.collection.t_death_hour")


# Create a table of interesting correlations (manually investigated)
#pcoi <- pcoi %>% filter(Condition %in% voi) %>%
#        arrange(desc(Adjusted_R_Sq))

# Plot each of the variables of interest (manually)
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "specimen.collection.t_death_bins.type",
                                      "specimen.collection.t_death_hour",
                                      "animal.registration.sex",
                                      "sample_key",
                                      "acute.test.t_complete",
                                      "acute.test.t_start",
                                      "animal.key.sacrificetime"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, 
                    color=animal.key.anirandgroup)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Liver (M & F)") +
        guides(color=guide_legend(title="animal.key.anirandgroup")) +
        scale_color_manual(values=ec_colors)


col_data %>% 
        filter(animal.key.anirandgroup == 'Exercise - 7 hr') %>%
        select(specimen.collection.t_death)


