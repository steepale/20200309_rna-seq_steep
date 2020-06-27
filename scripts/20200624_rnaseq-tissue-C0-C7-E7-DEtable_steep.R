'---
#' title: 'PASS1A Rat Tissue: -- C0, C7, E7 DE Comparisons'
#' author: 'Alec Steep & Jiayu Zhang' 
#' date: '20200624'
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
#' * Perform a differential gene expression analysis between:
#'     * Control 0 and Control 7
#'     * Contol 0 and Exercise 7
#'     * Control 7 and Exercise 7
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
#BiocManager::install('EnhancedVolcano')
#install.packages('tidyverse')

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr')
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

# Capture the Date and Author
################################################################################
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- 'steep'
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
rat_mouse_ortho <- function(x, column = 'ENSEMBL_RAT', direction = 'rat2mouse') {
        # Ensure 'x' is a data.frame
        if ( class(x) != 'data.frame' ) {
                stop('x must be a data frame', class.= FALSE)
        }
        if ( column %!in% c('ENSEMBL_RAT','ENSEMBL_MOUSE') ){
                stop('column must be either ENSEMBL_RAT or ENSEMBL_MOUSE', class.= FALSE)
        }
        if ( direction %!in% c('rat2mouse', 'mouse2rat')){
                stop('direction must be either rat2mouse or mouse2rat')
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart('ensembl', dataset='mmusculus_gene_ensembl')
        mart_rn_ens = useMart('ensembl', dataset='rnorvegicus_gene_ensembl')
        # Create ortholog table
        if (direction == 'rat2mouse'){
                ortho_df <- getLDS(attributes=c('ensembl_gene_id',
                                                'mmusculus_homolog_orthology_confidence'),
                                   filters='ensembl_gene_id', 
                                   values = x[[column]], 
                                   mart=mart_rn_ens,
                                   attributesL=c('ensembl_gene_id'), 
                                   martL=mart_mm_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                #ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_RAT') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_MOUSE = mapIds(org.Mm.eg.db, ENSEMBL_MOUSE, 'SYMBOL', 'ENSEMBL'))
                
        }else{
                ortho_df <- getLDS(attributes=c('ensembl_gene_id',
                                                'rnorvegicus_homolog_orthology_confidence'),
                                   filters='ensembl_gene_id', 
                                   values = x[[column]], 
                                   mart=mart_mm_ens,
                                   attributesL=c('ensembl_gene_id'), 
                                   martL=mart_rn_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                #ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
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
                                                   'SYMBOL', 'ENSEMBL'))
        }
        # Return the output
        x
}
#########################################################################################


# Function to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
mouse2rat_ortho <- function(x) {
        # Ensure 'x' is a data.frame
        if ( class(x) != 'data.frame' ) {
                stop('x must be a data frame', class.= FALSE)
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart('ensembl', dataset='mmusculus_gene_ensembl')
        mart_rn_ens = useMart('ensembl', dataset='rnorvegicus_gene_ensembl')
        # Create ortholog table
        ortho_df <- getLDS(attributes=c('ensembl_gene_id','rnorvegicus_homolog_orthology_confidence'),
                           filters='ensembl_gene_id', 
                           values = x$ENSEMBL_MOUSE, 
                           mart=mart_mm_ens,
                           attributesL=c('ensembl_gene_id'), 
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
        #x <- x[startsWith(rownames(x), 'ENSGALG'),]
        # Assign the HUGO symbols to a new column
        x <- left_join(x, ortho_df, by = 'ENSEMBL_MOUSE') %>%
                filter(!is.na(ENSEMBL_RAT)) %>%
                mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, 'SYMBOL', 'ENSEMBL'))
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

# Sine and cosine functions
################################################################################
SIN <- function(t) sin(2*pi*t/24)
COS <- function(t) cos(2*pi*t/24)
################################################################################
# Extract p-value from linear model TODO: Adjust this code to generate permutation based p-value
################################################################################
lmp <- function (modelobject) {
        if ('lm' %!in% class(modelobject)) stop('Not an object of class lm ')
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}
################################################################################

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
# Declare Tissue

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

# Generate an empty dataframe for final output
C0_C7_E7_df <- data.frame()

for(TISSUE in c('Hypothalamus', 'Liver', 'Kidney', 'Aorta', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Lung','Ovaries','Spleen','Testes', 'White Adipose')){
print(TISSUE)
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

if(F) {
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
        factor_cols <- c('labelid',
                         'vial_label',
                         'animal.registration.sex',
                         'animal.key.exlt4',
                         'X2D_barcode',
                         'BID',
                         'Seq_flowcell_lane',
                         'Seq_flowcell_run',
                         'Seq_end_type',
                         'Lib_UMI_cycle_num',
                         'pid',
                         'acute.test.staffid',
                         'acute.test.siteid',
                         'acute.test.versionnbr',
                         'acute.test.contactshock',
                         'animal.familiarization.staffid',
                         'animal.familiarization.siteid',
                         'animal.familiarization.versionnbr',
                         'animal.familiarization.compliant',
                         'animal.key.protocol',
                         'animal.key.agegroup',
                         'animal.key.batch',
                         'animal.key.intervention',
                         'animal.key.sitename',
                         'animal.registration.staffid',
                         'animal.registration.siteid',
                         'animal.registration.versionnbr',
                         'animal.registration.ratid',
                         'animal.registration.batchnumber',
                         'specimen.collection.bloodcomplete',
                         'specimen.collection.bloodtechid',
                         'specimen.collection.uterustype',
                         'specimen.collection.uterustechid',
                         'specimen.collection.deathtype',
                         'specimen.processing.versionnbr',
                         'specimen.processing.siteid',
                         'bid',
                         'specimen.processing.samplenumber',
                         'specimen.processing.techid',
                         'barcode',
                         'shiptositeid',
                         'receivedcas',
                         'receivestatuscas')
        for(fc in factor_cols){
                col_data[[fc]] <- as.factor(col_data[[fc]])
        }
        
        # To Dates: 03JUL2018
        date_cols <- c('acute.test.d_visit',
                       'acute.test.d_start',
                       'animal.familiarization.d_visit',
                       'animal.familiarization.d_treadmillbegin',
                       'animal.familiarization.d_treadmillcomplete',
                       'animal.registration.d_visit',
                       'animal.registration.d_arrive',
                       'animal.registration.d_reverselight',
                       'specimen.collection.d_visit',
                       'animal.registration.d_birth',
                       'Seq_date')
        for(dc in date_cols){
                col_data[[dc]] <- ymd(col_data[[dc]])
        }
        
        # From Dates: 2/14/2019
        date_cols <- c('RNA_extr_date',
                       'Lib_prep_date')
        for(dc in date_cols){
                col_data[[dc]] <- mdy(col_data[[dc]])
        }
        
        # To Times: 10:30:00
        time_cols <- c('acute.test.t_complete',
                       'specimen.collection.t_anesthesia',
                       'specimen.collection.t_bloodstart',
                       'specimen.collection.t_bloodstop',
                       'specimen.collection.t_edtafill',
                       'specimen.collection.uteruscomplete',
                       'specimen.collection.t_uterusstart',
                       'specimen.collection.t_uterusstop',
                       'specimen.collection.t_death',
                       'specimen.processing.t_collection',
                       'specimen.processing.t_edtaspin',
                       'specimen.processing.t_freeze',
                       'acute.test.howlongshock',
                       'acute.test.t_start')
        for(tc in time_cols){
                col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% parse_time()
        }
        
        # Releveling factors
        col_data$animal.key.anirandgroup <- as.character(col_data$animal.key.anirandgroup)
        col_data$animal.key.anirandgroup <- factor(col_data$animal.key.anirandgroup,
                                                   levels = ec_levels)
        
        # Create a variable for time post exercise
        col_data <- col_data %>%
                mutate(specimen.collection.t_exercise_hour = case_when(
                        animal.key.anirandgroup == 'Control - IPE' ~ -1,
                        animal.key.anirandgroup == 'Control - 7 hr' ~ 7,
                        animal.key.anirandgroup == 'Exercise - IPE' ~ 0,
                        animal.key.anirandgroup == 'Exercise - 0.5 hr' ~ 0.5,
                        animal.key.anirandgroup == 'Exercise - 1 hr' ~ 1,
                        animal.key.anirandgroup == 'Exercise - 4 hr' ~ 4,
                        animal.key.anirandgroup == 'Exercise - 7 hr' ~ 7,
                        animal.key.anirandgroup == 'Exercise - 24 hr' ~ 24,
                        animal.key.anirandgroup == 'Exercise - 48 hr' ~ 48))
        
        # Take the absolute value of the square root of seconds post exercise (consider negative numbers)
        # Make sure to Subtract 1 hour (3600s) from 'Control - IPE' groups to account for exercise effect
        col_data <- col_data %>%
                mutate(calculated.variables.deathtime_after_acute =
                               ifelse(animal.key.anirandgroup == 'Control - IPE',
                                      calculated.variables.deathtime_after_acute - 3600,
                                      calculated.variables.deathtime_after_acute))
        col_data <- col_data %>%
                mutate(specimen.collection.t_exercise_hour_sqrt = ifelse(
                        calculated.variables.deathtime_after_acute < 0,
                        (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)*(-1),
                        (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)))
        row.names(col_data) <- col_data$sample_key
        
        # Examine histograms
        col_data %>%
                filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
                ggplot(aes(x=calculated.variables.deathtime_after_acute)) +
                geom_histogram(bins = 68)
        col_data %>%
                filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
                ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt)) +
                geom_histogram(bins = 68)
        
        # Save data as an R objects
        # ################################################################################
        # To determine object size
        sl <- object.size(count_data)
        print(sl, units = 'auto')
        # Meta Data
        meta_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
        saveRDS(col_data, file = meta_file)
        
        # Count Data
        count_file <- paste0(WD, '/data/20200603_rnaseq-counts-pass1a-stanford-sinai-processed_steep.rds')
        saveRDS(count_data, file = count_file)
}

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

#' ### Differential Gene Expression: Control 0 hr vs Control 7 hr

#+ Differential Gene Expression: Control 0 hr vs Control 7 hr
################################################################################
###########Differential Gene Expression: Control 0 hr vs Control 7 hr  #########
################################################################################
if(TISSUE == c('Gastrocnemius')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr'))
}else if(TISSUE == c('Lung')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Lung') %>%
                filter(Seq_batch == 'MSSM_3') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr'))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr'))
}
rownames(tod_cols) <- tod_cols$sample_key
# Collect tissue specific counts
tod_counts <- count_data[,tod_cols$sample_key]

# Generate an annotation specific to this analysis
# Make sure the "control" or "untreated" level is first
tod_cols <- tod_cols %>%
        mutate(controls_binary = factor(case_when(animal.key.anirandgroup == 'Control - IPE' ~ 'Control_IPE',
                                           animal.key.anirandgroup == 'Control - 7 hr' ~ 'Control_7hr'),
               levels = c("Control_IPE","Control_7hr")))
row.names(tod_cols) <- tod_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
if(length(ADJ_VAR) == 2){
        design = as.formula(paste0('~',ADJ_VAR[1],' + ',ADJ_VAR[2],' + controls_binary'))
}else if(ADJ_VAR == 'None'){
        design = as.formula('~controls_binary')
}else{
        design = as.formula(paste0('~',ADJ_VAR[1],' + controls_binary'))
}

(title = paste0('Design: ',as.character(design)))
dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = design)

zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 5
keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
dds2 <- dds1[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
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
# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
# Generate a DESeq2 Model
dds <- DESeq(dds)
# Generate a results table and sort by fdr
res <- results(dds, alpha = 0.05,lfcThreshold=0.25)
res <- res[order(res$padj),]
# Generate a summary of the results
summary(res)
res_df <- as.data.frame(res)
# DE dataframe
C0C7_df <- data.frame(TISSUE = TISSUE,
                      ENSEMBL_RAT = row.names(res_df),
                      C0C7_log2FoldChange = res_df$log2FoldChange,
                      C0C7_padj = res_df$padj)

#' ### Differential Gene Expression: Control 0 hr vs Exercise 7 hr

#+ Differential Gene Expression: Control 0 hr vs Exercise 7 hr
################################################################################
###########Differential Gene Expression: Control 0 hr vs Exercise 7 hr  #########
################################################################################
if(TISSUE == c('Gastrocnemius')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Exercise - 7 hr'))
}else if(TISSUE == c('Lung')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Lung') %>%
                filter(Seq_batch == 'MSSM_3') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Exercise - 7 hr'))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - IPE', 'Exercise - 7 hr'))
}
rownames(tod_cols) <- tod_cols$sample_key
# Collect tissue specific counts
tod_counts <- count_data[,tod_cols$sample_key]
# Generate an annotation specific to this analysis
# Make sure the "control" or "untreated" level is first
tod_cols <- tod_cols %>%
        mutate(controls_binary = factor(case_when(animal.key.anirandgroup == 'Control - IPE' ~ 'Control_IPE',
                                                  animal.key.anirandgroup == 'Exercise - 7 hr' ~ 'Exercise_7hr'),
                                        levels = c("Control_IPE","Exercise_7hr")))
row.names(tod_cols) <- tod_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))
# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
if(length(ADJ_VAR) == 2){
        design = as.formula(paste0('~',ADJ_VAR[1],' + ',ADJ_VAR[2],' + controls_binary'))
}else if(ADJ_VAR == 'None'){
        design = as.formula('~controls_binary')
}else{
        design = as.formula(paste0('~',ADJ_VAR[1],' + controls_binary'))
}

(title = paste0('Design: ',as.character(design)))
dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = design)
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 5
keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
dds2 <- dds1[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
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
# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
# Generate a DESeq2 Model
dds <- DESeq(dds)
# Generate a results table and sort by fdr
res <- results(dds, alpha = 0.05,lfcThreshold=0.25)
res <- res[order(res$padj),]
# Generate a summary of the results
summary(res)
res_df <- as.data.frame(res)
# DE dataframe
C0E7_df <- data.frame(TISSUE = TISSUE,
                      ENSEMBL_RAT = row.names(res_df),
                      C0E7_log2FoldChange = res_df$log2FoldChange,
                      C0E7_padj = res_df$padj)

#' ### Differential Gene Expression: Control 7 hr vs Exercise 7 hr

#+ Differential Gene Expression: Control 7 hr vs Exercise 7 hr
################################################################################
###########Differential Gene Expression: Control 7 hr vs Exercise 7 hr  #########
################################################################################
if(TISSUE == c('Gastrocnemius')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - 7 hr', 'Exercise - 7 hr'))
}else if(TISSUE == c('Lung')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Lung') %>%
                filter(Seq_batch == 'MSSM_3') %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - 7 hr', 'Exercise - 7 hr'))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex)) %>%
                filter(sample_key %!in% OUTLIERS) %>%
                filter(animal.key.anirandgroup %in% c('Control - 7 hr', 'Exercise - 7 hr'))
}
rownames(tod_cols) <- tod_cols$sample_key
# Collect tissue specific counts
tod_counts <- count_data[,tod_cols$sample_key]
# Generate an annotation specific to this analysis
# Make sure the "control" or "untreated" level is first
tod_cols <- tod_cols %>%
        mutate(controls_binary = factor(case_when(animal.key.anirandgroup == 'Control - 7 hr' ~ 'Control_7hr',
                                                  animal.key.anirandgroup == 'Exercise - 7 hr' ~ 'Exercise_7hr'),
                                        levels = c("Control_7hr","Exercise_7hr")))
row.names(tod_cols) <- tod_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
if(length(ADJ_VAR) == 2){
        design = as.formula(paste0('~',ADJ_VAR[1],' + ',ADJ_VAR[2],' + controls_binary'))
}else if(ADJ_VAR == 'None'){
        design = as.formula('~controls_binary')
}else{
        design = as.formula(paste0('~',ADJ_VAR[1],' + controls_binary'))
}

(title = paste0('Design: ',as.character(design)))
dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = design)
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 5
keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
dds2 <- dds1[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
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
# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
# Generate a DESeq2 Model
dds <- DESeq(dds)
# Generate a results table and sort by fdr
res <- results(dds, alpha = 0.05,lfcThreshold=0.25)
res <- res[order(res$padj),]
# Generate a summary of the results
summary(res)
res_df <- as.data.frame(res)
# DE dataframe
C7E7_df <- data.frame(TISSUE = TISSUE,
                      ENSEMBL_RAT = row.names(res_df),
                      C7E7_log2FoldChange = res_df$log2FoldChange,
                      C7E7_padj = res_df$padj)

#' ### Join Results by Tissue and Gene

#+ Join Results by Tissue and Gene
################################################################################
########### Join Results by Tissue and Gene  ###################################
################################################################################

# Tissue specific
C0C7C0E7_df <- full_join(C0C7_df,C0E7_df, by = c('TISSUE','ENSEMBL_RAT'))
C0C7E7_df <- full_join(C0C7C0E7_df,C7E7_df, by = c('TISSUE','ENSEMBL_RAT'))

# Concatenate the dataframes
C0_C7_E7_df <- rbind(C0_C7_E7_df, C0C7E7_df)
}

# Save the final output table
table_file <- paste0(WD,'/data/20200624_rnaseq-tissue-C0-C7-E7-table_steep.txt')
#write.table(C0_C7_E7_df, file = table_file,sep = '\t',row.names = F,quote = F)




