'---
#' title: "PASS1A Rat Tissue: -- Apply Manova Results and SOM Clusters to All Genes in All Tissues"
#' author: "Alec Steep & Jiayu Zhang" 
#' date: "20200413"
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
#' * Model Circadian Rhythms
#'     * SIN/COS Linear Model
#'     * y = B_0 + B_1SIN(TOD) + B_2COS(TOD)
#'     * Apply Bollinger bands to Model
#' * Model Exercise Effects (68 df)
#'     * ns(x, 4)
#'     * Calculate the time (seconds) post exercise, and transform (Done)
#' * Investigate Effect of Covariants: Sex and Feeding on both models
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
#BiocManager::install("gapminder")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","tiff","splines","gam")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
counts <- DESeq2::counts
map <- purrr::map

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
        if ("lm" %!in% class(modelobject)) stop("Not an object of class 'lm' ")
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

# Create empty dataframes
manova_out <- data.frame()
cluster_out <- data.frame()

for(TISSUE in c('Spleen','Testes', 'White Adipose','Lung','Hypothalamus','Aorta','Kidney','Liver', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries')){
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
                #saveRDS(col_data, file = meta_file)
                
                # Count Data
                count_file <- paste0(WD, '/data/20200603_rnaseq-counts-pass1a-stanford-sinai-processed_steep.rds')
                #saveRDS(count_data, file = count_file)
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
        
        #' ## Collect Samples of Interest and Normalize
        
        #+ Collect Samples of Interest and Normalize
        ################################################################################
        #####     Collect Samples of Interest and Normalize      #######################
        ################################################################################
        if(TISSUE == c('Gastrocnemius')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Gastrocnemius') %>%
                        filter(Seq_batch == 'MSSM_1') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS)
        }else if(TISSUE == c('Lung')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Lung') %>%
                        filter(Seq_batch == 'MSSM_3') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS)
        }else{
                # Filter Samples (meta)
                tod_cols <- col_data %>%
                        filter(Tissue == TISSUE) %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS)
        }
        rownames(tod_cols) <- tod_cols$sample_key
        
        # Time post exercise
        tod_cols <- tod_cols %>%
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
        # Make sure to Subtract 1 hour (3600s) from "Control - IPE" groups to account for exercise effect
        tod_cols <- tod_cols %>%
                mutate(calculated.variables.deathtime_after_acute =
                               ifelse(animal.key.anirandgroup == 'Control - IPE', 
                                      calculated.variables.deathtime_after_acute - 3600,
                                      calculated.variables.deathtime_after_acute))
        tod_cols <- tod_cols %>%
                mutate(specimen.collection.t_exercise_hour_sqrt = ifelse(
                        calculated.variables.deathtime_after_acute < 0, 
                        (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)*(-1), 
                        (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)))
        
        row.names(tod_cols) <- tod_cols$sample_key
        # # Examine histograms
        # tod_cols %>%
        #         filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        #         ggplot(aes(x=calculated.variables.deathtime_after_acute)) +
        #         geom_histogram(bins = 68)
        # tod_cols %>%
        #         filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        #         ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt)) +
        #         geom_histogram(bins = 68)
        
        # Collect samples without NA values in TOD
        nona_sams <- tod_cols %>%
                filter(!is.na(specimen.collection.t_death_hour)) %>%
                filter(sample_key != OUTLIERS) %>%
                filter(!is.na(animal.registration.sex)) %>%
                select(sample_key) %>% unlist() %>% as.character()
        # Collect tissue specific counts
        tod_counts <- count_data[,nona_sams]
        
        #' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
        all(rownames(tod_cols) == colnames(tod_counts))
        
        # Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
        design = ~1 # Primary variable needs to be last.
        title = paste0('Design: ',as.character(design))
        dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                                       colData = tod_cols,
                                       design = design)
        # Reasoning from:
        #citation("PROPER")
        #dds
        #' #### We remove genes with an average sequencing depth of 10 or less
        #' Before Filtering
        # dds1
        zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% 
                nrow() %>% as.character()
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
        #' To see the reads per million for each sample
        sort(colSums(assay(dds)))/1e6
        
        # estimateSizeFactors gives us a robust estimate in sequencing depth
        dds <- estimateSizeFactors(dds)
        #' Size facotrs are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
        summary(sizeFactors(dds))
        # VST transform
        rld <- DESeq2::vst(dds, blind = F)
        # This command is redundent, but included for safety
        rs <- rowSums(counts(dds))
        
        #' #### TODO: Adjust Variance
        
        #+ Adjust Variance
        ################################################################################
        ################ Adjust Variance ###############################################
        ################################################################################
        
        # Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
        if(length(ADJ_VAR) == 2){
                assay(rld) <- limma::removeBatchEffect(assay(rld), rld[[ADJ_VAR[1]]])
                assay(rld) <- limma::removeBatchEffect(assay(rld), rld[[ADJ_VAR[2]]])
        }else if(ADJ_VAR == 'None'){
                assay(rld) <- assay(rld)
        }else{
                assay(rld) <- limma::removeBatchEffect(assay(rld), rld[[ADJ_VAR[1]]])
        }
        
        # Quick PCA Check to ensure proper variance adjustment
        pcaData <- DESeq2::plotPCA(rld, 
                                   intgroup=c('sample_key','animal.key.anirandgroup',
                                              'animal.registration.sex','animal.registration.cagenumber'),
                                   returnData=TRUE, ntop = 20000)
        percentVar <- round(100 * attr(pcaData, 'percentVar'))
        p <- pcaData %>%
                ggplot(aes_string('PC1', 'PC2', color='animal.key.anirandgroup')) +
                geom_point(size=3) +
                #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
                #geom_label_repel(aes(label=sample_misid_suspect),hjust=0, vjust=0) +
                xlab(paste0('PC1: ',percentVar[1],'% variance')) +
                ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
                coord_fixed() +
                ggtitle(paste0('PCA of ',TISSUE,' Gene Expression:\n Adjusted Variance')) +
                scale_color_manual(values=ec_colors) +
                theme(legend.title=element_blank())
        plot(p)
        
        #' #### Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group)
        
        #+ Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group)
        ################################################################################
        #Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group) #
        ################################################################################
        
        #' #### Steps in Analysis:
        #' * Orient Count Matrix by sample (group) for temporal visualization
        #' * Perform a variant test to remove genes that do not vary across groups
        #'     * Ensure that variant filter parameters are met if test is parametric
        #' * Scale and center data to have a mean of zero and a variance of 1
        #' * Perform SOM clustering
        
        # Order samples by exercise group
        group_order <- colData(rld) %>%
                as.data.frame() %>%
                arrange(factor(animal.key.anirandgroup, 
                               levels = c("Control - IPE", "Exercise - IPE", "Exercise - 0.5 hr",
                                          "Exercise - 1 hr", "Exercise - 4 hr", "Exercise - 7 hr",
                                          "Exercise - 24 hr", "Exercise - 48 hr", "Control - 7 hr")),
                        desc(specimen.collection.t_death), 
                        desc(animal.registration.sex)) %>%
                filter(animal.key.anirandgroup != "Control - 7 hr") %>%
                select(sample_key) %>% unlist() %>% as.character()
        # Order of exercise groups
        group_order2 <- colData(rld) %>%
                as.data.frame() %>%
                arrange(factor(animal.key.anirandgroup, 
                               levels = c("Control - IPE", "Exercise - IPE", "Exercise - 0.5 hr",
                                          "Exercise - 1 hr", "Exercise - 4 hr", "Exercise - 7 hr",
                                          "Exercise - 24 hr", "Exercise - 48 hr", "Control - 7 hr")),
                        desc(specimen.collection.t_death), 
                        desc(animal.registration.sex)) %>%
                filter(animal.key.anirandgroup != "Control - 7 hr") %>%
                select(animal.key.anirandgroup) %>% unlist() %>% as.character()
        
        # Arrange matrix in order of exercise and control groups
        # All expression values per gene are normalized to generate a level playing field. Means are subtracted by values and the resulting difference is divided by the standard deviation.
        som_xgroup <- assay(rld)[,group_order] %>% t() %>% data.frame()
        
        # Test for multiv
        #x <- as.numeric(MVN::mvn(som_xgroup, mvnTest = "hz")$univariateNormality$`p value`)
        #(p.adjust(x, method = "BH") > 0.05) %>% table()
        som_xgroup$animal.key.anirandgroup <- group_order2
        
        # MANOVA test
        # Generate a formula
        dependents <- colnames(som_xgroup)[colnames(som_xgroup) %!in% c("animal.key.anirandgroup")]
        form <- as.formula(paste0("cbind(",paste(dependents, collapse=","),")", "~animal.key.anirandgroup"))
        # Perform MANOVA (1-2 min per tissue)
        res_man <- manova(form, data = som_xgroup)
        res_sum <- summary.aov(res_man)
        
        # Organize the index of genes by pvalue
        df_x <- data.frame(1:(ncol(som_xgroup)-1))
        names(df_x) <- 'idx'
        pvec <- vector()
        for(nr in 1:(ncol(som_xgroup)-1)){
                pval <- res_sum[[nr]]$`Pr(>F)`[1]
                pvec <- c(pvec, pval)
        }
        df_x$pval <- pvec
        
        # Create a dataframe of genes and their manova p-value
        manova_df <- data.frame(ENSEMBL_RAT = names(som_xgroup)[df_x$idx],
                   MANOVA_PVAL = df_x$pval,
                   TISSUE = TISSUE)
        
        # Choose row index with significant p value
        c_idx <- df_x %>%
                filter(pval <= 0.05) %>%
                select(idx) %>% unlist() %>% as.numeric()
        # Select only genes that show significant variance across timepoints
        man_genes <- names(som_xgroup)[c_idx]
        som_xgroup <- som_xgroup[,c_idx]
        
        # Data should be normalized to have a mean zero and a variance of 1 (between -1 and 1)
        preproc1 <- preProcess(som_xgroup, method=c("center", "scale"))
        norm1 <- predict(preproc1, som_xgroup)
        som_mat <- t(norm1) %>% as.matrix()
        
        # Create a SOM (with som)
        set.seed(666)
        som_group <- som::som(som_mat,6,5)
        plot(som_group, ylim=c(-2,2))
        
        #som::som(som_xgroup,10,10) %>% plot()
        table(som_group$visual[,1:2])
        
        # Perform k-means clustering
        k6<-kmeans(som_mat,6)$cluster
        
        # Create an annotation for the heatmap
        ann_df <- col_data %>%
                filter(sample_key %in% colnames(som_mat)) %>%
                select(animal.key.anirandgroup)
        ann_df$animal.key.anirandgroup <- factor(ann_df$animal.key.anirandgroup, 
                                                 levels = c("Control - IPE", "Exercise - IPE", 
                                                            "Exercise - 0.5 hr","Exercise - 1 hr",
                                                            "Exercise - 4 hr", "Exercise - 7 hr",
                                                            "Exercise - 24 hr", "Exercise - 48 hr",
                                                            "Control - 7 hr"))
        ann_df <- ann_df %>% arrange(animal.key.anirandgroup) %>%
                filter(animal.key.anirandgroup != "Control - 7 hr")
        row.names(ann_df) <- colnames(som_mat)
        ann_colors = list(
                animal.key.anirandgroup = 
                        c("Control - IPE" = "steelblue1",
                          "Exercise - IPE" = "gold",
                          "Exercise - 0.5 hr" = "darkgoldenrod1",
                          "Exercise - 1 hr" = "orange",
                          "Exercise - 4 hr" = "darkorange",
                          "Exercise - 7 hr" = "darkorange2",
                          "Exercise - 24 hr" = "darkorange3",
                          "Exercise - 48 hr" = "darkorange4"))
        
        # Clusters are numbers in numerical order from 1 to 30 from bottom-left to top-right
        # Identify the circadian rhythm genes that belong to each cluster
        i <- 1
        cluster <- list()
        circ_cl <- list()
        circ_cldf <- list()
        ENSEMBL_RAT <- vector()
        for(yn in 0:4){
                for(xn in 0:5){
                        # Collect gene ids for each cluster
                        c_rn <- som_group$visual[(som_group$visual$x == xn & som_group$visual$y == yn),] %>%
                                row.names() %>% as.numeric()
                        cluster[[i]] <- som_group$data[c_rn,] %>% row.names()
                        # Collect all gene ids
                        ENSEMBL_RAT <- c(ENSEMBL_RAT, cluster[[i]])
                        # collect circadian gene ids per cluster
                        i <- i + 1
                }
        }
        
        # Create a geom_tile to mimic the som and visualize the kmeans clusters
        grad_plots <- vector('list', 6)
        kmeans_cluster <- list()
        for(kn in 1:6){
                # Kmeans cluster 1
                x <- 0:5
                y <- 0:4
                tile_df <- expand.grid(X=x,Y=y)
                # Kmeans cluster genes are in which som clusters
                kgenes <- names(k6[k6==kn])
                kmeans_cluster[[kn]] <- kgenes
                #circ_cl[[kn]] <- circ_kid %>%
                #        filter(ENSEMBL_RAT %in% kgenes) %>%
                #        select(ENSEMBL_RAT) %>% unlist(use.names = FALSE)
                # collect additional information for circadian genes in cluster
                #circ_cldf[[kn]] <- circ_kid %>%
                #        filter(ENSEMBL_RAT %in% kgenes) %>%
                #        select(NUM.TISSUE, ENSEMBL_RAT, SYMBOL_RAT)
                Z <- vector()
                for(i in 1:30){
                        Z <- c(Z, (kgenes[kgenes %in% cluster[[i]]] %>% length()))
                }
                tile_df$Z <- Z
                # Generate the kmeans plot
                # grad_plots[[kn]] <- local({
                #         kn <- kn
                #         kmsom <- ggplot(tile_df, aes(x = X, y = Y, fill =Z)) +
                #                 geom_raster(interpolate=TRUE) +
                #                 scale_fill_gradient2(low="navy", mid="white", high="red", 
                #                                      midpoint=30, limits=range(tile_df$Z),
                #                                      guide = FALSE) +
                #                 theme(axis.text        = element_blank(),
                #                       axis.ticks       = element_blank(),
                #                       axis.title       = element_blank(),
                #                       panel.background = element_blank())
                #         print(kmsom)
                # })
        }
        
        # Create a dataframe for the plot of clusters
        df1 <- data.frame(t(som_mat))
        df1$sample_key <- row.names(df1)
        df2 <- col_data %>%
                filter(sample_key %in% colnames(som_mat)) %>%
                select(animal.key.anirandgroup, specimen.collection.t_death, sample_key) %>%
                mutate(animal.key.anirandgroup = factor(animal.key.anirandgroup, 
                                                 levels = c("Control - IPE", "Exercise - IPE", 
                                                            "Exercise - 0.5 hr","Exercise - 1 hr",
                                                            "Exercise - 4 hr", "Exercise - 7 hr",
                                                            "Exercise - 24 hr", "Exercise - 48 hr",
                                                            "Control - 7 hr"))) %>%
                arrange(animal.key.anirandgroup) %>%
                filter(animal.key.anirandgroup != "Control - 7 hr")
        # Adjust the column names to be symbols for ease of plot interpretation
        df_plot <- left_join(df1, df2, by = "sample_key")
        # Melt the plot
        melt_plot <- reshape2::melt(df_plot, id.vars = c("specimen.collection.t_death", 
                                                         "animal.key.anirandgroup",
                                                         "sample_key"))
        # Adjust columns and factor levels
        melt_plot$value <- as.numeric(melt_plot$value)
        melt_plot <- melt_plot %>%
                mutate(animal.key.anirandgroup = factor(animal.key.anirandgroup, 
                                                levels = c("Control - IPE", "Exercise - IPE", 
                                                           "Exercise - 0.5 hr","Exercise - 1 hr",
                                                           "Exercise - 4 hr", "Exercise - 7 hr",
                                                           "Exercise - 24 hr", "Exercise - 48 hr")))
        names(melt_plot)[4] <- 'ENSEMBL_RAT'
        melt_plot <- melt_plot %>%
                mutate(CLUSTER = case_when(ENSEMBL_RAT %in% kmeans_cluster[[1]] ~ 1,
                                   ENSEMBL_RAT %in% kmeans_cluster[[2]] ~ 2,
                                   ENSEMBL_RAT %in% kmeans_cluster[[3]] ~ 3,
                                   ENSEMBL_RAT %in% kmeans_cluster[[4]] ~ 4,
                                   ENSEMBL_RAT %in% kmeans_cluster[[5]] ~ 5,
                                   ENSEMBL_RAT %in% kmeans_cluster[[6]] ~ 6))
        
        # Collect the median values for all points in a cluster and time point
        melt_plot$median_value <- 1000
        for(kn in c(1:6)){
                for(cohort in c('Control - IPE', 'Exercise - IPE', "Exercise - 0.5 hr",
                                "Exercise - 1 hr","Exercise - 4 hr", "Exercise - 7 hr",
                                "Exercise - 24 hr", "Exercise - 48 hr")){
                        median_val <- melt_plot %>% 
                                filter(CLUSTER == kn & animal.key.anirandgroup == cohort) %>%
                                select(value) %>% unlist() %>% as.numeric() %>% median()
                        melt_plot <- melt_plot %>%
                                mutate(median_value = if_else(
                                        CLUSTER == kn & animal.key.anirandgroup == cohort,
                                       median_val, median_value))
                }
        }
        melt_plot <- melt_plot %>%
                filter(!is.na(median_value))
        # Split-apply-combine to annotate quantile
        split_dt <- lapply(split(melt_plot, melt_plot$CLUSTER), 
                           transform, 
                           QUANTILE = as.numeric(
                                   cut(median_value, breaks = quantile(median_value, probs = seq(0,1,.25)), include.lowest = T)))
        melt_plot <- unsplit(split_dt, melt_plot$CLUSTER)
        # Generate an annotation for early/mid/late response genes up/down
        # Find peaks
        peak_df <- melt_plot %>%
                as_tibble() %>%
                group_by(CLUSTER) %>% 
                dplyr::slice(which.max(median_value)) %>%
                mutate(PEAK = animal.key.anirandgroup) %>%
                select(CLUSTER, PEAK)
        # Find bowls
        bowl_df <- melt_plot %>%
                as_tibble() %>%
                group_by(CLUSTER) %>% 
                dplyr::slice(which.min(median_value)) %>%
                mutate(BOWL = animal.key.anirandgroup) %>%
                select(CLUSTER, BOWL)
        melt_df <- melt_plot %>%
                left_join(peak_df, by = 'CLUSTER') %>%
                left_join(bowl_df, by = 'CLUSTER')
        p <- melt_df %>%
                 select(CLUSTER, BOWL, PEAK, CLUSTER) %>% unique() %>%
                 arrange(CLUSTER)
        print(p)
        # Plot the clusters from k-means
        p <- melt_df %>%
                select(median_value, CLUSTER, animal.key.anirandgroup) %>%
                unique() %>%
                ggplot(aes(x = as.integer(animal.key.anirandgroup), 
                           y = median_value),
                       group = "1") +
                geom_point(alpha = 1) +
                #geom_line(aes(group = "1")) +
                #stat_smooth(alpha = 1, se = F) +
                geom_line(stat = "smooth", alpha = 0.5, 
                          method = "lm", formula = y ~ ns(x, 5), se = T) +
                ylab("rlog Transformed Expression") +
                xlab("Exercise/Control Groups") +
                scale_x_continuous(breaks = seq(1,8,by=1),
                                   labels=c("Control - IPE", "Exercise - IPE", 
                                            "Exercise - 0.5 hr","Exercise - 1 hr",
                                            "Exercise - 4 hr", "Exercise - 7 hr",
                                            "Exercise - 24 hr", "Exercise - 48 hr")) +
                theme(legend.position = "none") +
                facet_wrap(~ CLUSTER)
        pdf(paste0(WD,'/plots/20200413_rnaseq-',TISSUE,'-kmeans-clusters_steep.pdf'),width=24,height=12)
        plot(p)
        dev.off()
        # Annotate the clusters
        melt_df <- melt_df %>%
                mutate(CLUSTER_ANN = case_when(
                        ( PEAK %in% c('Exercise - IPE', 'Exercise - 0.5 hr') &
                                BOWL %in% c('Control - IPE',  'Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Early_Up',
                        ( PEAK == 'Exercise - 1 hr' &
                                  BOWL %in% c('Control - IPE',  'Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Middle_Up',
                        ( PEAK %in% c('Exercise - 4 hr','Exercise - 7 hr') &
                                  BOWL %in% c('Control - IPE','Exercise - 0.5 hr','Exercise - IPE',  'Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Late_Up',
                        ( PEAK %in% c('Exercise - 24 hr') &
                                  BOWL %in% c('Control - IPE')) ~ 'Delayed_Up',
                        ( BOWL %in% c('Exercise - IPE', 'Exercise - 0.5 hr') &
                                  PEAK %in% c('Control - IPE', 'Exercise - 4 hr','Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Early_Down',
                        ( BOWL %in% c('Exercise - 1 hr') &
                                  PEAK %in% c('Control - IPE', 'Exercise - 4 hr', 'Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Middle_Down',
                        ( BOWL %in% c('Exercise - 4 hr', 'Exercise - 7 hr') &
                                  PEAK %in% c('Control - IPE','Exercise - IPE','Exercise - 0.5 hr', 'Exercise - 7 hr', 'Exercise - 24 hr', 'Exercise - 48 hr')) ~ 'Late_Down',
                        ( BOWL %in% c('Exercise - 24 hr') &
                                  PEAK %in% c('Control - IPE')) ~ 'Delayed_Down'))
        melt_df$CLUSTER_ANN <- factor(melt_df$CLUSTER_ANN, 
                                      levels = c('Early_Up','Middle_Up','Late_Up','Delayed_Up',
                                                 'Early_Down','Middle_Down','Late_Down','Delayed_Down'))
        # table(melt_df$CLUSTER_ANN)
        # To re-examine SOM
        # plot(som_group, ylim=c(-2,2))
        
        # Plot the circadian genes for cluster of choice
        # Show all plots they are enriched with and without in powerpoint
        #plot(grad_plots[[2]])
        #heat_plots[[1]]
        p <- melt_df %>%
                select(median_value, CLUSTER, CLUSTER_ANN, animal.key.anirandgroup) %>%
                unique() %>%
                ggplot(aes(x = as.integer(animal.key.anirandgroup), 
                           y = median_value),
                       group = "1") +
                geom_point(alpha = 1) +
                #geom_line(aes(group = "1")) +
                #stat_smooth(alpha = 1, se = F) +
                geom_line(stat = "smooth", alpha = 0.5, 
                          method = "lm", formula = y ~ ns(x, 5), se = T) +
                ylab("rlog Transformed Expression") +
                xlab("Exercise/Control Groups") +
                scale_x_continuous(breaks = seq(1,8,by=1),
                                   labels=c("Control - IPE", "Exercise - IPE", 
                                                     "Exercise - 0.5 hr","Exercise - 1 hr",
                                                     "Exercise - 4 hr", "Exercise - 7 hr",
                                                     "Exercise - 24 hr", "Exercise - 48 hr")) +
                theme(legend.position = "none") +
                facet_wrap(~ CLUSTER_ANN)
        pdf(paste0(WD,'/plots/20200413_rnaseq-',TISSUE,'-annotated-clusters_steep.pdf'),width=24,height=12)
        plot(p)
        dev.off()
        # Examine if clusters were not annotated
        p <- melt_df %>%
                select(CLUSTER, BOWL, PEAK, CLUSTER_ANN) %>% unique() %>%
                arrange(CLUSTER_ANN)
        print(p)
        # Add the tissue
        melt_df$TISSUE <- TISSUE
        melt_df$MANOVA <- TISSUE
        cluster_df <- melt_df %>%
                select(TISSUE, ENSEMBL_RAT, CLUSTER_ANN, CLUSTER)
        # Concatenate the manova dataframes
        manova_out <- rbind(manova_out, manova_df)
        # Save the final output table
        manova_file <- paste0(WD,'/data/20200413_rnaseq-tissue-manova-table_steep.txt')
        write.table(manova_out, file = manova_file,sep = '\t',row.names = F,quote = F)
        
        # Concatentate the cluster dataframes
        cluster_out <- rbind(cluster_out, cluster_df)
        # Save the final output table
        cluster_file <- paste0(WD,'/data/20200413_rnaseq-tissue-cluster-table_steep.txt')
        write.table(cluster_out, file = cluster_file,sep = '\t',row.names = F,quote = F)
}

