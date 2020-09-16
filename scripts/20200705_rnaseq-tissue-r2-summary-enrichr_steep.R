'---
#' title: "PASS1A Rat Tissue: -- Summary R2 Pathway: EDA"
#' author: "Alec Steep & Jiayu Zhang" 
#' date: "20200705"
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
#' * Perform Pathway Enrichment Analysis on Gene of Interest
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

# Load dependencies
pacs...man <- c("tidyverse",'ggplot2',"GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","tiff","splines","gam")
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
                # ortho_df <- ortho_df %>%
                #         select(-CONFIDENCE)
                # Assign the symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_RAT') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_MOUSE = mapIds(org.Mm.eg.db, ENSEMBL_MOUSE, 'SYMBOL', 'ENSEMBL') %>% as.character())
                
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
                # ortho_df <- ortho_df %>%
                #         select(-CONFIDENCE)
                # Assign the HUGO symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_MOUSE') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, 
                                                   ENSEMBL_RAT, 
                                                   'SYMBOL', 'ENSEMBL') %>% as.character())
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

#' #### Retrieve Circadian Genes Associated with Tissue
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script:

# Circadian Genes 
#in_file <- paste0(WD,'/data/20200503_rnaseq-circadian-',TIS,'-mouse-rat-ortho_steep-yan.txt')
#circ_df <- read.table(file=in_file, sep = '\t', header = TRUE)

#' ## Place Genes in Genomic Ranges
#' #### Reference Genome and Annotation: Rnor_6.0 (GCA_000001895.4) assembly from Ensembl database (Release 96)
#' Found at: http://uswest.ensembl.org/Rattus_norvegicus/Info/Index.
#' 
#' FASTA: Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
#' 
#' GTF: Rattus_norvegicus.Rnor_6.0.96.gtf.gz ftp://ftp.ensembl.org/pub/release-96/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.96.gtf.gz

#+ Annotate Genes by Chromosome

################################################################################
#####     Annotate Genes by Chromosome       ###################################
################################################################################

### Determine which control samples are male and female
# Get the list of genes on the W chromosome

# Construct your own personal galgal5 reference genome annotation
# Construct from gtf file from Ensembl (same file used in mapping)
#ens_gtf <- paste0(WD,'/data/Rattus_norvegicus.Rnor_6.0.96.gtf')
#Rn_TxDb <- makeTxDbFromGFF(ens_gtf,
# format=c('gtf'),
# dataSource='Ensembl_Rattus6_gtf',
# organism='Rattus norvegicus',
# taxonomyId=NA,
# circ_seqs=DEFAULT_CIRC_SEQS,
# chrominfo=NULL,
# miRBaseBuild=NA,
# metadata=NULL)
# Save the Rat Genomic Ranges Object
#gf_file <- paste0(WD,'/data/20200603_Rnor-6.0.96-GRanges_steep.sqlite')
#saveDb(Rn_TxDb, file=gf_file)

# To load the annotation
gf_file <- paste0(WD,'/data/20200603_Rnor-6.0.96-GRanges_steep.sqlite')
Rn_TxDb <- loadDb(gf_file)
# Define Female specific sex genes (X chromosome)
# To examine chromosome names
seqlevels(Rn_TxDb)[1:23]
# Extract genes as GRanges object, then names
X_genes_gr <- genes(Rn_TxDb, columns = 'TXCHROM', filter = list(tx_chrom=c('X')))
# Collect ensembl gene ids for female specific genes
X_ens_id <- names(X_genes_gr)
# Examine the gene symbols
X_sym <- mapIds(org.Rn.eg.db, names(X_genes_gr), 'SYMBOL', 'ENSEMBL')
# Extract genes as GRanges object, then names
Y_genes_gr <- genes(Rn_TxDb, columns = 'TXCHROM', filter = list(tx_chrom=c('Y')))
# Collect ensembl gene ids for female specific genes
Y_ens_id <- names(Y_genes_gr)
sex_ens_id <- c(X_ens_id,Y_ens_id)
# Examine the gene symbols
Y_sym <- mapIds(org.Rn.eg.db, names(Y_genes_gr), 'SYMBOL', 'ENSEMBL')

#' ## Combine Annotations Tables

#+ Combine Annotations Tables
################################################################################
#####     Combine Annotations Tables       ##################################
################################################################################
# Load the R2 values from modeling
models_file <- paste0(WD,'/data/20200603_rnaseq-tissue-models-residuals-r2-table_steep.txt')
models_df <- read.table(file = models_file ,sep = '\t', header = T, check.names = F) %>%
        as_tibble()

# Load the DE values from different DE comparisons
de_file <- paste0(WD,'/data/20200624_rnaseq-tissue-all-de-table_steep.txt')
de_df <- read.table(file = de_file ,sep = '\t', header = T, check.names = F) %>%
        as_tibble()

# Load the Manova file
manova_file <- paste0(WD,'/data/20200413_rnaseq-tissue-manova-table_steep.txt')
manova_df <- read.table(file = manova_file ,sep = '\t', header = T, check.names = F) %>%
        as_tibble()

# Load the SOM and K-means cluster file
clusters_file <- paste0(WD,'/data/20200413_rnaseq-tissue-cluster-table_steep.txt')
clusters_df <- read.table(file = clusters_file ,sep = '\t', header = T, check.names = F) %>%
        as_tibble()

# Join the tables
bin_df <- full_join(models_df, de_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
        full_join(manova_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
        full_join(clusters_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
        mutate(BIN_TYPE = case_when(SIN1_R2 > GAM1_R2 + 0.15 ~ 'Circadian',
                                    GAM1_R2 > SIN1_R2 + 0.15 ~ 'Exercise',
                                    ((SIN1_R2 <= GAM1_R2 + 0.15) & (GAM1_R2 <= SIN1_R2 + 0.15)) ~ 'Ambiguous'))

#' ## Additional Annotations
#' 
#+ Additional Annotations
################################################################################
#####     Additional Annotations       #########################################
################################################################################

# Adds R2 Distance and Rat SYmbol
bin_df <- bin_df %>%
        mutate(R2_DIST = abs(SIN1_R2 - GAM1_R2)) %>%
        mutate(ENSEMBL_RAT = as.character(ENSEMBL_RAT)) %>%
        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))


# Find mouse orthologs for enrichr
orth_df <- bin_df %>%
        select(ENSEMBL_RAT) %>% 
        unique() %>% as.data.frame() %>%
        rat_mouse_ortho(column = 'ENSEMBL_RAT', direction = 'rat2mouse') %>%
        as_tibble()
bin_df <- left_join(bin_df, orth_df, by='ENSEMBL_RAT')

#' ## Enrichr: EDA
#' 
#+ Enrichr: EDA
################################################################################
#####     Enrichr: EDA       ###################################################
################################################################################

# Analysis performed per Tissue
enrichr_final_df <- data.frame()
for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
print(TISSUE)
for(gene_group in c('Circadian', 'Ambiguous','Exercise')){
print(gene_group)
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

# Select the "Circadian Modeled" Genes
df <- bin_df %>%
        filter(TISSUE == TISSUE1) %>%
        filter(BIN_TYPE == gene_group) %>%
        select(TISSUE, ENSEMBL_RAT, CONFIDENCE, ENSEMBL_MOUSE, 
               SYMBOL_MOUSE, SIN1_R2, GAM1_R2, MANOVA_PVAL, BIN_TYPE,R2_DIST) %>%
        unique() %>%
        arrange(desc(R2_DIST))

# Pathway Analysis
###############################

# Collect gene symbols
symbols <- df %>%
        filter(!is.na(SYMBOL_MOUSE)) %>%
        select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()

# Collect the desired pathways to query for enrichment
desired <- c('KEGG_2019_Mouse')
enriched <- enrichr(symbols, desired)

n<-1
enrichr_df <- data.frame()
for(db in desired){
        # Collapse the data into flat format
        if(!is.null(nrow(enriched[[desired]]))){
                enriched[[n]]$Data.Base <- db
                enrichr_df <- rbind(enrichr_df,enriched[[n]])
        }
        n <- n + 1
}
# Generate a new column for the number of genes per enriched term
if(nrow(enrichr_df) > 0){
        enrichr_df  <- enrichr_df %>% separate(Overlap, sep = '/', into = c('Overlap.N','Total.N'))
        enrichr_df$FDR_nlog <- -log(enrichr_df$Adjusted.P.value)
        # Filter the results
        enrichr_plot <- enrichr_df %>%
                filter(Total.N <= 500) %>%
                filter(Total.N >= 10) %>%
                filter(Adjusted.P.value <= 0.1) %>%
                arrange(desc(FDR_nlog)) %>%
                dplyr::slice(1:40) %>%
                group_by(Term) %>%
                arrange(FDR_nlog) %>%
                ungroup() %>%
                mutate(TISSUE = TISSUE1) %>%
                mutate(BIN_TYPE = gene_group)
        
        enrichr_final_df <- rbind(enrichr_final_df, enrichr_plot)
        # Create a special column that will allow for proper ordering of Pathways in Plot:
        # Ordering will have 2 priorities in this order: Up/Down regulation & -Log(FDR)
        # Visualize results with barplot
        #pdf(paste0(WD,'/plots/20200426_rnaseq-',TIS,'-pathways-sexmod-C0vsC7_steep.pdf'),
        #    width = 12, height = 6)
        p <- ggplot(enrichr_plot, 
                    aes(x=reorder(Term, FDR_nlog), y=FDR_nlog)) +
                geom_bar(stat='identity') +
                coord_flip() +
                ylab('-Log(FDR)') +
                xlab('') +
                theme_linedraw() +
                guides(fill=guide_legend(title='Gene Expression')) +
                theme(panel.background = element_blank(),
                      plot.background = element_blank(),
                      strip.background = element_blank(),
                      axis.text = element_text(size=12, colour = 'black'),
                      axis.ticks = element_line(colour = 'black'),
                      axis.title=element_text(size=12,face='bold'),
                      strip.text = element_blank()) +
                theme(panel.grid.major = element_blank()) +
                ggtitle(paste0(TISSUE1,' ',gene_group,' 
Pathway Enrichment (',desired,')'))
        pdf(paste0(WD,'/plots/20200705_rnaseq-',TISSUE1,'-',gene_group,'-pathways_steep.pdf'),
            width = 12, height = 6)
        plot(p)
        dev.off()
}
}
}

enrichr_final_df$TISSUE <- as.factor(enrichr_final_df$TISSUE)

# Circadian Plot
#########################################

enrichr_circ_df <- enrichr_final_df %>%
        filter(Total.N <= 500) %>%
        filter(Total.N >= 10) %>%
        filter(Adjusted.P.value <= 0.1) %>%
        arrange(desc(FDR_nlog)) %>%
        filter(BIN_TYPE == 'Circadian') %>%
        group_by(TISSUE) %>%
        dplyr::slice(1:10)

enrichr_circ_df %>%
        filter(BIN_TYPE == 'Circadian') %>%
        ggplot(aes(x=reorder(Term, FDR_nlog), y=FDR_nlog)) +
        geom_bar(stat='identity') +
        coord_flip() +
        ylab('-Log(FDR)') +
        xlab('') +
        theme_linedraw() +
        guides(fill=guide_legend(title='Gene Expression')) +
        # theme(panel.background = element_blank(),
        #        plot.background = element_blank()) +
        # #       strip.background = element_blank(),
        # #       axis.text = element_text(size=12, colour = 'black'),
        # #       axis.ticks = element_line(colour = 'black'),
        # #       axis.title=element_text(size=12,face='bold'),
        # #       strip.text = element_blank()) +
        theme(panel.grid.major = element_blank()) +
        ggtitle(paste0('Circadian-Modeled Genes
Pathway Enrichment (',desired,')')) +
        facet_wrap(vars(TISSUE))

