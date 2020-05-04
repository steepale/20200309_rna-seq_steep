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

#+ Circadian Gene Sets
################################################################################
#####     Circadian Gene Sets       ############################################
################################################################################

# Circadian Genes
in_file <- paste0(WD,'/data/20080516_mouse-tissue-circadian-genes_yan.txt')
circ_df <- read.table(in_file,sep = '\t', header = TRUE,check.names = FALSE)
# Adjust gene symbol names
names(circ_df)[1] <- "SYMBOL_MOUSE"

# Collect the tissues in a vector
yan_tissues <- names(circ_df)[names(circ_df) %!in% c('SYMBOL_MOUSE','NUM.TISSUE','RANGE.P','PEAK.MEAN') ]
#yan_tissues <- yan_tissues[9:14]
#tis <- 'SCN'
circ_tis <- list()
circ_tis2 <- list()
for(tis in (yan_tissues)){
        print(tis)
        circ_tis[[tis]] <- circ_df %>%
                select("SYMBOL_MOUSE", 
                       all_of(tis), 
                       "NUM.TISSUE", "RANGE.P", "PEAK.MEAN") %>% 
                filter(!is.na(!!rlang::sym(tis))) %>%
                filter(NUM.TISSUE >= 1)
        input_n <- circ_tis$SYMBOL_MOUSE %>% unique() %>% length() %>% as.character()
        # Some genes were annotated with synonyms by Yan et al. We manually convert these synonyms to gene symbols
        circ_tis[[tis]]$SYMBOL_MOUSE <- as.character(circ_tis[[tis]]$SYMBOL_MOUSE)
        circ_tis[[tis]]$ENSEMBL_MOUSE <- mapIds(org.Mm.eg.db, 
                                                circ_tis[[tis]]$SYMBOL_MOUSE, 
                                                "ENSEMBL", "SYMBOL")
        # Useful trick, search symbol in aliases and grab symbol
        circ_tis[[tis]]$ALIAS_MOUSE <- mapIds(org.Mm.eg.db, 
                               circ_tis[[tis]]$SYMBOL_MOUSE, "SYMBOL", "ALIAS",
                               multiVals = "first")
        # With a slow for loop, iterate through rows of dataframe and adjust nomenclature scheme as needed
        for( r in 1:nrow(circ_tis[[tis]])){
                ENSEMBL <- circ_tis[[tis]][r,'ENSEMBL_MOUSE']
                ALIAS <- circ_tis[[tis]][r,'ALIAS_MOUSE']
                if(is.na(ENSEMBL)){
                        if(!is.na(ALIAS)){
                                # Remember that ALIAS is actually from Ensembl's SYMBOL column
                                ENSEMBL_RPL <- mapIds(org.Mm.eg.db, ALIAS, "ENSEMBL", "SYMBOL")
                                circ_tis[[tis]][r,'ENSEMBL_MOUSE'] <- ENSEMBL_RPL
                                circ_tis[[tis]][r,'SYMBOL_MOUSE'] <- ALIAS
                        }
                }
        }
        circ_tis[[tis]] <- circ_tis[[tis]] %>%
                filter(!is.na(ENSEMBL_MOUSE)) %>%
                select(-ALIAS_MOUSE)
        tryCatch({
        circ_tis2[[tis]] <- mouse2rat_ortho(circ_tis[[tis]])
        output_n <- circ_tis2[[tis]]$ENSEMBL_RAT %>% 
                unique() %>% 
                length() %>% as.character()
        #' ##### High confidence ortholgos were collected from Yan et. al. (Supplementary Table 2) and converted to high confidence rat orthologs with Ensembl.
        #' ###### Mouse annotation: GRCm38.p6
        #' ###### Rat annotation: Rnor_6.0
        #' ###### Stats:
        #' * Mouse genes input: `r input_n`
        #' * High confidence Rat orthologs output: `r output_n`
        #' ###### Steps in ortholog selection:
        #' * Demonstration of circadian oscilliations in mouse kidney (from Yan et. al.)
        #' * Demonstration of circadian oscilliations in 1 (min Kidney) or more tissues (from Yan et. al.)
        #' * Mouse genes were removed that did not have Ensembl gene symbol (or alias) to gene id 
        #'     * First choice of 1:many were included
        #' * Genes were removed if they did not have a high orthology confidence score between mouse gene id and rat gene id  (binary value 0|1)
        #' * Duplicate orthologs were removed
        #'     * First choice of 1:many were included
        
        # Cleanup
        circ_tis[[tis]] <- circ_tis2[[tis]]
        #rm(circ_tis2[[tis]])
        # Collect kidney circadian kidney genes
        circ_tis_ens <- circ_tis[[tis]] %>%
                select(ENSEMBL_RAT) %>% unlist()
        #' #### Data Save: Mouse and Rat Orthologus Circadian Genes
        
        #+ Data Save: Mouse and Rat Orthologus Circadian Genes 
        ########################################################################
        ######## Data Save: Mouse and Rat Orthologus Circadian Genes  ##########
        ########################################################################
        
        # Circadian Genes (Kidney)
        out_file <- paste0(WD,'/data/20200503_rnaseq-circadian-',tis,'-mouse-rat-ortho_steep-yan.txt')
        write.table(circ_tis[[tis]], file=out_file, row.names = FALSE, quote = FALSE, sep = '\t')})
}


## Make one custom adjustment (Biomart did not catch this annotation)
#circ_tis <- circ_tis %>%
#        mutate(SYMBOL_RAT = ifelse(
#                ENSEMBL_RAT %in% c('ENSRNOG00000060956','ENSRNOG00000047309',
#                                   'ENSRNOG00000000419','ENSRNOG00000018536'), 
#                SYMBOL_MOUSE, SYMBOL_RAT))




