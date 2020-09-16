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
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr','rlang')
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

#' ### Differential Gene Expression

#+ Differential Gene Expression
################################################################################
###########Differential Gene Expression  #######################################
################################################################################
# Generate an empty dataframe for final output
tissue_de_df <- data.frame()
# TISSUE <- 'Hypothalamus'
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
        
        de_df <- data.frame()
        #G1G2 <- 'C0C7'
# for(G1G2 in c('C0C7','C0E0','C0E0.5','C0E1','C0E4','C0E7','C0E24','C0E48',
#               'C7E0','C7E0.5','C7E1','C7E4','C7E7','C7E24','C7E48')){
for(G1G2 in c('C0C7','C7E7','C0E0')){
        if(G1G2 == 'C0C7'){
                cohorts <- c('Control - IPE', 'Control - 7 hr')
        }else if(G1G2 == 'C0E0'){
                cohorts <- c('Control - IPE', 'Exercise - IPE')
        }else if(G1G2 == 'C0E0.5'){
                cohorts <- c('Control - IPE', 'Exercise - 0.5 hr')
        }else if(G1G2 == 'C0E1'){
                cohorts <- c('Control - IPE', 'Exercise - 1 hr')
        }else if(G1G2 == 'C0E4'){
                cohorts <- c('Control - IPE', 'Exercise - 4 hr')
        }else if(G1G2 == 'C0E7'){
                cohorts <- c('Control - IPE', 'Exercise - 7 hr')
        }else if(G1G2 == 'C0E24'){
                cohorts <- c('Control - IPE', 'Exercise - 24 hr')
        }else if(G1G2 == 'C0E48'){
                cohorts <- c('Control - IPE', 'Exercise - 48 hr')
        }else if(G1G2 == 'C7E0'){
                cohorts <- c('Control - 7 hr', 'Exercise - IPE')
        }else if(G1G2 == 'C7E0.5'){
                cohorts <- c('Control - 7 hr', 'Exercise - 0.5 hr')
        }else if(G1G2 == 'C7E1'){
                cohorts <- c('Control - 7 hr', 'Exercise - 1 hr')
        }else if(G1G2 == 'C7E4'){
                cohorts <- c('Control - 7 hr', 'Exercise - 4 hr')
        }else if(G1G2 == 'C7E7'){
                cohorts <- c('Control - 7 hr', 'Exercise - 7 hr')
        }else if(G1G2 == 'C7E24'){
                cohorts <- c('Control - 7 hr', 'Exercise - 24 hr')
        }else if(G1G2 == 'C7E48'){
                cohorts <- c('Control - 7 hr', 'Exercise - 48 hr')
        }
        print(G1G2)
        if(TISSUE == c('Gastrocnemius')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Gastrocnemius') %>%
                        filter(Seq_batch == 'MSSM_1') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS) %>%
                        filter(animal.key.anirandgroup %in% cohorts)
        }else if(TISSUE == c('Lung')){
                tod_cols <- col_data %>%
                        filter(Tissue == 'Lung') %>%
                        filter(Seq_batch == 'MSSM_3') %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS) %>%
                        filter(animal.key.anirandgroup %in% cohorts)
        }else{
                # Filter Samples (meta)
                tod_cols <- col_data %>%
                        filter(Tissue == TISSUE) %>%
                        filter(!is.na(animal.registration.sex)) %>%
                        filter(sample_key %!in% OUTLIERS) %>%
                        filter(animal.key.anirandgroup %in% cohorts)
        }
        rownames(tod_cols) <- tod_cols$sample_key
        # Collect tissue specific counts
        tod_counts <- count_data[,as.character(tod_cols$sample_key)]
        
        # Generate an annotation specific to this analysis
        # Make sure the "control" or "untreated" level is first
        
        tod_cols <- tod_cols %>%
                mutate(binary = factor(
                        case_when(animal.key.anirandgroup == cohorts[1] ~ 
                                          str_replace(cohorts[1],' - ','_') %>%
                                          str_replace(' ',''),
                                  animal.key.anirandgroup == cohorts[2] ~ 
                                          str_replace(cohorts[2],' - ','_') %>% 
                                          str_replace(' ','')),
                                                levels = c(str_replace(cohorts[1],' - ','_') %>% str_replace(' ',''),
                                                           str_replace(cohorts[2],' - ','_') %>% str_replace(' ',''))))
        
        row.names(tod_cols) <- as.character(tod_cols$sample_key)
        #' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
        print(all(rownames(tod_cols) == colnames(tod_counts)))
        
        # Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
        if(length(ADJ_VAR) == 2){
                design = as.formula(paste0('~',ADJ_VAR[1],' + ',ADJ_VAR[2],' + binary'))
        }else if(ADJ_VAR == 'None'){
                design = as.formula('~binary')
        }else{
                design = as.formula(paste0('~',ADJ_VAR[1],' + binary'))
        }
        
        #design <- formula(paste0('~animal.registration.sex + animal_time_last_fed + binary'))
        (title = paste0('Design: ',as.character(design)))
        dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                                       colData = tod_cols,
                                       design = design)
        
        zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
        reads_n <- 1
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
        # rld <- DESeq2::vst(dds, blind = F)
        # pc_cor_df <- cor_PC_1_4(rld = rld, ntop = 20000, intgroups = names(tod_cols)) %>%
        #         #filter(Adjusted_R_Sq > 0.2) %>%
        #         arrange(PC)
        # # Examine Variables of interest for Pcs 1 - 4
        # pc_cor_df %>%
        #         filter(PC %in% c(1))
        # # animal.key.anirandgroup, animal.registration.sex, animal_time_last_fed,
        # # calculated.variables.deathtime_after_acute, RIN
        # v <- 'animal.registration.d_birth'
        # p <- DESeq2::plotPCA(rld, intgroup =v, ntop = 20000) +
        #         ggtitle(paste0(TISSUE,': ',v)) +
        #         coord_equal()
        # plot(p)
        # # The batch effect can only be removed with limma
        # # https://support.bioconductor.org/p/76099/ (See Michael Love's Comment)
        # assay(rld) <- limma::removeBatchEffect(assay(rld), rld[[pri_var]])
        # Generate a DESeq2 Model
        dds <- DESeq(dds)
        # Generate a results table and sort by fdr
        res <- results(dds, alpha = 0.05,lfcThreshold=0.25)
        res <- res[order(res$padj),]
        # Generate a summary of the results
        summary(res)
        res_df <- as.data.frame(res)
        # DE dataframe
        G1G2_df <- data.frame("TISSUE" = TISSUE,
                              "ENSEMBL_RAT" = row.names(res_df))
        G1G2_df[[paste0(G1G2,'_log2FoldChange')]] <- res_df$log2FoldChange
        G1G2_df[[paste0(G1G2,'_padj')]] <- res_df$padj
        
        #' ### Join Results by Tissue and Gene
        
        #+ Join Results by Tissue and Gene
        ################################################################################
        ########### Join Results by Tissue and Gene  ###################################
        ################################################################################
        
        # Tissue specific
        if(!nrow(de_df) > 0){
                de_df <- rbind(de_df, G1G2_df)
        }else{
                de_df <- full_join(de_df,G1G2_df, by = c('TISSUE','ENSEMBL_RAT'))
        }
}
        # Concatenate the dataframes
        tissue_de_df <- rbind(tissue_de_df , de_df)
        # Save the final output table
        table_file <- paste0(WD,'/data/20200624_rnaseq-tissue-all-de-table_steep.txt')
        write.table(tissue_de_df, file = table_file,sep = '\t',row.names = F,quote = F)
}

#' ### Determine DE Genes

#+ Determine DE Genes
################################################################################
########### Determine DE Genes  ###################################
################################################################################

# Load the DE values from different DE comparisons
de_file <- paste0(WD,'/data/20200624_rnaseq-tissue-all-de-table_steep.txt')
de_df <- read.table(file = de_file ,sep = '\t', header = T, check.names = F) %>%
        as_tibble()

# Iterate through each possibility
for(G1G2 in c('C0C7','C0E0','C0E0.5','C0E1','C0E4','C0E7','C0E24','C0E48',
              'C7E0','C7E0.5','C7E1','C7E4','C7E7','C7E24','C7E48')){
        lfc <- paste0(G1G2, '_log2FoldChange')
        padj <- paste0(G1G2, '_padj')
        sig <- paste0(G1G2, '_SIG')
        de_df <- de_df %>%
                mutate(!!sym(sig) := ifelse(((!!sym(lfc) >= 0.25 | !!sym(lfc) <= -0.25) 
                                            & !!sym(padj) <= 0.05),'SIG','NON-SIG'))
}
table_file <- paste0(WD,'/data/20200624_rnaseq-tissue-all-de-table_steep.txt')
write.table(de_df, file = table_file,sep = '\t',row.names = F,quote = F)

