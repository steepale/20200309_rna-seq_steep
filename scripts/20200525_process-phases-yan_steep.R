#'---
#' title: "Data Processing: Process Phase Values of Circadian Genes from "
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200525"
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
#' * Convert genes to high confidence rat orthologs
#'     * Low confidence orthologs were included in this analysis
#' * Save orthologs and peak values
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
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","org.Hs.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","vcd")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
counts <- DESeq2::counts

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
direction <- "mouse2rat"
column <- "ENSEMBL_MOUSE"
x <- data.frame(x)
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
                #ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
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
                                                   "SYMBOL", "ENSEMBL"))
        }
        # Return the output
        x
}
#########################################################################################

#' ## Load & Clean Data
#' ##### Data files to load:
#' * Phase Shift Information from TODO: Determine study

#+ Load the Data
################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Phases
in_file <- paste0(WD,'/data/20200525_Yan-pcbi-1000193-s005-trim_junzli.txt')
phase_df <- read.table(in_file,sep = '\t', header = TRUE,check.names = FALSE) %>% as_tibble()
names(phase_df)[1] <- "SYMBOL_MOUSE"

# Let's take a look at those phases
phase_df %>% 
        filter(RANGE.P <= 0.05) %>%
        ggplot(aes(x=PEAK.MEAN)) +
        geom_density(alpha=0.5) +
        #geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
        labs(title="Density of Phases: \nPhases of Circadian Rhythms",x="Phases", y = "Frequency") +
        #xlim(0,24) +
        theme_classic()

# Collect ensembl gene ids for mouse and convert to rat 
phase_df$SYMBOL_MOUSE <- as.character(phase_df$SYMBOL_MOUSE)
# (gene id 2 gene id conversion is adventageous over symbol to symbol)
phase_df$ENSEMBL_MOUSE <- mapIds(org.Mm.eg.db, phase_df$SYMBOL_MOUSE, 
                                 "ENSEMBL", "SYMBOL", 
                                 multiVals = "first") %>% as.character()

# Collect the orthologs
phase_df <- rat_mouse_ortho(as.data.frame(phase_df), 
                             column = "ENSEMBL_MOUSE", 
                             direction = 'mouse2rat') %>% as_tibble()

# Ensure that the distribution remains the same
phase_df %>% 
        #filter(is.na(SYMBOL_RAT)) %>% 
        filter(RANGE.P <= 0.05) %>% 
        ggplot(aes(x=PEAK.MEAN)) +
        geom_density(alpha=0.5) +
        #geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
        labs(title="Density of Phases: \nPhases of Circadian Rhythms",x="Phases", y = "Frequency") +
        #xlim(0,24) +
        theme_classic()

# Save the phase values
out_file <- paste0(WD,'/data/20200525_phase-yan_steep.txt')
write.table(phase_df, file = out_file, quote = FALSE, 
            row.names = FALSE, col.names = TRUE, sep ='\t')




