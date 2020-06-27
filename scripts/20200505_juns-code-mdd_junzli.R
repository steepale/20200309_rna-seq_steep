#'---
#' title: "Jun's MDD ANalysis"
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
#' * Recreate Jun's ANalysis
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
#' * 

#+ Load the Data
################################################################################
#####     Load & Clean Data      ###############################################
################################################################################


# Count matrix
# 1. read in data and gender check
################################################################################
in_file <- paste0(WD,'/data/20200505_dlpfc-avg_junzli.txt')
dlpfc<-as.matrix(read.delim(in_file,header=T,row.names=1,sep="\t", check.names = F))
colnames(dlpfc) %>% head()
in_file <- paste0(WD, '/data/20200505_dlpfc-avg-sample-subjects_junzli.txt')
attributes <- read.table(in_file, sep = "\t", header = TRUE, fill = TRUE)
row.names(attributes) <- attributes$SUBJECT

dlpfc.sample <- attributes[colnames(dlpfc), c("SUBJECT", "COHORT", "SITE","DIAGNOSIS", "DIAGNOSIS_SSRI", "AGONAL_FACTOR", "PH", "GENDER","RACE", "AGE", "ACI", "SAMPLE_TYPE", "ORIGIN", "CHIP_TYPE","BLOCK","BLOCK_FULLNAME")]

# Ordering and normalization
################################################################################
dlpfc<-dlpfc[,order(dlpfc.sample$SITE,dlpfc.sample$COHORT,dlpfc.sample$ORIGIN, dlpfc.sample$SUBJECT)]
dlpfc.sample<-dlpfc.sample[order(dlpfc.sample$SITE,dlpfc.sample$COHORT,dlpfc.sample$ORIGIN, dlpfc.sample$SUBJECT),]

library(limma)
dlpfc.norm <- normalizeQuantiles(dlpfc)
dlpfc.cor<-cor(dlpfc.norm)
image(dlpfc.cor,col=redblue(100))

range<-max(dlpfc.cor)-0.8
tmp<-unique(as.factor(dlpfc.sample$COHORT))
cohort<-matrix(NA,dim(dlpfc.sample)[1],5)
i <- 1
for (i in 1:length(tmp)){
        cohort[ dlpfc.sample$COHORT==tmp[i],1]<-0.8+range*(i-1)/(length(tmp)-1)
}
cohort[,2]<- cohort[,3]<- cohort[,4]<- cohort[,5]<- cohort[,1]

tmp<-c(unique(as.factor(dlpfc.sample$SITE)), 'UCDavis')
site<-matrix(NA,dim(dlpfc.sample)[1],5)
for (i in 1:length(tmp)){
        site[ dlpfc.sample$SITE==tmp[i],1]<-0.8+range*(i-1)/(length(tmp)-1)
}
site[,2]<- site[,3]<- site[,4]<- site[,5]<- site[,1]

tmp<-unique(as.factor(dlpfc.sample$ORIGIN))
origin<-matrix(NA,dim(dlpfc.sample)[1],5)
for (i in 1:length(tmp)){
        origin[ dlpfc.sample$ORIGIN==tmp[i],1]<-0.8+range*(i-1)/(length(tmp)-1)
}
origin[,2]<- origin[,3]<- origin[,4]<- origin[,5]<- origin[,1]

image(cbind(site,cohort,origin,dlpfc.cor),col=redblue(100),zlim=c(0.8,1))

plot(apply(dlpfc.cor,1,median))

pca.dlpfc<-prcomp(t(dlpfc.norm))
plot(pca.dlpfc$x[,1],pca.dlpfc$x[,2])
points(pca.dlpfc$x[apply(dlpfc.cor,1,median)<0.90,1],pca.dlpfc$x[apply(dlpfc.cor,1,median)<0.90,2],pch=19,col=2)
points(pca.dlpfc$x[apply(dlpfc.cor,1,median)<0.85,1],pca.dlpfc$x[apply(dlpfc.cor,1,median)<0.85,2],pch=19,col=3)






