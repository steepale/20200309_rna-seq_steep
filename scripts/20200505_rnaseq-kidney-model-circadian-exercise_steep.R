#'---
#' title: "PASS1A Rat All Tissue: -- Modeling of RNASeq Data Expression"
#' author: "Alec Steep & Jiayu Zhang (Code copied and adapted from Jun Z. Li)" 
#' date: "20200505"
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
#' * Model Circadian Rhythms (5 major time points)
#'     * SIN/COS Linear Model
#'     * y = B_1SIN(TOD) + B_2COS(TOD) + B_0
#' * Model Exercise Effects (7 time points)
#'     * Comparison of linear models:
#'     * Cubic: y = ax^3 + bx^2 + cx + d
#'     * y = ax^4 + bx^3 + cx^2 + dx + e
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
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","tiff")
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

# Sine and cosine functions
################################################################################
SIN <- function(t) sin(2*pi*t/24)
COS <- function(t) cos(2*pi*t/24)
################################################################################

# I’m writing a custom function, it returns PVE (explained below).  We need to 
# tweak the codes to return the beta’s. 
################################################################################
#expr <- tod_counts[i,]
fit.FR <- function(expr,TOD) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TOD)),SIN(TOD),COS(TOD)))
        #this is the code for the formula above, where “solve” is the matrix inversion operation - a nice trick to have.
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr) 
        SSz <- t(beta) %*% t(M) %*% expr - (sum(expr)^2) / length(expr)
        SSy <- sum(expr^2)-(sum(expr)^2)/ length(expr)
        #The SSz is the variance explained by the fitted function, SSy is the total variance in the data (in this case, the variance of the 55 values for this gene). PVE is Percent Variance Explained.  It measures the goodness of fit.
        PVE <- SSz / SSy 
        PVE
}
################################################################################

# I’m writing a custom function, it returns PVE (explained below).  We need to 
# tweak the codes to return the beta’s. 
################################################################################
#expr <- tod_counts[i,]
fit.3_pve <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3)) 
        #this is the code for the formula above, where “solve” is the matrix 
        #inversion operation - a nice trick to have.
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr) 
        SSz <- t(beta) %*% t(M) %*% expr - (sum(expr)^2) / length(expr)
        SSy <- sum(expr^2)-(sum(expr)^2)/ length(expr)
        #The SSz is the variance explained by the fitted function, SSy is the total variance in the data (in this case, the variance of the 55 values for this gene). PVE is Percent Variance Explained.  It measures the goodness of fit.
        PVE <- SSz / SSy
        BV <- list(beta[1], beta[2], beta[3], beta[4], PVE)
        PVE
}
################################################################################
################################################################################
#expr <- tod_counts[i,]
fit.4_pve <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        #this is the code for the formula above, where “solve” is the matrix 
        #inversion operation - a nice trick to have.
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr) 
        SSz <- t(beta) %*% t(M) %*% expr - (sum(expr)^2) / length(expr)
        SSy <- sum(expr^2)-(sum(expr)^2)/ length(expr)
        #The SSz is the variance explained by the fitted function, SSy is the total variance in the data (in this case, the variance of the 55 values for this gene). PVE is Percent Variance Explained.  It measures the goodness of fit.
        PVE <- SSz / SSy
        BV <- list(beta[1], beta[2], beta[3], beta[4], PVE)
        PVE
}
fit.4_B1 <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr)
        beta[1]
}
fit.4_B2 <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr)
        beta[2]
}
fit.4_B3 <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr)
        beta[3]
}
fit.4_B4 <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr)
        beta[4]
}
fit.4_B5 <- function(expr,TPE) {
        expr <- as.vector(expr) %>% as.numeric()
        M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
        beta <- solve(t(M)%*%M)%*%(t(M)%*%expr)
        beta[5]
}
################################################################################

# Function to calculate the proportion of variance explained
################################################################################
pve.3<-function(data,TPE) {
        tmp<- apply(data,1,function(x) fit.3_pve(x,TPE))
        return(tmp)
}
################################################################################
################################################################################
pve.4<-function(data,TPE) {
        tmp<- apply(data,1,function(x) fit.4_pve(x,TPE))
        return(tmp)
}
################################################################################
# Function to calculate the proportion of variance explained
################################################################################
pve.FR<-function(data,TOD) {
        tmp<- apply(data,1,function(x) fit.FR(x,TOD))
        return(tmp)
}
################################################################################


#Some codes for producing an empirical p value for the PVE by permuting TOD many times, say 100 or 1000 times, and compare the observed PVE with the null distribution of simulated PVE’s
# A permutation test that tests the null distribution of PVEs
#data <- tod_counts
#iter = 100
################################################################################
p.test.cos<- function(data,TOD,iter,every=10) {
        pval <- rep(0,dim(data)[1])
        cat(iter, "permutations in progress\n")
        PVE<-pve.FR(data,TOD)
        for(i in 1:iter) {
                # Shuffles the samples
                new.data<-data[,sample(1:dim(data)[2])]
                # Calculate the proportion of variance explained
                tmp.PVE <- pve.FR(new.data,TOD)
                # Generates a running count of simulated PVE 
                # that are larger or equal to the actual PVE
                pval <- pval + as.numeric(tmp.PVE >= PVE)
                if(i %% every == 0) cat(every)
                if(i %% (10*every) == 0) cat("\n")
        }
        # Calculate the p value by dividing the running count by the number of iterations
        pval <- pval / iter
}
################################################################################
################################################################################
p.test.3 <- function(data,TPE,iter,every=10) {
        pval <- rep(0,dim(data)[1])
        cat(iter, "permutations in progress\n")
        PVE<-pve.3(data,TPE)
        for(i in 1:iter) {
                # Shuffles the samples
                new.data<-data[,sample(1:dim(data)[2])]
                # Calculate the proportion of variance explained
                tmp.PVE <- pve.3(new.data,TPE)
                # Generates a running count of simulated PVE 
                # that are larger or equal to the actual PVE
                pval <- pval + as.numeric(tmp.PVE >= PVE)
                if(i %% every == 0) cat(every)
                if(i %% (10*every) == 0) cat("\n")
        }
        # Calculate the p value by dividing the running count by the number of iterations
        pval <- pval / iter
}
################################################################################
################################################################################
p.test.4 <- function(data,TPE,iter,every=10) {
        pval <- rep(0,dim(data)[1])
        cat(iter, "permutations in progress\n")
        PVE<-pve.4(data,TPE)
        for(i in 1:iter) {
                # Shuffles the samples
                new.data<-data[,sample(1:dim(data)[2])]
                # Calculate the proportion of variance explained
                tmp.PVE <- pve.4(new.data,TPE)
                # Generates a running count of simulated PVE 
                # that are larger or equal to the actual PVE
                pval <- pval + as.numeric(tmp.PVE >= PVE)
                if(i %% every == 0) cat(every)
                if(i %% (10*every) == 0) cat("\n")
        }
        # Calculate the p value by dividing the running count by the number of iterations
        pval <- pval / iter
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

#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script 20200413_rnaseq-kidney-female-temporal-clusters_steep.R

# Circadian Genes (Kidney)
in_file <- paste0(WD,'/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt')
circ_kid <- read.table(file=in_file, sep = '\t', header = TRUE)

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

#' ## Collect Samples of Interest and Normalize

#+ Collect Samples of Interest and Normalize
################################################################################
#####     Collect Samples of Interest and Normalize      #######################
################################################################################

# Filter Female Kidney Samples (meta)
tod_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(animal.registration.sex == 'Female') %>%
        filter(animal.key.anirandgroup %!in% c('Control - IPE', 'Control - 7 hr'))
rownames(tod_cols) <- tod_cols$sample_key

# Collect samples without NA values in TOD
nona_sams <- tod_cols %>%
        filter(!is.na(specimen.collection.t_death_hour)) %>%
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
#' To see the reads per million for each sample
sort(colSums(assay(dds)))/1e6

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)
#' Size facotrs are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
summary(sizeFactors(dds))

#rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for(n in 1){
        start_time <- Sys.time()
        rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

# Select the normailzed counts
tod_counts <- assay(rld) 

#' ## Model Circadian Rhythms by Time of Death (TOD)

#+ Model Circadian Rhythms by Time of Death (TOD)
################################################################################
#####     Model Circadian Rhythms by Time of Death (TOD)      ##################
################################################################################

# Time of death
TOD <- tod_cols %>%
        filter(sample_key %in% nona_sams) %>%
        select(specimen.collection.t_death_hour) %>% 
        unlist() %>% as.numeric()

# This sets up the “design matrix”, one row for each subject, with one column of 1’s, followed by two columns: a sine function and a cosine function. The goal is to estimate the three coefficients in y~B0+B1*sin(x)+B2*cos(x).
M <- as.matrix(cbind(rep(1,length(TOD)),SIN(TOD),COS(TOD))) 
# In case when we need to set up a cubic function: y~B0+B1*x+B2*x^2+B3*x^3.
plot(M[,2],M[,3],type="n",axes=F, xlab=" ",,ylab=" ",main="Time of Death")
text(M[,2],M[,3], labels = round(TOD),cex=0.8)

#The traditional approach is to run linear regression and loop through each of the p genes.
#The alternative approach is to solve for the three B’s (beta values) by matrix algebra where the vector beta can be calculated as B = (X^TX)^-1 * X^T * Y

#To get pve’s for all genes without looping p times:
tmp <- apply(tod_counts,1,function(x) fit.FR(x,TOD))

#If you want to verify it by looping (1ms per iter)
pve2 <- vector()
for (i in 1:nrow(tod_counts)) {
        pve2[i]<- fit.FR(tod_counts[i,],TOD)
}

# Calculate the p value empirically with permutation tests: shuffle time points randomly
set.seed(100)
tod_counts_pval<-p.test.cos(tod_counts,TOD,iter=100,every=10)

# Get the phase
t<-c(1:24)
sin_phase<- matrix(NA,length(TOD),24)
for (i in 1:24) {
        sin_phase[,i]<- sin(2*pi*(TOD+t[i])/24)
}

# Collect the maximum phase
phase_tod_counts <- matrix(NA,nrow(tod_counts),24)
for (i in 1:nrow(tod_counts)) {
        for (j in 1:24) {
                phase_tod_counts[i,j]<-cor(tod_counts[i,],sin_phase[,j])
        }
}
phase_tod_counts_max <-rep(NA,nrow(tod_counts))
for (i in 1:nrow(tod_counts)) {
        phase_tod_counts_max[i] <- c(1:24)[rank(phase_tod_counts[i,])==24]
}

# Store results in a dataframe
sinmod_df <- data.frame("ENSEMBL_RAT" = row.names(tod_counts),
                        "PVAL_SIN" = tod_counts_pval,
                        "PVE_SIN" = pve2,
                        "PHASE_SIN" = phase_tod_counts_max)
sinmod_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(sinmod_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")
# Filter data frame and arrange
sinmod_df <- sinmod_df %>%
        filter(PVAL_SIN <= 0.05) %>%
        arrange(desc(PVE_SIN))
dim(sinmod_df)
#' ## Model Exercise Response by Time Post Exercsie (TPE)

#+ Model Exercise Response by Time Post Exercsie
################################################################################
#####     Model Exercise Response by Time Post Exercsie      ###################
################################################################################

# Time post exercise
TPE <- tod_cols %>%
        filter(sample_key %in% nona_sams) %>%
        mutate(specimen.collection.t_exercise_hour = case_when(
                animal.key.anirandgroup == 'Exercise - IPE' ~ 0,
                animal.key.anirandgroup == 'Exercise - 0.5 hr' ~ 0.5,
                animal.key.anirandgroup == 'Exercise - 1 hr' ~ 1,
                animal.key.anirandgroup == 'Exercise - 4 hr' ~ 4,
                animal.key.anirandgroup == 'Exercise - 7 hr' ~ 7,
                animal.key.anirandgroup == 'Exercise - 24 hr' ~ 24,
                animal.key.anirandgroup == 'Exercise - 48 hr' ~ 0)) %>%
        select(specimen.collection.t_exercise_hour) %>% 
        unlist() %>% as.numeric()

# Cubic polynomial
################################################################################
# This sets up the “design matrix”, one row for each subject, with one column of 1’s, followed by two columns: a sine function and a cosine function. The goal is to estimate the four coefficients in y~B0+B1*x+B2*x^2+B3*x^3.
M_TPE <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3)) 
table(TPE)
# In case when we need to set up a cubic function: 
plot(M_TPE[,2],M_TPE[,3],type="n",axes=F, xlab=" ",,ylab=" ",main="Hours Post Exercise")
text(M_TPE[,2],M_TPE[,3], labels = round(TPE),cex=0.8)

#The traditional approach is to run linear regression and loop through each of the p genes.
#The alternative approach is to solve for the three B’s (beta values) by matrix algebra where the vector beta can be calculated as B = (X^TX)^-1 * X^T * Y

#To get pve’s for all genes without looping p times:
pve1 <- apply(tod_counts,1,function(x) fit.3_pve(x,TPE))

#If you want to verify it by looping (1ms per iter)
pve2 <- vector()
for (i in 1:nrow(tod_counts)) {
        pve2[i]<- fit.3_pve(tod_counts[i,],TPE)
}
all(pve1 == pve2)

# Calculate the p value empirically with permutation tests: shuffle time points randomly
set.seed(100)
tpe_counts_pval<-p.test.3(tod_counts,TPE,iter=100,every=10)

# Store results in a dataframe
cubmod_df <- data.frame("ENSEMBL_RAT" = row.names(tod_counts),
                        "PVAL_CUB" = tpe_counts_pval,
                        "PVE_CUB" = pve2)
cubmod_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(cubmod_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")
# Filter data frame and arrange
cubmod_df <- cubmod_df %>%
        filter(PVAL_CUB <= 0.05) %>%
        arrange(desc(PVE_CUB))

# 4 Degree Polynomial
################################################################################
# This sets up the “design matrix”, one row for each subject, with one column of 1’s, followed by two columns: a sine function and a cosine function. The goal is to estimate the four coefficients in y~B0+B1*x+B2*x^2+B3*x^3.
M <- as.matrix(cbind(rep(1,length(TPE)),TPE,TPE^2,TPE^3,TPE^4)) 
# In case when we need to set up a cubic function: 
#plot(M[,2],M[,3],type="n",axes=F, xlab=" ",,ylab=" ",main="Time of Death")
#text(M[,2],M[,3], labels = round(TPE),cex=0.8)

#The traditional approach is to run linear regression and loop through each of the p genes.
#The alternative approach is to solve for the three B’s (beta values) by matrix algebra where the vector beta can be calculated as B = (X^TX)^-1 * X^T * Y

#To get pve’s for all genes without looping p times:
pve1 <- apply(tod_counts,1,function(x) fit.4_pve(x,TPE))

#If you want to verify it by looping (1ms per iter)
pve2 <- vector()
for (i in 1:nrow(tod_counts)) {
        pve2[i]<- fit.4_pve(tod_counts[i,],TPE)
}
all(pve1 == pve2)
# Collect the betas
B1 <- apply(tod_counts,1,function(x) fit.4_B1(x,TPE))
B2 <- apply(tod_counts,1,function(x) fit.4_B2(x,TPE))
B3 <- apply(tod_counts,1,function(x) fit.4_B3(x,TPE))
B4 <- apply(tod_counts,1,function(x) fit.4_B4(x,TPE))
B5 <- apply(tod_counts,1,function(x) fit.4_B5(x,TPE))

# Calculate the p value empirically with permutation tests: shuffle time points randomly
set.seed(100)
tpe_counts_pval<-p.test.4(tod_counts,TPE,iter=100,every=10)
B1["ENSRNOG00000051980"]
B1[19088]
tod_counts[ensembl,]
tod_counts[19088,]
row_n <- which(row.names(tod_counts) == ensembl)
dim(tod_counts)
tod_counts[1:4,1:4]
length(B1)

# Store results in a dataframe
p4mod_df <- data.frame("ENSEMBL_RAT" = row.names(tod_counts),
                       "PVAL_P4" = tpe_counts_pval,
                       "PVE_P4" = pve2,
                       "B1_P4" = B1,
                       "B2_P4" = B2,
                       "B3_P4" = B3,
                       "B4_P4" = B4,
                       "B5_P4" = B5)
p4mod_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(p4mod_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")
# Filter data frame and arrange
p4mod_df <- p4mod_df %>%
        filter(PVAL_P4 <= 0.05) %>%
        arrange(desc(PVE_P4))
dim(p4mod_df)
#' ## Compare Models

#+ Compare Models
################################################################################
###################     Compare Models      ####################################
################################################################################

# Compare the polynomial models
######################################

# Select the shared genes between models
shared_mod_df <- p4mod_df %>%
        inner_join(cubmod_df, by = "ENSEMBL_RAT")

# Generate a scatterplot of the models' PVE's
shared_mod_df %>%
        ggplot(aes(x = PVE_P4, y= PVE_CUB)) +
        geom_point(alpha = 0.3) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("PVE: Polynomial (Degree 4)") +
        ylab("PVE: Polynomial (Degree 3)")

shared_mod_df %>%
        select(ENSEMBL_RAT, PVE_CUB, PVE_P4) %>%
        gather(MODEL, PVE, -ENSEMBL_RAT) %>%
        ggplot(aes(MODEL,PVE)) +
        geom_violin() +
        geom_boxplot(alpha = 0.1) +
        scale_x_discrete(labels=c("PVE_CUB" = "Polynomial (Degree 3)", 
                                  "PVE_P4" = "Polynomial (Degree 4)"))

shared_mod_df$PHASE_SIN

# Compare best polynomial model to SIN model
######################################
# Select the shared genes between models
shared_mod_df <- sinmod_df %>%
        inner_join(p4mod_df, by = "ENSEMBL_RAT")

# Generate a venn diagram of the models
# Venn DIagram
venn.plot <- venn.diagram(
        x = list(
                SIN = sinmod_df$ENSEMBL_RAT,
                POLY4 = p4mod_df$ENSEMBL_RAT),
        filename = paste0(WD,"/plots/20200505_rnaseq-kidney-p4vssin-venn_steep.tiff"),
        height = 900,
        width = 900,
        resolution = 300,
        col = "black",
        fill = c("dodgerblue", "goldenrod1"),
        alpha = 0.50,
        #cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
        #        1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
        #cex = c(rep(0.8, 4), rep(0.4,120)),
        #cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 0.8,
        cat.fontface = "bold",
        margin = 0.05
)
# Display the Venn
dev.off()
tiff_file <- paste0(WD,"/plots/20200505_rnaseq-kidney-p4vssin-venn_steep.tiff")
img <- readTIFF(tiff_file)
grid::grid.raster(img)

# Venn DIagram (now with CIRC Genes)
venn.plot <- venn.diagram(
        x = list(
                SIN = sinmod_df$ENSEMBL_RAT,
                POLY4 = p4mod_df$ENSEMBL_RAT,
                CIRC = circ_kid$ENSEMBL_RAT),
        filename = paste0(WD,"/plots/circ_steep.tiff"),
        height = 900,
        width = 900,
        resolution = 300,
        col = "black",
        fill = c("dodgerblue", "goldenrod1", "seagreen3"),
        alpha = 0.50,
        #cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
        #        1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
        #cex = c(rep(0.8, 4), rep(0.4,120)),
        #cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.cex = 0.8,
        cat.fontface = "bold",
        margin = 0.05
)


# Generate a scatterplot of the models' PVE's
shared_mod_df %>%
        ggplot(aes(x = PVE_P4, y= PVE_SIN)) +
        geom_point(alpha = 0.3) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("PVE: Polynomial (Degree 4)") +
        ylab("PVE: SIN")
# Supporting violin plot
shared_mod_df %>%
        select(ENSEMBL_RAT, PVE_SIN, PVE_P4) %>%
        gather(MODEL, PVE, -ENSEMBL_RAT) %>%
        ggplot(aes(MODEL,PVE)) +
        geom_violin() +
        geom_boxplot(alpha = 0.1) +
        scale_x_discrete(labels=c("PVE_SIN" = "SIN", 
                                  "PVE_P4" = "Polynomial (Degree 4)")) +
        scale_y_continuous(breaks=seq(0,1,0.2))

# Visualize the models on CIRC Genes
##############################
shared_mod_circ <- shared_mod_df %>%
        filter(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT) %>%
        mutate(SYMBOL_RAT = SYMBOL_RAT.x) %>%
        select(-SYMBOL_RAT.y, -SYMBOL_RAT.x) %>%
        arrange(desc(PVE_P4))
# Scatterplot
shared_mod_circ %>%
        ggplot(aes(x = PVE_P4, y= PVE_SIN)) +
        geom_point(alpha = 0.3) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("PVE: Polynomial (Degree 4)") +
        ylab("PVE: SIN")
# Supporting violin plot
shared_mod_circ %>%
        select(ENSEMBL_RAT, PVE_SIN, PVE_P4) %>%
        gather(MODEL, PVE, -ENSEMBL_RAT) %>%
        ggplot(aes(MODEL,PVE)) +
        geom_violin() +
        geom_boxplot(alpha = 0.1) +
        scale_x_discrete(labels=c("PVE_SIN" = "SIN", 
                                  "PVE_P4" = "Polynomial (Degree 4)")) +
        scale_y_continuous(breaks=seq(0,1,0.2))
dim(shared_mod_circ)
ensembl_genes <- shared_mod_circ$ENSEMBL_RAT %>% as.character()
ensembl <- ensembl_genes[4]
shared_mod_circ %>%
        select(SYMBOL_RAT,PVE_SIN,PVE_P4,PHASE_SIN)
ensembl <- "ENSRNOG00000047213"
for(ensembl in ensembl_genes){
        tod_counts_y <- tod_counts[ensembl,]
        range <- max(tod_counts[ensembl,])-min(tod_counts[ensembl,])
        med<- median(tod_counts_y)
        ordered_TOD<-TOD[order(TOD)]
        ordered_TPE<-TPE[order(TPE)]
        # Collect the phase
        phase <- shared_mod_circ %>%
                filter(ENSEMBL_RAT == ensembl) %>%
                select(PHASE_SIN) %>% unlist() %>% as.numeric()
        # Add the polynomial as line
        # Collect the betas
        B <- shared_mod_circ %>%
                filter(ENSEMBL_RAT == ensembl) %>%
                select(B1_P4,B2_P4,B3_P4,B4_P4,B5_P4) %>% unlist() %>% as.numeric()
        pred_y <- B[1] + (B[2]*(ordered_TPE)) + 
                (B[3]*(ordered_TPE)^2) +
                (B[4]*(ordered_TPE)^3) +
                (B[5]*(ordered_TPE)^4)
        SYMBOL <-  shared_mod_circ %>%
                filter(ENSEMBL_RAT == ensembl) %>%
                select(SYMBOL_RAT) %>% unlist() %>% as.character()
        # Plot the TOD counts
        plot(ordered_TOD,tod_counts_y, main = SYMBOL, ylab = "Expression", xlab = "Time of Death")
        # Add the sine function as line
        lines(ordered_TOD, med+0.4*range*sin(2*pi*(ordered_TOD + 3)/24),col=2)
        # Add the polynomial as line
        #lines(ordered_TPE, pred_y ,col=3)
        # Plot the TPE counts
        plot(ordered_TPE,tod_counts_y, main = SYMBOL, ylab = "Expression", xlab = "Hours Post Exercise")
        #Add the polynomial as line
        lines(ordered_TPE, pred_y ,col=3)
        # Add the sine function as line
        #lines(ordered_TOD, med+0.4*range*sin(2*pi*(ordered_TOD + 3)/24),col=2)
}
(shared_mod_circ$PVE_SIN < shared_mod_circ$PVE_P4) %>% table()






dim(circ_kid)


ggplot(df, aes(x = TOD, y= Exp)) +
        geom_point() +
        stat_smooth()
        







pval.cos.dlpfc2[pval.cos.dlpfc2 <= 0.05]



