#'---
#' title: "PASS1A Rat Tissue: -- Modeling of RNASeq Data Expression"
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
#'     * y = B_0 + B_1SIN(TOD) + B_2COS(TOD)
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
TISSUE <- "Kidney"
TIS <- "KID"
# Declare Outliers
if(TISSUE == "Kidney"){
        OUTLIERS <- c('90042016803_SF1')
}

#### Tissue:
print(TISSUE)
#### Outliers:
print(OUTLIERS)

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

#' #### Retrieve Circadian Genes Associated with Tissue
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script:

# Circadian Genes 
in_file <- paste0(WD,'/data/20200503_rnaseq-circadian-',TIS,'-mouse-rat-ortho_steep-yan.txt')
circ_df <- read.table(file=in_file, sep = '\t', header = TRUE)

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

# Filter Samples (meta)
tod_cols <- col_data %>%
        filter(Tissue == TISSUE) %>%
        filter(sample_key != OUTLIERS) %>%
        filter(!is.na(animal.registration.sex))
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

# Convert to seconds and take the square root (consider negative numbers)
tod_cols$specimen.collection.t_exercise_seconds <- 
        tod_cols$specimen.collection.t_exercise_hour * 60 * 60
tod_cols <- tod_cols %>%
        mutate(specimen.collection.t_exercise_hour_sqrt = ifelse(
                specimen.collection.t_exercise_seconds < 0, 
                (sqrt(abs(specimen.collection.t_exercise_seconds))/60/60)*(-1), 
                (sqrt(abs(specimen.collection.t_exercise_seconds))/60/60)))
tod_cols$specimen.collection.t_exercise_hour_sqrt_jit <- 
        jitter(tod_cols$specimen.collection.t_exercise_hour_sqrt, 
               factor = 0.1)
# Examine histograms
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=specimen.collection.t_exercise_hour)) +
        geom_histogram(bins = 68)
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=specimen.collection.t_exercise_seconds)) +
        geom_histogram(bins = 68)
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt)) +
        geom_histogram(bins = 68)
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt_jit)) +
        geom_histogram(bins = 68)

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

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Naive Model (~ 1)") +
        guides(color=guide_legend(title="animal.key.anirandgroup")) +
        scale_color_manual(values=ec_colors)

#' ### Adjust for Between Sex Variance

#+ Adjust for Between Sex Variance
################################################################################
########### Adjust for Between Sex Variance  ###################################
################################################################################

# "To adjust for batch effects, we median- centered the expression levels of each transcript within each batch and confirmed, using the correlation matrices, that the batch effects were removed after the adjustment." 
#~ Li, J. Z. et al. Circadian patterns of gene expression in the human brain and disruption in major depressive disorder. Proc. Natl. Acad. Sci. U. S. A. 110, 9950–9955 (2013).

# Here we have 2 Groups: Control - IPE and Control 7 hr; we'll median center these groups to combine the sexes.

M_samples <- col_data %>%
        filter(Tissue == TISSUE) %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Male') %>%
        filter(sample_key != OUTLIERS) %>%
        #filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
F_samples <- col_data %>%
        filter(Tissue == TISSUE) %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(sample_key != OUTLIERS) %>%
        filter(animal.registration.sex == 'Female') %>%
        #filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
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

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Median Centered by Sex") +
        guides(color=guide_legend(title="animal.key.anirandgroup")) +
        scale_color_manual(values=ec_colors)

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

# Order samples by exercise group (polynomial model)
group_order_poly <- tod_cols %>%
        filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
        arrange(factor(specimen.collection.t_exercise_hour), 
                desc(specimen.collection.t_exercise_hour), 
                desc(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Order of exercise groups
group_order_poly2 <- tod_cols %>%
        filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
        arrange(factor(specimen.collection.t_exercise_hour), 
                desc(specimen.collection.t_exercise_hour), 
                desc(animal.registration.sex)) %>%
        select(specimen.collection.t_exercise_hour) %>% unlist() %>% as.character()

# Arrange matrix in order of exercise and control groups
# All expression values per gene are normalized to generate a level playing field. Means are subtracted by values and the resulting difference is divided by the standard deviation.
som_xgroup <- assay(rld)[,group_order_poly]
som_xgroup <- som_xgroup %>% t() %>% data.frame()

# Test for multiv
#x <- as.numeric(MVN::mvn(som_xgroup, mvnTest = "hz")$univariateNormality$`p value`)
#(p.adjust(x, method = "BH") > 0.05) %>% table()
som_xgroup$specimen.collection.t_exercise_hour <- group_order_poly2

# MANOVA test
# Generate a formula
dependents <- colnames(som_xgroup)[colnames(som_xgroup) %!in% c("specimen.collection.t_exercise_hour")]
form <- as.formula(paste0("cbind(",paste(dependents, collapse=","),")", "~specimen.collection.t_exercise_hour"))
# Perform MANOVA
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

# Choose row index with significant p value
c_idx <- df_x %>%
        filter(pval <= 0.05) %>%
        select(idx) %>% unlist() %>% as.numeric()

# Select only genes that show significant variance across timepoints
som_xgroup <- som_xgroup[,c_idx]
dim(som_xgroup)

# Data should be normalized to have a mean zero and a variance of 1 (between -1 and 1)
preproc1 <- preProcess(som_xgroup, method=c("center", "scale"))
norm1 <- predict(preproc1, som_xgroup)
som_mat <- t(norm1) %>% as.matrix()

# Create a SOM (with som)
set.seed(100)
som_group <- som::som(som_mat,6,5)
plot(som_group, ylim=c(-2,2))

#som::som(som_xgroup,10,10) %>% plot()
table(som_group$visual[,1:2])

# Perform k-means clustering
k6<-kmeans(som_mat,6)$cluster

# Create an annotation for the heatmap
ann_df <- tod_cols %>%
        filter(sample_key %in% colnames(som_mat)) %>%
        select(specimen.collection.t_exercise_hour)
ann_df$specimen.collection.t_exercise_hour <- factor(ann_df$specimen.collection.t_exercise_hour)
ann_df <- ann_df %>% arrange(specimen.collection.t_exercise_hour)
row.names(ann_df) <- colnames(som_mat)
ann_colors = list(
        specimen.collection.t_exercise_hour = 
                c("-1" = "gold",
                  "0" = "darkgoldenrod1",
                  "0.5" = "orange",
                  "1" = "darkorange",
                  "4" = "darkorange2",
                  "7" = "darkorange3",
                  "24" = "steelblue1",
                  "48" = "steelblue4"))

# Create heatmaps from kmeans
heat_plots <- vector('list', 6)
heat_plots[[1]] <- pheatmap(som_mat[k6==1,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[2]] <- pheatmap(som_mat[k6==2,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[3]] <- pheatmap(som_mat[k6==3,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[4]] <- pheatmap(som_mat[k6==4,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[5]] <- pheatmap(som_mat[k6==5,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[6]] <- pheatmap(som_mat[k6==6,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)

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
        circ_cl[[kn]] <- circ_df %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(ENSEMBL_RAT) %>% unlist(use.names = FALSE)
        # collect additional information for circadian genes in cluster
        circ_cldf[[kn]] <- circ_df %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(NUM.TISSUE, ENSEMBL_RAT, SYMBOL_RAT)
        Z <- vector()
        for(i in 1:30){
                Z <- c(Z, (kgenes[kgenes %in% cluster[[i]]] %>% length()))
        }
        tile_df$Z <- Z
        # Generate the kmeans plot
        grad_plots[[kn]] <- local({
                kn <- kn
                kmsom <- ggplot(tile_df, aes(x = X, y = Y, fill =Z)) +
                        geom_raster(interpolate=TRUE) +
                        scale_fill_gradient2(low="navy", mid="white", high="red", 
                                             midpoint=30, limits=range(tile_df$Z),
                                             guide = FALSE) +
                        theme(axis.text        = element_blank(),
                              axis.ticks       = element_blank(),
                              axis.title       = element_blank(),
                              panel.background = element_blank())
                print(kmsom)
        })
}

# Create a dataframe for logistic regression
logreg_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
logreg_df <- mutate(logreg_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_df$ENSEMBL_RAT, 1, 0))
# Generate 6 columns for each kmeans cluster status
for(i in 1:6){
        new_col_name <- paste0('C',i)
        logreg_df <- logreg_df %>% 
                mutate(!!sym(new_col_name) := ifelse(ENSEMBL_RAT %in% kmeans_cluster[[i]], 1, 0))
}

# Remove gene ids and set them to rownames
row.names(logreg_df) <- logreg_df$ENSEMBL_RAT
logreg_df <-  logreg_df %>% select(-ENSEMBL_RAT)

# Adjust column objects
factor_cols <- colnames(logreg_df)
for(fc in factor_cols){
        logreg_df[[fc]] <- as.factor(logreg_df[[fc]])
}
str(logreg_df)

# Fit a logistic regression model
mod <- glm(formula = CIRC ~ C1 + C2 + C3 + C4 + C5 + C6, 
           family = binomial(link = logit), 
           data = logreg_df)
summary(mod)

# Perform a Chi-square analysis
#######################
# Create a proper data structure for chi squared
x2_list <- list()
for(i in 1:6){
        n_clust <- length(kmeans_cluster[[i]])
        n_circ <- length(circ_cl[[i]])
        n_no_circ <- n_clust - n_circ
        x2_list[[i]] <- c(n_circ, n_no_circ)
}

x2_df <- as.data.frame(x2_list)
row.names(x2_df) <- c('CIRC','NON-CIRC')
for(i in 1:6){
        new_col_name <- paste0('C',i)
        colnames(x2_df)[i] <- new_col_name
}

# Create the proper data structure for visualization of chi square with mosaic plot
mosaic_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
mosaic_df <- mutate(mosaic_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_df$ENSEMBL_RAT, 'CIRC', "NON-CIRC"))
# Generate a Cluster columcluster status
mosaic_df <- mosaic_df %>%
        mutate(CLUSTER = case_when(ENSEMBL_RAT %in% kmeans_cluster[[1]] ~ 1,
                                   ENSEMBL_RAT %in% kmeans_cluster[[2]] ~ 2,
                                   ENSEMBL_RAT %in% kmeans_cluster[[3]] ~ 3,
                                   ENSEMBL_RAT %in% kmeans_cluster[[4]] ~ 4,
                                   ENSEMBL_RAT %in% kmeans_cluster[[5]] ~ 5,
                                   ENSEMBL_RAT %in% kmeans_cluster[[6]] ~ 6))
# Generate a table of results
mosaic_tbl <- table(mosaic_df$CIRC, mosaic_df$CLUSTER, dnn=c("CIRC","CLUSTER"))

# Perform Chi-squared
#sink('')
(( c=chisq.test(x2_df, simulate.p.value = TRUE) ))
(( fisher.test(x2_df, simulate.p.value = TRUE) ))

#sink()
c$observed
round(c$expected)

# Visualize the chisquare analysis with a mosaic plot
mosaic(~ CIRC + CLUSTER, data = mosaic_tbl,
       shade = TRUE, legend = TRUE)

grad_plots[[4]]

# Clusters enriched WITH circadian genes
# 4
grad_plots[[4]]
heat_plots[[4]]
mypar()
plot(som_group, ylim=c(-2,2))
plot(som_group)

# Cluster enriched WITHOUT circadian genes
# 5
grad_plots[[5]]
heat_plots[[5]]
plot(som_group, ylim=c(-2,2))

# Plot the circadian genes in cluster 4
# Create a dataframe for the plot
df1 <- data.frame(t(som_mat))
df1$sample_key<- row.names(df1)
df2 <- tod_cols %>%
        select(specimen.collection.t_death, 
               specimen.collection.t_exercise_hour)
df2$sample_key <- row.names(df2)
# Adjust the column names to be symbols for ease of plot interpretation
df_plot <- left_join(df1, df2, by = "sample_key")
genes <- colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))]
symbols <- mapIds(org.Rn.eg.db, genes, "SYMBOL", "ENSEMBL")
colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))] <- symbols

# Melt the plot
melt_plot <- reshape2::melt(df_plot, id.vars = c("specimen.collection.t_death", 
                                                 "specimen.collection.t_exercise_hour"))
# Adjust columns and factor levels
melt_plot$value <- as.numeric(melt_plot$value)
melt_plot$specimen.collection.t_exercise_hour <- factor(melt_plot$specimen.collection.t_exercise_hour)

# To re-examine SOM
plot(som_group, ylim=c(-2,2))

# Plot the circadian genes for cluster of choice
# Show all plots they are enriched with and without in powerpoint
plot(grad_plots[[6]])
heat_plots[[6]]

melt_plot %>%
        filter(variable %in% circ_cldf[[3]]$SYMBOL_RAT) %>%
        ggplot(aes(
                x = as.integer(as.character(specimen.collection.t_exercise_hour)), 
                y = value, 
                color = variable),
               group = "1") +
        geom_point(alpha = 0.01) +
        #geom_smooth(alpha = 0.2, se = F, method = "loess") +
        #geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
        geom_line(stat = "smooth", alpha = 0.5, 
                  method = "lm", formula = y ~ ns(x, 4), se = T) +
        ylab("rlog Transformed Expression") +
        xlab("Hours Pre/Post Exercise") +
        theme(legend.position = "none")
        #geom_line(aes(group = "1")) +
        #scale_x_continuous(breaks = seq(-1,48,by=10))

#' #### Annotate Data for Modeling By Cluster

#+ Annotate Data for Modeling By Cluster
################################################################################
################ Annotate Data for Modeling By Cluster #########################
################################################################################

# Select the normailzed counts
tod_counts <- assay(rld) 
t_counts <- setNames(melt(tod_counts), 
                     c('ENSEMBL_RAT', 'sample_key', 'count'))

# Join the dataframes and nest
by_gene_df <- tod_cols %>%
        left_join(t_counts, by = "sample_key") %>%
        filter(animal.key.anirandgroup %!in% c('Control - 7 hr')) %>%
        group_by(ENSEMBL_RAT) %>%
        arrange(sample_key) %>%
        nest()
# Add Cluster and Circ Status
by_gene_df <- by_gene_df %>%
        left_join(mosaic_df, by = "ENSEMBL_RAT") %>%
        filter(!is.na(CIRC)) # This filter removes genes that did not pass variance filter (MANOVA)
# Add the gene symbol
by_gene_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")

# Join the dataframes and nest
by_gene_df7 <- tod_cols %>%
        left_join(t_counts, by = "sample_key") %>%
        filter(animal.key.anirandgroup %in% c('Control - 7 hr')) %>%
        group_by(ENSEMBL_RAT) %>%
        arrange(sample_key) %>%
        nest()
# Add Cluster and Circ Status
by_gene_df7 <- by_gene_df7 %>%
        left_join(mosaic_df, by = "ENSEMBL_RAT") %>%
        filter(!is.na(CIRC)) # This filter removes genes that did not pass variance filter (MANOVA)
# Add the gene symbol
by_gene_df7$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df7$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")

# Cleanup
#rm(count_data, col_data, counts_centered, dds ,dds1,dds2,F_centered,F_counts,
#   M_centered,M_counts,pcaData,t_counts)
# TODO: The entire dataframe cannot be unnested. It may be a memory issue or it may be a corrupted row -- come back and troubleshoot. For now, randomly select 2K rows to continue analysis.
#by_gene_df2 <- by_gene_df
#by_gene_df72 <- by_gene_df7
by_gene_df <- by_gene_df2
by_gene_df7 <- by_gene_df72

# Load in the DE genes between C0 and C7 for appropriate tissues
in_file <- paste0(WD, '/data/20200426_rnaseq-kidneyMF-C0C7DE-phaseMDD_steep.txt')
de_df <- read.table(in_file, sep = '\t', header = T)
# Take all the differentially expressed genes
de_all <- de_df %>%
        arrange(padj) %>%
        filter(padj <= 0.05) %>%
        select(ENSEMBL_RAT) %>%
        unlist() %>% as.character()
# Take the top 20 differentially expressed genes
de_top <- de_df %>%
        arrange(padj) %>%
        select(ENSEMBL_RAT) %>%
        dplyr::slice(1:50) %>%
        unlist() %>% as.character()
de_top20 <- de_df %>%
        arrange(padj) %>%
        select(ENSEMBL_RAT) %>%
        dplyr::slice(1:20) %>%
        unlist() %>% as.character()

#randomRows = sample(1:nrow(by_gene_df[,1]), 500, replace=F)
#by_gene_df <- by_gene_df[randomRows, ]
#by_gene_df7 <- by_gene_df7[randomRows, ]
by_gene_df <- by_gene_df %>% filter(ENSEMBL_RAT %in% de_all)
by_gene_df7 <- by_gene_df7 %>% filter(ENSEMBL_RAT %in% de_all)
# Must be true
all(by_gene_df$ENSEMBL_RAT == by_gene_df7$ENSEMBL_RAT)

# Generate model functions for the dataframes
gam_mod <- function(df) {
        lm(count ~ ns(specimen.collection.t_exercise_hour_sqrt_jit, df = 4), data = df)
}
# Generate a model function for the dataframes
sin_mod <- function(df) {
        lm(count ~ SIN(specimen.collection.t_death_hour) + 
                   COS(specimen.collection.t_death_hour),
           data = df)
}

# Generalized Additive Models
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
        mutate(gam_model = map(data, gam_mod))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
        mutate(gam_ANOVA = map(gam_model, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
        mutate(gam_resid = map2(data, gam_model, add_residuals))
# Examine the model metrics
gam_metrics <- by_gene_df %>%
        mutate(gam_metrics = map(gam_model, broom::glance)) %>%
        unnest(gam_metrics)
# Examine the model metrics
by_gene_df <- by_gene_df %>%
        mutate(gam_metrics = map(gam_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
        mutate(gam_summary = map(gam_model, summary))

# SIN/COS Model
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
        mutate(sin_model = map(data, sin_mod))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
        mutate(sin_ANOVA = map(sin_model, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
        mutate(sin_resid = map2(data, sin_model, add_residuals))
# Examine the model metrics
sin_metrics <- by_gene_df %>%
        mutate(sin_metrics = map(sin_model, broom::glance)) %>%
        unnest(sin_metrics)
# Examine the model metrics
by_gene_df <- by_gene_df %>%
        mutate(sin_metrics = map(sin_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
        mutate(sin_summary = map(sin_model, summary))

# Arrange the dataframe by model R2
row_order <- by_gene_df %>%
        unnest(gam_metrics) %>%
        arrange(desc(adj.r.squared)) %>%
        select(ENSEMBL_RAT) %>% 
        unlist() %>% as.character()
by_gene_df <- by_gene_df %>%
        ungroup() %>%
        arrange(factor(ENSEMBL_RAT, levels = row_order))

# Add the predictions (GAM model)
genes <- by_gene_df$ENSEMBL_RAT %>% as.character() %>% unique()
# Must be true
all(genes == by_gene_df$ENSEMBL_RAT)
gam_pred <- list()
i <- 1
for( g in genes){
        # Subset the dataframe by gene
        sub_data <- (by_gene_df %>% 
                             filter(ENSEMBL_RAT == g) %>%
                             ungroup() %>%
                             select(data))[[1]] %>% as.data.frame()
        # Generate a grid
        grid <- data.frame(specimen.collection.t_exercise_hour_sqrt_jit = 
                                   seq_range(
                                           sub_data$specimen.collection.t_exercise_hour_sqrt_jit, n = 68))
        #grid$ENSEMBL_RAT <- g
        mod <- (by_gene_df %>% 
                        filter(ENSEMBL_RAT == g) %>% ungroup() %>% 
                        select(gam_model))[[1]][[1]]
        summary(mod)
        grid <- add_predictions(grid, mod, "pred") %>% as_tibble()
        names(grid)[1] <- "grid_t_exercise_hour_sqrt_jit"
        grid$specimen.collection.t_exercise_hour_sqrt_jit <- 
                sub_data$specimen.collection.t_exercise_hour_sqrt_jit
        grid$count <- sub_data$count
        gam_pred[[i]] <- grid
        i <- i + 1
}
by_gene_df$gam_pred <- gam_pred

# Add the predictions (SIN Model)
genes <- by_gene_df$ENSEMBL_RAT %>% as.character() %>% unique()
# Must be true
all(genes == by_gene_df$ENSEMBL_RAT)
sin_pred <- list()
i <- 1
for( g in genes){
        # Subset the dataframe by gene
        sub_data <- (by_gene_df %>% 
                             filter(ENSEMBL_RAT == g) %>%
                             ungroup() %>%
                             select(data))[[1]] %>% as.data.frame()
        # Generate a grid
        grid <- data.frame(specimen.collection.t_death_hour = 
                                   seq_range(
                                           sub_data$specimen.collection.t_death_hour, n = 68))
        #grid$ENSEMBL_RAT <- g
        mod <- (by_gene_df %>% 
                        filter(ENSEMBL_RAT == g) %>% ungroup() %>% 
                        select(sin_model))[[1]][[1]]
        summary(mod)
        grid <- add_predictions(grid, mod, "pred") %>% as_tibble()
        names(grid)[1] <- "grid_t_death_hour"
        grid$grid_t_death_hour <- round(grid$grid_t_death_hour, digits = 1)
        grid$specimen.collection.t_death_hour <- 
                sub_data$specimen.collection.t_death_hour
        grid$count <- sub_data$count
        sin_pred[[i]] <- grid
        i <- i + 1
}
by_gene_df$sin_pred <- sin_pred

gene_n <- 5
# Visualize the GAM model by gene
gam_pred_df <- by_gene_df %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data, CIRC, CLUSTER, gam_model, gam_pred) %>%
        ungroup() %>%
        filter(row_number() == gene_n) %>%
        unnest(gam_pred)
gam_gene <- unique(gam_pred_df$ENSEMBL_RAT) %>% as.character()
# Visualize the SIN model by gene
sin_pred_df <- by_gene_df %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data, CIRC, CLUSTER, sin_model, sin_pred) %>%
        ungroup() %>%
        filter(row_number() == gene_n) %>%
        unnest(sin_pred)
sin_gene <- unique(sin_pred_df$ENSEMBL_RAT) %>% as.character()
# Must be true
gam_gene == sin_gene

# Collect the Control 7 hr data points
gam_hr7_df <- by_gene_df7 %>%
        filter(ENSEMBL_RAT == gam_gene) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data) %>%
        ungroup() %>%
        unnest(data) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, specimen.collection.t_exercise_hour_sqrt_jit, count)

# Collect the hours of death that occur for each exercise group
sac_hrs <- sort(unique(round(tod_cols$specimen.collection.t_death_hour, digits = 1)))
#names(tod_cols)
#e_df <- tod_cols %>%
#        select(specimen.collection.t_death_hour, specimen.collection.t_exercise_hour_sqrt) %>%
#        arrange(specimen.collection.t_death_hour)
#e_df$specimen.collection.t_death_hour <- round(e_df$specimen.collection.t_death_hour, digits = 1)

# Perform a second prediction (for hour post exercise -- predictions from SIN Model)
sin_pred_hour <- sin_pred_df %>%
        filter(round(grid_t_death_hour, digits = 1) %in% sac_hrs) %>%
        mutate(grid_t_exercise_hour = case_when(grid_t_death_hour == 9.9 ~ -0.01666667,
                                                grid_t_death_hour == 10.0 ~ -0.01666667,
                                                grid_t_death_hour == 10.2 ~ -0.01666667,
                                                grid_t_death_hour == 10.3 ~ -0.01666667,
                                                grid_t_death_hour == 10.6 ~ 0.08164966,
                                                grid_t_death_hour == 10.9 ~ 0.08164966,
                                                grid_t_death_hour == 11.2 ~ 0.11547005,
                                                grid_t_death_hour == 11.6 ~ 0.11547005,
                                                grid_t_death_hour == 11.9 ~ 0.01178511,
                                                grid_t_death_hour == 12.2 ~ 0.01178511,
                                                grid_t_death_hour == 12.3 ~ 0.01178511,
                                                grid_t_death_hour == 13.2 ~ 0.00000000,
                                                grid_t_death_hour == 13.3 ~ 0.00000000,
                                                grid_t_death_hour == 13.6 ~ 0.00000000,
                                                grid_t_death_hour == 13.9 ~ 0.01666667,
                                                grid_t_death_hour == 14.2 ~ 0.01666667,
                                                grid_t_death_hour == 14.3 ~ 0.01666667,
                                                grid_t_death_hour == 14.6 ~ 0.03333333,
                                                grid_t_death_hour == 14.9 ~ 0.03333333,
                                                grid_t_death_hour == 16.9 ~ 0.04409586,
                                                grid_t_death_hour == 17.2 ~ 0.04409586,
                                                grid_t_death_hour == 17.3 ~ 0.04409586,
                                                grid_t_death_hour == 17.6 ~ 0.04409586,
                                                grid_t_death_hour == 17.9 ~ 0.04409586))

# Visualize the raw counts, model predictions, and control 7 counts (x=Hours Post Exercise)
gam_pred_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt_jit, count), 
                   color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        geom_point(data = gam_hr7_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt_jit, count),
                   color = "red") +
        geom_point(data = sin_pred_hour, 
                  mapping = aes(x = grid_t_exercise_hour, y = pred), 
                  size = 3, alpha = 1, color = "orange") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(gam_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Hours Post Acute Exercise (Transformed)") +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

# Visualize the raw counts, model predictions, and control 7 counts (x=Hours Post Exercise)
# Collect the Control 7 hr data points
sin_hr7_df <- by_gene_df7 %>%
        filter(ENSEMBL_RAT == sin_gene) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data) %>%
        ungroup() %>%
        unnest(data) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, specimen.collection.t_death_hour, count)
# Visualize the raw counts, model predictions, and control 7 counts (x=Hours Post Exercise)
sin_pred_df %>%
        ggplot(aes(specimen.collection.t_death_hour, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_point(data = sin_hr7_df, 
                   mapping = aes(specimen.collection.t_death_hour, count),
                   color = "red") +
        #geom_point(data = sin_pred_hour, 
        #           mapping = aes(x = grid_t_death_hour, y = pred), 
        #           size = 1, alpha = 1, color = "orange") +
        geom_line(data = sin_pred_df, 
                  aes(grid_t_death_hour, pred), 
                  size = 1, alpha = 0.8, color = "orange") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(sin_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Time of Death")

# Visualize the gam model metrics
(by_gene_df %>%
        dplyr::slice(gene_n) %>%
        select(gam_model) %>%
        ungroup %>% 
        select(gam_model))[[1]][[1]] %>%
        anova()
by_gene_df[gene_n,"gam_metrics"][[1]][[1]]
by_gene_df[gene_n,"sin_metrics"][[1]][[1]]
# Visualize the model ANOVA summary
#by_gene_df$gam_ANOVA[[1]]
#by_gene_df$sin_ANOVA[[1]]
################################################################################

# Compare the R2 values and p values from models
sin_metrics$r.squared.sin <- sin_metrics$r.squared
sin_metrics$p.value.sin <- sin_metrics$p.value  
gam_metrics$r.squared.gam <- gam_metrics$r.squared
gam_metrics$p.value.gam <- gam_metrics$p.value 
model_metrics <- sin_metrics %>%
        select(ENSEMBL_RAT, r.squared.sin, p.value.sin) %>%
        left_join(gam_metrics, by = "ENSEMBL_RAT") %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, CIRC, CLUSTER, r.squared.sin, p.value.sin,
               r.squared.gam, p.value.gam)

# Join the log fold change for sake of graph annotation
de_join <- de_df %>% 
        select(ENSEMBL_RAT, log2FoldChange)
# Annotate the top 20 de genes
model_metrics <- model_metrics %>%
        ungroup() %>% 
        left_join(de_join, by = "ENSEMBL_RAT") %>%
        mutate(Expression = factor(ifelse(log2FoldChange > 0, "UP", "DOWN"),
                                   levels = c("UP", "DOWN"))) %>%
        mutate(top_de = ifelse(ENSEMBL_RAT %in% de_top, "TOP GENE", "NOT TOP GENE"))
top_metrics <- model_metrics %>%
        filter(ENSEMBL_RAT %in% de_top)
top20_metrics <- model_metrics %>%
        filter(ENSEMBL_RAT %in% de_top20)
circ_metrics <- model_metrics %>%
        filter(CIRC == "CIRC")

# Compare the R2 between plots
# First show circadian genes
ggplot(model_metrics, aes(r.squared.gam, r.squared.sin)) +
        geom_point(alpha = 0.1) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Natural Spline Model (Exercise)") +
        ylab("R^2 SIN/COS Model (Circadian)") +
        ggtitle("R2 Comparisons Between Models:\nDifferentially Expressed Genes (C0 -> C7)")
ggplot(model_metrics, aes(r.squared.gam, r.squared.sin)) +
        geom_point(alpha = 0.1) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Natural Spline Model (Exercise)") +
        ylab("R^2 SIN/COS Model (Circadian)") +
        ggtitle("R2 Comparisons Between Models:\nDifferentially Expressed Genes (C0 -> C7)") +
        geom_point(data = circ_metrics,
                   mapping = aes(r.squared.gam, r.squared.sin), 
                   alpha = 1)

# Map the top DE genes
ggplot(model_metrics, aes(r.squared.gam, r.squared.sin)) +
        geom_point(alpha = 0.1) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Natural Spline Model (Exercise)") +
        ylab("R^2 SIN/COS Model (Circadian)") +
        geom_point(data = top_metrics,
                   mapping = aes(r.squared.gam, r.squared.sin, color = Expression), 
                   alpha = 1)
# Label the top DE genes
ggplot(model_metrics, aes(r.squared.gam, r.squared.sin)) +
        geom_point(alpha = 0.1) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Natural Spline Model (Exercise)") +
        ylab("R^2 SIN/COS Model (Circadian)") +
        geom_point(data = top_metrics,
                   mapping = aes(r.squared.gam, r.squared.sin, color = Expression), 
                   alpha = 1)+
        geom_label_repel(data = top20_metrics,
                         mapping = aes(label=SYMBOL_RAT), alpha = 0.8, 
                         hjust=0, vjust=0)























# COmpare the p value between models
ggplot(model_metrics, aes(-log(p.value.poly), -log(p.value.circ))) +
        geom_point(alpha = 0.05) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("-Log(p-value) Polynomial Spline Model") +
        ylab("-Log(p-value) Circadian Model")
ggplot(model_metrics, 
       aes(-log(p.value.poly), -log(p.value.circ), color = CIRC)) +
        geom_point(alpha = 0.3) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("-Log(p-value) Polynomial Spline Model") +
        ylab("-Log(p-value) Circadian Model") +
        scale_color_manual(breaks = c("CIRC", "NON-CIRC"),
                           values=c("red", "white"))

names(poly_resid_df)[grepl("mod",names(poly_resid_df))]

# Visualize the model
poly_resid_df %>%
        data_grid(count, 
                  specimen.collection.t_exercise_hour = 
                          seq_range(specimen.collection.t_exercise_hour, n = 8)) %>%
        add_predictions(poly_model)
seq_range(poly_resid_df$specimen.collection.t_exercise_hour, n=100)

# Now we can plot the residuals
# As a grand model
poly_resid_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour, resid)) +
        geom_line(aes(group = ENSEMBL_RAT), alpha = 1/3) +
        geom_smooth(se = F)
# Facetted by k means cluster
poly_resid_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour, resid, 
                   group = ENSEMBL_RAT)) +
        geom_line(alpha = 1/3) +
        facet_wrap(~CLUSTER)
# Demonstrate the residual for a particular gene of interest
genes <- poly_resid_df %>%
        select(ENSEMBL_RAT) %>%
        unlist() %>% as.character()
poly_resid_df %>%
        filter(ENSEMBL_RAT == genes[1]) %>%
        ggplot(aes(specimen.collection.t_exercise_hour, resid)) +
        geom_line(aes(group = ENSEMBL_RAT), alpha = 1/3) +
        geom_smooth(se = F)

grid <- poly_pred_df[1:8,] %>%
        data_grid(
                specimen.collection.t_exercise_hour = 
                        seq_range(specimen.collection.t_exercise_hour, 8))
grid <- add_predictions(grid, poly_pred_df[1:8,]$poly_model[[1]], "count")

poly_pred_df[1:8,] %>%
        ggplot(aes(specimen.collection.t_exercise_hour, count)) +
        geom_point() +
        geom_line(data = grid, color = "red", size = 1)



        
by_gene_df$data[[1]]

############
# Investigation of poly model


# Examine the model metrics
circ_metrics <- by_gene_df %>%
        mutate(circ_metrics = map(circ_model, broom::glance)) %>%
        #select(-data,-circ_model,-poly_ns4_model,
        #-circ_resid,-poly_ns4_resid) %>%
        unnest(circ_metrics)

dual_metrics <- by_gene_df %>%
        mutate(dual_metrics = map(dual_model, broom::glance)) %>%
        #select(-data,-circ_model,-poly_ns4_model,
        #-circ_resid,-poly_ns4_resid) %>%
        unnest(dual_metrics)

# Examine the R_2 for CIRC and NON-CIRC genes
# CIRC Model
circ_metrics %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_jitter(width = 0.5, alpha = 0.2)
# For a random set of noncirc genes equal in n to circ gene set
circ_set <- circ_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character()
noncirc_set <- circ_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'NON-CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character() %>%
        sample(length(circ_set))
circ_metrics_sub <- circ_metrics %>%
        filter(ENSEMBL_RAT %in% c(circ_set, noncirc_set))
# Examine the R_2 for CIRC and NON-CIRC genes (SUBSET of NON-CIRC)
circ_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_violin() +
        geom_boxplot(alpha = 0.4)

# Examine the R_2 for CIRC and NON-CIRC genes
# Polynomial model
poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_jitter(width = 0.5, alpha = 0.2)
# For a random set of noncirc genes equal in n to circ gene set
circ_set <- poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character()
noncirc_set <- poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'NON-CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character() %>%
        sample(length(circ_set))
poly_ns4_metrics_sub <- poly_ns4_metrics %>%
        filter(ENSEMBL_RAT %in% c(circ_set, noncirc_set))
# Examine the R_2 for CIRC and NON-CIRC genes (SUBSET of NON-CIRC)
poly_ns4_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_violin() +
        geom_boxplot(alpha = 0.4)

# Examine the R_2 for CIRC and NON-CIRC genes
# Dual model
dual_metrics %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_jitter(width = 0.5, alpha = 0.2)
# For a random set of noncirc genes equal in n to circ gene set
circ_set <- dual_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character()
noncirc_set <- dual_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'NON-CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character() %>%
        sample(length(circ_set))
dual_metrics_sub <- dual_metrics %>%
        filter(ENSEMBL_RAT %in% c(circ_set, noncirc_set))
# Examine the R_2 for CIRC and NON-CIRC genes (SUBSET of NON-CIRC)
dual_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_violin() +
        geom_boxplot(alpha = 0.4)

# Examine the R_2 for CIRC and NON-CIRC genes
# Dual
poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_jitter(width = 0.5, alpha = 0.2)
# For a random set of noncirc genes equal in n to circ gene set
circ_set <- poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character()
noncirc_set <- poly_ns4_metrics %>%
        filter(p.value <= 0.017) %>%
        filter(CIRC == 'NON-CIRC') %>%
        select(ENSEMBL_RAT) %>%
        unique() %>% unlist() %>% as.character() %>%
        sample(length(circ_set))
poly_ns4_metrics_sub <- poly_ns4_metrics %>%
        filter(ENSEMBL_RAT %in% c(circ_set, noncirc_set))
# Examine the R_2 for CIRC and NON-CIRC genes (SUBSET of NON-CIRC)
poly_ns4_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        as.data.frame() %>% as_tibble() %>%
        ggplot(aes(CIRC, r.squared)) +
        geom_violin() +
        geom_boxplot(alpha = 0.4)

# Unnest the dataframe to collect the residuals
circ_resid_df <- circ_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        unnest(circ_resid)
# Plot the residuals and facet
circ_resid_df %>%
        ggplot(aes(specimen.collection.t_death_hour, resid)) +
        geom_line(aes(group = ENSEMBL_RAT), alpha = 0.1) +
        geom_smooth(se = FALSE) +
        facet_wrap(~CIRC)

# Unnest the dataframe to collect the residuals
poly_ns4_resid_df <- poly_ns4_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        unnest(poly_ns4_resid)
# Plot the residuals and facet
poly_ns4_resid_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour, resid)) +
        geom_line(aes(group = ENSEMBL_RAT), alpha = 0.1) +
        geom_smooth(se = FALSE) +
        facet_wrap(~CIRC)

# Unnest the dataframe to collect the residuals
dual_resid_df <- dual_metrics_sub %>%
        filter(p.value <= 0.017) %>%
        unnest(dual_resid)
# Plot the residuals and facet
dual_resid_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour, resid)) +
        geom_line(aes(group = ENSEMBL_RAT), alpha = 0.1) +
        geom_smooth(se = FALSE) +
        facet_wrap(~CIRC)

# Compare the R2 values and p values from models
poly_ns4_metrics$r.squared.poly <- poly_ns4_metrics$r.squared
poly_ns4_metrics$p.value.poly <- poly_ns4_metrics$p.value  
circ_metrics$r.squared.circ <- circ_metrics$r.squared
circ_metrics$p.value.circ <- circ_metrics$p.value 
model_metrics <- poly_ns4_metrics %>%
        select(ENSEMBL_RAT, r.squared.poly, p.value.poly) %>%
        left_join(circ_metrics, by = "ENSEMBL_RAT") %>%
        select(ENSEMBL_RAT, CIRC, r.squared.poly, p.value.poly,
               r.squared.circ, p.value.circ)
# Compare the R2 between plots
ggplot(model_metrics, aes(r.squared.poly, r.squared.circ)) +
        geom_point(alpha = 0.05) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Polynomial Spline Model") +
        ylab("R^2 Circadian Model")
ggplot(model_metrics, aes(r.squared.poly, r.squared.circ, color = CIRC)) +
        geom_point(alpha = 0.25) +
        xlim(0,1) + ylim(0,1) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("R^2 Polynomial Spline Model") +
        ylab("R^2 Circadian Model") +
        scale_color_manual(breaks = c("CIRC", "NON-CIRC"),
                           values=c("black", "grey"))

# COmpare the p value between models
ggplot(model_metrics, aes(-log(p.value.poly), -log(p.value.circ))) +
        geom_point(alpha = 0.05) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("-Log(p-value) Polynomial Spline Model") +
        ylab("-Log(p-value) Circadian Model")
ggplot(model_metrics, 
       aes(-log(p.value.poly), -log(p.value.circ), color = CIRC)) +
        geom_point(alpha = 0.3) +
        geom_abline(intercept = 0, slope = 1) +
        xlab("-Log(p-value) Polynomial Spline Model") +
        ylab("-Log(p-value) Circadian Model") +
        scale_color_manual(breaks = c("CIRC", "NON-CIRC"),
                           values=c("red", "white"))

################################

qqnorm( -log(model_metrics$p.value.poly), main='Splines Model')
qqline( -log(model_metrics$p.value.poly) )

qqnorm( -log(model_metrics$p.value.circ), main='Circadian Model')
qqline( -log(model_metrics$p.value.circ) )






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

# Save the phase values
out_file <- paste0(WD,'/data/20200505_rnaseq-bothsexes-liver-phases_steep.txt')
#write.table(sinmod_df, file = out_file, quote = FALSE, 
#            row.names = FALSE, col.names = TRUE, sep ='\t')

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
                animal.key.anirandgroup == 'Exercise - 48 hr' ~ -1)) %>%
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
        filename = paste0(WD,"/plots/20200505_rnaseq-liver-p4vssin-venn_steep.tiff"),
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
tiff_file <- paste0(WD,"/plots/20200505_rnaseq-liver-p4vssin-venn_steep.tiff")
img <- readTIFF(tiff_file)
grid::grid.raster(img)

# Venn DIagram (now with CIRC Genes)
venn.plot <- venn.diagram(
        x = list(
                SIN = sinmod_df$ENSEMBL_RAT,
                POLY4 = p4mod_df$ENSEMBL_RAT,
                CIRC = circ_df$ENSEMBL_RAT),
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
        filter(ENSEMBL_RAT %in% circ_df$ENSEMBL_RAT) %>%
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






dim(circ_df)


ggplot(df, aes(x = TOD, y= Exp)) +
        geom_point() +
        stat_smooth()








pval.cos.dlpfc2[pval.cos.dlpfc2 <= 0.05]



