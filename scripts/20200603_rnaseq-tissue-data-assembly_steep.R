'---
#' title: "PASS1A Rat Tissue: -- Data Assembly"
#' author: "Alec Steep & Jiayu Zhang" 
#' date: "20200604"
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
#' * Assembly Data
#' * Clean Outliers
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
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("EnhancedVolcano")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","tiff","splines","gam","EnhancedVolcano")
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

# TISSUE: Hypothalamus, Liver, Kidney, Aorta, Adrenal, Brown Adipose, Cortex, Gastrocnemius, Heart, Hippocampus,Lung,Ovaries,PaxGene,Spleen,Testes, White Adipose

# SCN: Suprachiasmatic nucleus -- Hypothalamus
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
TISSUE <- "Adrenal"
TIS <- "ADG"
OUTLIERS <- c('None')
MIS_ID <- c('None')
MIS_ID_Fw2Mc <- c('None')
# Declare Outliers
if(TISSUE == "Kidney"){
        #OUTLIERS <- c('90042016803_SF1','90109015902_SN1')
        #OUTLIERS <- c('90109015902_SN1')
        OUTLIERS <- c('None')
}else if(TISSUE == "Liver"){
        OUTLIERS <- c('90042016803_SF1')
}else if(TISSUE == "Hypothalamus"){
        OUTLIERS <- c('None')
        MIS_ID <- c('90010015402_SF2','90011015402_SF2')
        PRI_VAR <- c("animal.registration.sex")
}else if(TISSUE == "Aorta"){
        OUTLIERS <- c('None')
}else if(TISSUE == "Gastrocnemius"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Heart"){
        OUTLIERS <- c('90052015802_SN1')
}else if(TISSUE == "Adrenal"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Brown Adipose"){
        OUTLIERS <- c('None')
}else if(TISSUE == "White Adipose"){
        OUTLIERS <- c('90127017003_SF1','90033017003_SF1')
}else if(TISSUE == "Cortex"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Hippocampus"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Lung"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Ovaries"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Testes"){
        OUTLIERS <- c('90127016302_SN2','90015016302_SN2')
}else if(TISSUE == "Spleen"){
        OUTLIERS <- c('ABC')
}

#### Tissue:
print(TISSUE)
#### Outliers:
print(OUTLIERS)
print(MIS_ID)
print(PRI_VAR)

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
circ_df <- read.table(file=in_file, sep = '\t', header = TRUE) %>% as_tibble()

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
# Make sure to Subtract 1 hour (3600s) from "Control - IPE" groups to account for exercise effect
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
################################################################################
# To determine object size
sl <- object.size(count_data)
print(sl, units = "auto") 
# Meta Data
out_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
#saveRDS(col_data, file = out_file)
# Count Data
out_file <- paste0(WD, '/data/20200603_rnaseq-meta-pass1a-stanford-sinai-processed_steep.rds')
#saveRDS(count_data, file = out_file)

# Restore the object
#readRDS(file = "my_data.rds")

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
ens_gtf <- paste0(WD,'/data/Rattus_norvegicus.Rnor_6.0.96.gtf')
Rn_TxDb <- makeTxDbFromGFF(ens_gtf,
                           format=c("gtf"),
                           dataSource="Ensembl_Rattus6_gtf",
                           organism="Rattus norvegicus",
                           taxonomyId=NA,
                           circ_seqs=DEFAULT_CIRC_SEQS,
                           chrominfo=NULL,
                           miRBaseBuild=NA,
                           metadata=NULL)

# Define Female specific sex genes (X chromosome)
# To examine chromosome names
seqlevels(Rn_TxDb)[1:23]
# Extract genes as GRanges object, then names
X_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("X")))
# Collect ensembl gene ids for female specific genes
X_ens_id <- names(X_genes_gr)
# Examine the gene symbols
X_sym <- mapIds(org.Rn.eg.db, names(X_genes_gr), "SYMBOL", "ENSEMBL")
# Extract genes as GRanges object, then names
Y_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("Y")))
# Collect ensembl gene ids for female specific genes
Y_ens_id <- names(Y_genes_gr)
sex_ens_id <- c(X_ens_id,Y_ens_id)
# Examine the gene symbols
Y_sym <- mapIds(org.Rn.eg.db, names(Y_genes_gr), "SYMBOL", "ENSEMBL")

#' ## Collect Samples of Interest and Normalize

#+ Collect Samples of Interest and Normalize
################################################################################
#####     Collect Samples of Interest and Normalize      #######################
################################################################################

# Filter Samples (meta)
tod_cols <- col_data %>%
        filter(Tissue == TISSUE) %>%
        #filter(sample_key != OUTLIERS) %>%
        filter(!is.na(animal.registration.sex))
rownames(tod_cols) <- tod_cols$sample_key

# Collect samples without NA values in TOD
nona_sams <- tod_cols %>%
        filter(!is.na(specimen.collection.t_death_hour)) %>%
        #filter(sample_key != OUTLIERS) %>%
        filter(!is.na(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Collect tissue specific counts
tod_counts <- count_data[,nona_sams]

#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)

design = ~ animal.registration.sex + animal.key.anirandgroup # Primary variable needs to be last.
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

#' ## Qualitative Assessment of Variance in Tissue

#+ Qualitative Assessment of Variance in Tissue
################################################################################
#####     Qualitative Assessment of Variance in Tissue      ####################
################################################################################

# Examine the Annotations most strongly associated with PCs
################################################################################
# Function to grab topp 500 variance genes and find correlated variables
###########
cor_PC_1_4 <- function(rld = rld, ntop = 500, intgroups){
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  out_df <- data.frame(matrix(ncol = 3, nrow = 0))
  names(out_df) <- c("PC","Adjusted_R_Sq","Condition")
  for(intgroup in intgroups) {
    intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
    group <- colData(rld)[[intgroup]]
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],
                    group = group, 
                    intgroup.df, name = colnames(rld))
    # Perform linear model on PC vs
    if(length(unique(group[!is.na(group)])) >=2){
      r2_pc1 <- (lm(PC1 ~ group, data = d) %>% summary())$adj.r.squared
      r2_pc2 <- (lm(PC2 ~ group, data = d) %>% summary())$adj.r.squared
      r2_pc3 <- (lm(PC3 ~ group, data = d) %>% summary())$adj.r.squared
      r2_pc4 <- (lm(PC4 ~ group, data = d) %>% summary())$adj.r.squared
      int_df <- data.frame(PC = 1:4, 
                           Adjusted_R_Sq = c(r2_pc1, r2_pc2,r2_pc3, r2_pc4),
                           Condition = intgroup) %>% as_tibble()
      out_df <- rbind(out_df, int_df)
    }else{
      NA
    }
  }
  out_df %>%
    arrange(desc(Adjusted_R_Sq))
  
}
#########
# Find variables associated with PCs
pc_cor_df <- cor_PC_1_4(rld = rld, ntop = 500, intgroups = names(tod_cols)) %>%
  filter(Adjusted_R_Sq > 0.3) %>%
  arrange(PC)
# Examine Variables of interest for Pcs 1 and 2
pc_cor_df %>%
  filter(PC %in% c(1,2))

# Visualize variables in PC's 1 and 2
voi <- pc_cor_df %>%
  filter(PC %in% c(1,2)) %>%
  select(Condition) %>% unlist() %>% as.character
# Plot the top variables associated with PCs 1 or 2
for( v in voi){
  # Quick investigation of variables
  plot(DESeq2::plotPCA(rld, intgroup =v) +
    guides(color=guide_legend(title=v)))
}

#' #### Do Tissue Samples cluster by sex.
#' Grey samples represent "reference samples".
DESeq2::plotPCA(rld, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="Sex"))

################################################################################
###### Manual Annotation Opportunity: Fill in PRI_VAR variable #################
################################################################################

#for(pri_var in PRI_VAR){
  pri_var <- "animal.registration.sex"
#}

# Examine a labled PCA Plot
################################################################################
#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c(pri_var, "sample_key","animal.key.anirandgroup"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf(paste0(WD,"/plots/20200505_rnaseq-",TIS,"-PCA-naive-modeling_steep.pdf"),
#    width = 6, height = 4)
pcaData %>%
  ggplot(aes_string("PC1", "PC2", color=pri_var)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\nNaive Model (~ 1)")) +
        guides(color=guide_legend(title=pri_var)) +
        theme(legend.title=element_blank())
#dev.off()

# Plot Sex-Specific Plots
################################################################################
if(pri_var == "animal.registration.sex"){
  
}
# Variables of interest
male_tis <- (tod_cols %>% 
                     filter(animal.registration.sex == 'Male'))$sample_key %>% 
        as.character()
female_tis <- (tod_cols %>% 
                       filter(animal.registration.sex == 'Female'))$sample_key %>% 
        as.character()
Y_genes <- Y_ens_id[Y_ens_id %in% row.names(assay(rld))]
X_genes <- X_ens_id[X_ens_id %in% row.names(assay(rld))]
sex <- tod_cols$animal.registration.sex
group <- tod_cols$animal.registration.sex

#' #### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.
chryexp <- colMeans(assay(rld)[Y_genes,])
chryexp_df <- data.frame("counts" = chryexp)
chryexp_df$sample <- row.names(chryexp_df) %>% as.factor()
chryexp_df <- chryexp_df %>%
        mutate(sex = ifelse(sample %in% male_tis, 'Male', 'Female'))

# Visualize genes on the Y Chromosome
##############################################
#' ##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution.
mypar()
chryexp_df %>%
  mutate(plot_labs = ifelse(sample %in% MIS_ID, as.character(sample), '')) %>%
  ggplot(aes(counts, colour = sex)) +
  geom_freqpoly(alpha = 0.6) +
  xlab("Median counts on Y genes (normalized)") +
  ylab("Frequency") +
  ggtitle("Frequency Polyplot & Sample Dotplot: \nSamples plotted by Y gene coverage (median)") +
  geom_dotplot(aes(fill = sex), alpha = 0.5) +
  geom_text_repel(aes(label = plot_labs), box.padding = unit(2, 'lines'), y = 0)

################################################################################
###### Manual Annotation Opportunity: Fill in MIS_ID variable ##################
################################################################################

#' ## Sample Classification: High Dimension Outlier Prediction + LDA

#+ Sample Classification: High Dimension Outlier Prediction + LDA 
################################################################################
#####     Sample Classification: High Dimension Outlier Prediction + LDA #######
################################################################################

# (iii) standardization of the data so that the expression measures for each array have mean 0 and variance 1 across genes.
edata <- scale(assay(rld))

# check that we get mean of 0 and sd of 1
apply(edata , 2, mean)
apply(edata , 2, sd)

#' ## Reassign Annotation Groups to MisIdentified Samples: Outlier Detection (Per Gene)

#+ Reassign Annotation Groups to MisIdentified Samples
################################################################################
#####     Reassign Annotation Groups to MisIdentified Samples      #############
################################################################################
# TODO: Add a conditional statement if samples we misidentieid or not

# identify multivariate outliers
x.out <- pcout(edata, makeplot = T)

# Create a dataframe with 3 columns: gene, weight_combined, gene_index
df_out <- data.frame(names(x.out$wfinal), x.out$wfinal) %>% as_tibble()
names(df_out) <- c("gene", 'weight_combined')
df_out$index <- row.names(df_out)

# Arrange the dataframe by weight
df_out <- df_out %>%
  arrange(weight_combined)

# Plot the Outliers by combined weight
df_out %>%
  sample_n(300) %>%
  ggplot(aes(index, weight_combined)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0.1)

# Number of strong outlier genes
genes_out <- df_out %>%
  filter(weight_combined < 0.10) %>%
  select(gene) %>% unlist() %>% as.character()
  
# LDA Analysis, Predict Sex with outlier gene set
set.seed(100)
rld_out <- rld[names(rld) %in% genes_out,]
# Remove the Sample in question
rld_mis <- rld_out[,rld_out$sample_key %in% MIS_ID]
rld_out <- rld_out[,rld_out$sample_key %!in% MIS_ID]

# Create a train and test set
index = createDataPartition(y=rld_out$animal.registration.sex, p=0.7, list=FALSE)
sample_train <- as.character(rld$sample_key[index])
sample_test <- as.character(rld$sample_key[-index])
rld_train = rld_out[,rld_out$sample_key %in% sample_train]
rld_test = rld_out[,rld_out$sample_key %in% sample_test]
# Reorganize the train and test DFs
join_df <- tod_cols %>%
  select(sample_key, animal.registration.sex)
# Train
train_df <- data.frame(t(assay(rld_train)))
train_df$sample_key <- row.names(train_df)
train_df <- train_df %>%
  left_join(join_df, by = c("sample_key")) %>%
  as_tibble() %>%
  select(-sample_key)
# Test
test_df <- data.frame(t(assay(rld_test)))
test_df$sample_key <- row.names(test_df)
test_df <- test_df %>%
  left_join(join_df, by = c("sample_key")) %>%
  as_tibble() %>%
  select(-sample_key)
# Misidentified
mis_df <- data.frame(t(assay(rld_mis)))
mis_df$sample_key <- row.names(mis_df)
mis_df <- mis_df %>%
  left_join(join_df, by = c("sample_key")) %>%
  as_tibble() %>%
  select(-sample_key)

# genes to remove from model
x = nearZeroVar(train_df, saveMetrics = TRUE)
str(x, vec.len=2)
# The Zero variance predictors
x[x[,"zeroVar"] > 0, ]
# The near Zero variance predictors
x[x[,"zeroVar"] + x[,"nzv"] > 0, ]

# prepare training scheme
control <- trainControl(method="repeatedcv", number = 10, repeats = 10)
#control <- trainControl(method="LOOCV")
# Fit the LDA
lda.fit = train(animal.registration.sex ~ ., data=train_df, method="lda",
                trControl = control)

# Random Forest
######################
# Set up the random forest tuning grid
grid_rf <- expand.grid(mtry = c(1000,2000,5000,10000))
set.seed(100)
rf.fit <- train(animal.registration.sex ~ ., data=train_df,
                method = "rf",
                ntree = 1000,
                trControl = trainControl(method="repeatedcv", number = 10, repeats = 10),
                tuneGrid = grid_rf)

# Examine model
rf.fit

# Visualize feature importance
rfImp <- varImp(rf.fit, scale = FALSE)
plot(rfImp, top = 25, scales = list(y = list(cex = .95)))

lda.fit
summary(lda.fit)

pred_sex = predict(rf.fit, test_df)
table(pred_sex, test_df$animal.registration.sex)

pred.accuracy = round(mean(pred_sex == test_df$animal.registration.sex)*100,2)
pred.accuracy

pred_sex = predict(lda.fit, test_df)
table(pred_sex, test_df$animal.registration.sex)

pred.accuracy = round(mean(pred_sex == test_df$animal.registration.sex)*100,2)
pred.accuracy

true_sex = predict(rf.fit, mis_df)
MIS_ID
mis_df$animal.registration.sex

#' ## Correct for Sex Specific Variance: Examine Genes Driving Sex-Specific Variance

#+ Examine Genes Driving Sex-Specific Variance
################################################################################
#########     Examine Genes Driving Sex-Specific Variance      #################
################################################################################

# DE Test for Sex (Ensure males and females are ordered)
design = ~ animal.registration.sex # Primary variable needs to be last.
dds_vp <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = design)
dds_vp <- dds_vp[rowSums(counts(dds_vp))/ncol(dds_vp) >= 1,]
dds_vp <- estimateSizeFactors(dds_vp)





#dds_vp_bk<- dds_vp
#assay(dds_vp) <- MC_DEg_counts
#?DESeq

res_vp <- DESeq(dds_vp) %>%
  results(alpha = 0.05,lfcThreshold=0.5)
# Adjust the result dataframe
#res_vp <- res_vp_bk
res_vp <- res_vp %>%
  as.data.frame() %>%
  mutate(ENSEMBL_RAT = row.names(res_vp)) %>%
  mutate(SIG = ifelse((abs(log2FoldChange) >= 0.5 & padj <= 0.05), 'SIG','NS')) %>%
  mutate(CHROM = factor(case_when(ENSEMBL_RAT %in% Y_ens_id ~ "Y",
                           ENSEMBL_RAT %in% X_ens_id ~ "X",
                           ENSEMBL_RAT %!in% c(X_ens_id,Y_ens_id) ~ "AUTO"),
                        levels = c("X", "Y", "AUTO"))) %>%
  arrange(padj) %>%
  as_tibble()
res_vp 

# Generate a volcano plot
res_vp %>%  
  ggplot(aes(log2FoldChange, -log10(padj), color = CHROM, size = CHROM)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = -0.5, linetype="dashed") +
  ggtitle("DE Genes from Females to Males:
Adjusted p-value <= 0.05; |Log2FC| >= 0.5") +
  ylab(bquote('-'~Log[10]~ 'Adjusted p-value (Bonferroni)')) +
  scale_color_manual(values=c('blue', 'red', 'grey')) +
  scale_size_manual(values=c(3, 3, 1)) +
  theme_bw()

# Generate an accompanying pvalue histogram
res_vp %>%
  ggplot(aes(x = padj, fill = CHROM)) +
  #geom_freqpoly() +
  ggtitle("DE Genes Females to Males (adjusted p-value <= 0.05)") +
  geom_dotplot(aes(x = padj), dotsize = 0.5, alpha = 0.5) +
  geom_vline(xintercept = 0.05, linetype="dashed") +
  xlab('Adjusted p-value (Bonferroni)') +
  ylab("Density") +
  scale_fill_manual(values=c('blue', 'red', 'grey')) +
  theme_bw()

################################################################################
###### Manual Annotation Opportunity: Fill in VAR_SEX_CRT variable #############
################################################################################
# Determine that Strategy for Sex-Specific Variance Adjustment

#' ### Adjust for Between Sex Variance

#+ Adjust for Between Sex Variance
################################################################################
########### Adjust for Between Sex Variance  ###################################
################################################################################
# TODO: Generate a decision tree for which Sex-Specific Strategy to utilize
# TODO: Adjust rld.test to just rld
# Key:
# Remove - Rm
# Median center - Mc
# Model - Mo
# Sex Genes - Sg
# Autosomal Genes - Au
# DE Genes (Sex) - DEg
# All Genes - Ag

VAR_SEX_CRT <- "MC_DEg"

# Collect male and female samples
M_samples <- col_data %>%
  filter(Tissue == TISSUE) %>%
  filter(animal.registration.sex == 'Male') %>%
  dplyr::select(sample_key) %>% unlist() %>% as.character()
F_samples <- col_data %>%
  filter(Tissue == TISSUE) %>%
  filter(animal.registration.sex == 'Female') %>%
  dplyr::select(sample_key) %>% unlist() %>% as.character()
SIG_genes <- res_vp %>%
  filter(SIG == "SIG") %>%
  dplyr::select(ENSEMBL_RAT) %>% unlist() %>% as.character()
NS_genes <- res_vp %>%
  filter(SIG == "NS") %>%
  dplyr::select(ENSEMBL_RAT) %>% unlist() %>% as.character()
# Select the counts
MSIG_counts <- assay(rld[SIG_genes, M_samples])
FSIG_counts <- assay(rld[SIG_genes, F_samples])
MNS_counts <- assay(rld[NS_genes, M_samples])
FNS_counts <- assay(rld[NS_genes, F_samples])
  
if(VAR_SEX_CRT == "MC_DEg"){
  # Median Center Significant Genes
  # Collects median of each row, then subtracts by row medians
  M_medians <- apply(MSIG_counts,1,median)
  MSIG_counts <- MSIG_counts - M_medians
  F_medians <- apply(FSIG_counts,1,median)
  FSIG_counts <- FSIG_counts - F_medians
  SIG_counts <- cbind(MSIG_counts, FSIG_counts)
  
  # Non-Significant Genes
  NS_counts <- cbind(MNS_counts, FNS_counts)
  MC_DEg_counts <- rbind(SIG_counts, NS_counts)
  MC_DEg_counts <- MC_DEg_counts[, colnames(assay(rld))]
  #rld.test <- rld
  assay(rld.test) <- MC_DEg_counts
}else if(VAR_SEX_CRT == "MC_Ag"){
  # Median Center Significant Genes
  # Collects median of each row, then subtracts by row medians
  M_medians <- apply(MSIG_counts,1,median)
  MSIG_counts <- MSIG_counts - M_medians
  F_medians <- apply(FSIG_counts,1,median)
  FSIG_counts <- FSIG_counts - F_medians
  SIG_counts <- cbind(MSIG_counts, FSIG_counts)
  
  # Non-Significant Genes
  M_medians <- apply(MNS_counts,1,median)
  MNS_counts <- MNS_counts - M_medians
  F_medians <- apply(FNS_counts,1,median)
  FNS_counts <- FNS_counts - F_medians
  NS_counts <- cbind(MNS_counts, FNS_counts)
  MC_Ag_counts <- rbind(SIG_counts, NS_counts)
  MC_Ag_counts <- MC_Ag_counts[, colnames(assay(rld))]
  #rld.test <- rld
  assay(rld.test) <- MC_Ag_counts
}
  
# Find variables associated with PCs
pc_cor_df <- cor_PC_1_4(rld = rld.test, ntop = 500, intgroups = names(colData(rld.test))) %>%
  filter(Adjusted_R_Sq > 0.1) %>%
  arrange(PC)
# Examine Variables of interest for Pcs 1 and 2
pc_cor_df %>%
  filter(PC %in% c(1,2))
# Visualize variables in PC's 1 and 2
voi <- pc_cor_df %>%
  filter(PC %in% c(1,2)) %>%
  dplyr::select(Condition) %>% unlist() %>% as.character() %>% unique()
# Plot the top variables associated with PCs 1 or 2
for( v in voi){
  # Quick investigation of variables
  plot(DESeq2::plotPCA(rld.test, intgroup =v) +
         guides(color=guide_legend(title=v)))
}
#' #### We see just how well duplicate samples correlate regardless of sequencing batch
DESeq2::plotPCA(rld.test, intgroup ="animal.registration.sex") +
  guides(color=guide_legend(title="Sex"))

mypar()
pcaData <- DESeq2::plotPCA(rld.test, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key", unique(voi)), 
                           returnData=TRUE, ntop = 10000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
  geom_point(size=3) +
  geom_label_repel(aes(label=RNA_extr_conc),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\ny ~ sex + cohort")) +
  guides(color=guide_legend(title="animal.key.anirandgroup")) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())

#' ## Outlier Detection (By Sample)

#+ Outlier Detection (By Sample)
################################################################################
#####     Outlier Detection (By Sample)     ####################################
################################################################################
pri_var <- "animal.registration.sex"
mypar()
pcaData <- DESeq2::plotPCA(rld.test, 
                           intgroup=c(pri_var, "sample_key",
                                      "animal.key.anirandgroup"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- pcaData %>%
  ggplot(aes_string("PC1", "PC2", color=pri_var)) +
  geom_point(size=3) +
  #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\nNaive Model (~ 1)")) +
  guides(color=guide_legend(title=pri_var)) +
  theme(legend.title=element_blank()) +
  #scale_color_manual(values=ec_colors) +
  scale_shape_manual(values=c(3, 19))
plot(p)

# (iii) standardization of the data so that the expression measures for each array have mean 0 and variance 1 across genes.
edata <- t(assay(rld))
rv <- rowVars(assay(rld))
top_genes <- order(rv, decreasing = TRUE)[seq_len(min(100, length(rv)))]
edata <- scale(t(assay(rld)[top_genes, ]))
dim(edata)

# check that we get mean of 0 and sd of 1
apply(edata , 2, mean)
apply(edata , 2, sd)
edata.mad=apply(edata,2,mad)

# Sometimes genes will demonstrate Median Absolute Deviation (mad) == 0, this causes pcout() to stop. We have manually examined these genes, remove them.
remove_genes <- names(edata.mad[edata.mad == 0])
edata <- edata[,colnames(edata)[colnames(edata) %!in% remove_genes]]

# identify multivariate outliers
x.out <- pcout(edata, makeplot = T)

# Create a dataframe with 3 columns: gene, weight_combined, gene_index
df_out <- data.frame(names(x.out$wfinal), x.out$wfinal) %>% as_tibble()
names(df_out) <- c("sample_key", 'weight_combined')
df_out$index <- row.names(df_out)
# Add all other annotations
df_out <- left_join(df_out, as_tibble(colData(rld)), by = c("sample_key"))

# Arrange the dataframe by weight
df_out <- df_out %>%
  arrange(weight_combined)

# Plot the Outliers by combined weight (annotation 1)
df_out %>%
  mutate(plot_label = ifelse(sample_key %in% 
                               c('90009017003_SF1','90039017003_SF1',
                                 '90027017003_SF1','90047017003_SF1'),
                             sample_key, '')) %>%
  ggplot(aes(index, weight_combined, color = animal.key.anirandgroup)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_label_repel(aes(label=plot_label),hjust=0, vjust=0) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  guides(color=guide_legend(title="animal.key.anirandgroup")) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())

# Number of strong outlier samples
samples_out <- df_out %>%
  filter(weight_combined < 0.05) %>%
  dplyr::select(sample_key) %>% unlist() %>% as.character()
# Adjust the annotation of the rld object
colData(rld)$sample_outlier <- factor(ifelse(colData(rld)$sample_key %in% samples_out, 
                                             'Outlier', 'Normal'),
                                      levels = c("Normal", 'Outlier'))

# Determine what is correlated with outliers
# Function to grab topp 500 variance genes and find correlated variables
###########
cor_outlier <- function(rld = rld, ntop = 500, intgroups = names(colData(rld))){
  # Create out dataframe
  out_df <- data.frame(matrix(ncol = 2, nrow = 0))
  names(out_df) <- c("Adjusted_R_Sq","Condition")
  for(intgroup in intgroups) {
    #intgroup <- "animal.key.anirandgroup"
    group <- colData(rld)[[intgroup]]
    d <- colData(rld) %>% as_tibble()
    d <- d %>% mutate(sample_outlier = ifelse(sample_outlier == 'Outlier', 1, 0))
    # Perform linear model on PC vs
    if(length(unique(group[!is.na(group)])) >=2){
      int_r2 <- (lm(d[["sample_outlier"]] ~ d[[intgroup]]) %>% summary())$adj.r.squared
      int_df <- data.frame(Adjusted_R_Sq = int_r2,
                           Condition = intgroup) %>% as_tibble()
      out_df <- rbind(out_df, int_df)
    }else{
      NA
    }
  }
  out_df %>%
    arrange(desc(Adjusted_R_Sq)) %>%
    filter(Adjusted_R_Sq >= 0.1) %>%
    filter(Condition  != "sample_outlier")
}
###########
# Examine metadata varibales most correlated with outlier samples to determine if outliers are genuine
cor_out_df <- cor_outlier(rld = rld, ntop = 500, intgroups = names(colData(rld)))
cor_out_df <- cor_out_df %>% 
  filter(Condition != "outlier")

# PCA that examines all variables associated with outlier status
voi <- c(as.character(cor_out_df$Condition), "animal.key.anirandgroup")
for( v in voi){
  mypar()
  pcaData <- DESeq2::plotPCA(rld, 
                             intgroup=c(pri_var, "sample_key",
                                        "animal.key.anirandgroup","sample_outlier",v), 
                             returnData=TRUE, ntop = 500)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  p <- pcaData %>%
    ggplot(aes_string("PC1", "PC2", color=v, shape = "sample_outlier")) +
    geom_point(size=3) +
    geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\nNaive Model (~ 1)")) +
    guides(color=guide_legend(title=pri_var)) +
    theme(legend.title=element_blank()) +
    #scale_color_manual(values=ec_colors) +
    scale_shape_manual(values=c(3, 19))
  plot(p)
}

################################################################################
###### Manual Annotation Opportunity: Fill in OUTLIER variable #################
################################################################################
# Note: If Samples Labled as Outliers are in Acute exercise groups, then they likely do not meet our criteria for outliers.


