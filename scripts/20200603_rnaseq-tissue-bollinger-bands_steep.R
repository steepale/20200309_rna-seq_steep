'---
#' title: "PASS1A Rat Tissue: -- Bollinger Bands on Circadian Model"
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
TISSUE <- "Hypothalamus"
TIS <- "SCN"
OUTLIERS <- c('')
MIS_ID_Mw2Fc <- c('')
# Declare Outliers
if(TISSUE == "Kidney"){
        #OUTLIERS <- c('90042016803_SF1','90109015902_SN1')
        #OUTLIERS <- c('90109015902_SN1')
        OUTLIERS <- c('')
}else if(TISSUE == "Liver"){
        OUTLIERS <- c('90042016803_SF1')
}else if(TISSUE == "Hypothalamus"){
        OUTLIERS <- c('')
        MIS_ID_Mw2Fc <- c('90010015402_SF2')
}else if(TISSUE == "Aorta"){
        OUTLIERS <- c('')
}else if(TISSUE == "Gastrocnemius"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Heart"){
        OUTLIERS <- c('90052015802_SN1')
}else if(TISSUE == "Adrenal"){
        OUTLIERS <- c('ABC')
}else if(TISSUE == "Brown Adipose"){
        OUTLIERS <- c('')
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
print(MIS_ID_BY_SEX)

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
        filter(sample_key != OUTLIERS) %>%
        filter(!is.na(animal.registration.sex))
rownames(tod_cols) <- tod_cols$sample_key

# Adjust for the sex ID of misidentified samples

names(tod_cols)

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
# Examine histograms
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=calculated.variables.deathtime_after_acute)) +
        geom_histogram(bins = 68)
tod_cols %>%
        filter(animal.key.anirandgroup != 'Control - 7 hr') %>%
        ggplot(aes(x=specimen.collection.t_exercise_hour_sqrt)) +
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

#' #### Do Tissue Samples cluster by sex.
#' Grey samples represent "reference samples".
DESeq2::plotPCA(rld, intgroup ="animal.registration.sex") +
    guides(color=guide_legend(title="Sex"))

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf(paste0(WD,"/plots/20200505_rnaseq-",TIS,"-PCA-naive-modeling_steep.pdf"),
#    width = 6, height = 4)
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex)) +
    geom_point(size=3) +
    geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    #coord_fixed() +
    ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\nNaive Model (~ 1)")) +
    guides(color=guide_legend(title="animal.key.anirandgroup")) +
    theme(legend.title=element_blank())
#dev.off()


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
group <- tod_cols$animal.key.anirandgroup

#' #### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.
chryexp <- colMeans(assay(rld)[Y_genes,])

chryexp_df <- data.frame("counts" = chryexp)
chryexp_df$sample <- row.names(chryexp_df) %>% as.factor()
chryexp_df <- chryexp_df %>%
    mutate(sex = ifelse(sample %in% male_tis, 'Male', 'Female'))

#' ##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution. Indeed, we do in kidney.
mypar()
hist(chryexp_df$counts, breaks = 200, 
     main = "Histogram: \nY Chromosome Genes across Samples",
     xlab = "Total reads on Y chromosome genes per sample")
ggplot(chryexp_df, aes(counts, colour = sex)) +
    geom_freqpoly() +
    xlab("Total reads on Y chromosome genes per sample") +
    ylab("Frequency") +
    ggtitle("Frequency Polygon: \nY Chromosome Genes across Samples")
# summary(chryexp)

#' The distribution of sex by group
table(group, sex) # <- Bad idea.

#' #### Generate a heatmap of expression from 3 sets of genes:
#' * Genes from the Y chromosome
#' * The top and bottom 25 genes (50 total) associated with sex
#' * Randomly selected genes

#' ##### Males and Females demonstrate distinctly different gene expression profiles.
#' * Genes on the Y chromosome are a good predictor of sex in Kidney mRNA measures (Figure 1)
#' * Male and female samples show distinct correlation to one another (Figure 2)

# Ensure males and females are ordered
kidney_counts <- assay(rld)[,c(male_tis,female_tis,ref_kidneys)] %>% as.matrix()

# T-test of expression associated with sex
tt <- rowttests(kidney_counts,sex)

# Take genes from the Y chromosome
# Y_genes
# Take the top and bottom 25 genes associated with variable of interest (remove any genes in Y chromosome)
top <- row.names(tt[order(-tt$dm),][1:25,])
bot <- row.names(tt[order(tt$dm),][1:25,])
top_n_bot <- setdiff(c(top,bot), Y_genes)

# Randomly select 50 genes not in prior sets 
set.seed(123)
randos <- setdiff(row.names(tt[sample(seq(along=tt$dm),50),]), c(Y_genes,top_n_bot))
geneindex <- c(randos,top_n_bot,Y_genes)
# Generate the heatmap and support with a plot of a correlation matrix
mat <- kidney_counts[geneindex,]
mat <- mat -rowMeans(mat)
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
mypar(1,2)
image(t(mat),xaxt="n",yaxt="n",col=icolors)
y <- kidney_counts - rowMeans(kidney_counts)
image(1:ncol(y),1:ncol(y),cor(y),col=icolors,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
axis(2,1:ncol(y),sex,las=2)
axis(1,1:ncol(y),sex,las=2)

#' #### A naive t-test and genes with q values less than or equal to 0.05
#' ##### The left figure represents a histogram of p values from a naive t-test (all genes). We see that a significant proportion of genes correlate with sex. To investigate if these genes are located on sex chromosomes or autosomal chromosomes, we contruct a volcano plot on the right. To our suprise, again, genes on the Y chromosome do not show signifcant correlation to sex. Rather some genes on the X chromosome demonstrate significance, however, these genes do not make up the majority of significantly assocaited genes.
#' ## The variance in gene expression is not dicated by differential expression of genes on sex chromosomes.
mypar(1,2)
# Histogram of p values associated with ttest
hist(tt$p.value,main="",ylim=c(0,1300), breaks = 100)
plot(tt$dm,-log10(tt$p.value))
points(tt[X_genes,]$dm,-log10(tt[X_genes,]$p.value),col=1,pch=16)
points(tt[Y_genes,]$dm,-log10(tt[Y_genes,]$p.value),col=2,pch=16, xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("X","Y"),col=1:2,pch=16)
p <- tt$p.value
qvals <- qvalue(tt$p.value)$qvalue
index <- which(qvals<=0.05)
abline(h=-log10(max(tt$p.value[index])))



#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(paste0(WD,"/plots/20200505_rnaseq-",TIS,"-PCA-naive-modeling_steep.pdf"),
    width = 6, height = 4)
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\nNaive Model (~ 1)")) +
        guides(color=guide_legend(title="animal.key.anirandgroup")) +
        scale_color_manual(values=ec_colors) +
        theme(legend.title=element_blank())
dev.off()

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
pdf(paste0(WD,"/plots/20200426_rnaseq-",TIS,"-PCA-sexmod-modeling_steep.pdf"),
    width = 6, height = 4)
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle(paste0("PCA of ",TISSUE," Gene Expression:\ny ~ sex + cohort")) +
        guides(color=guide_legend(title="animal.key.anirandgroup")) +
        scale_color_manual(values=ec_colors) +
        theme(legend.title=element_blank())
dev.off()

#' #### Annotate Data for Modeling By Cluster

#+ Annotate Data for Modeling By Cluster
################################################################################
################ Annotate Data for Modeling By Cluster #########################
################################################################################

time_cols <- tod_cols %>%
        filter(sample_key %in% nona_sams)
time_cols %>%
        ggplot(aes(x = specimen.collection.t_death_hour, fill = animal.key.anirandgroup)) + 
        geom_bar(breaks = seq(0, 24), width = 2, colour = "grey") + 
        coord_polar(start = 0) + 
        theme_minimal() + 
        ggtitle("Death of Experimental Groups by Time of Day") + 
        scale_x_continuous("", limits = c(0, 24), 
                           breaks = seq(0, 24), 
                           labels = seq(0, 24)) +
        theme(legend.title = element_blank()) +
        scale_fill_manual(values=ec_colors[1:9]) +
        theme(axis.title.y =element_blank(),
              axis.text.y =element_blank(),
              axis.ticks.y=element_blank())


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
        mutate(CIRC = ifelse(ENSEMBL_RAT %in% circ_df$ENSEMBL_RAT, 'CIRC', 'NON-CIRC'))

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
#by_gene_df <- by_gene_df2
#by_gene_df7 <- by_gene_df72

# Load in the DE genes between C0 and C7 for appropriate tissues
in_file <- paste0(WD, '/data/20200426_rnaseq-',TIS,'-DEgenes-sexmod-C0vsC7_steep.txt')
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
by_gene_df <- by_gene_df %>%
        mutate(gam_metrics = map(gam_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
        mutate(gam_summary = map(gam_model, summary))
# Save the model metrics
gam_metrics <- by_gene_df %>%
        mutate(gam_metrics = map(gam_model, broom::glance)) %>%
        unnest(gam_metrics)


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
        unnest(sin_metrics) %>%
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
                                           sub_data$specimen.collection.t_exercise_hour_sqrt_jit, n = length(sub_data$specimen.collection.t_exercise_hour_sqrt_jit)))
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
                                           sub_data$specimen.collection.t_death_hour, 
                                           n = length(sub_data$specimen.collection.t_exercise_hour_sqrt_jit)))
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

gene_n <- 1
# Visualize the GAM model by gene
gam_pred_df <- by_gene_df %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data, CIRC, gam_model, gam_pred) %>%
        ungroup() %>%
        filter(row_number() == gene_n) %>%
        unnest(gam_pred)
gam_gene <- unique(gam_pred_df$ENSEMBL_RAT) %>% as.character()
# Visualize the SIN model by gene
sin_pred_df <- by_gene_df %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data, CIRC, sin_model, sin_pred) %>%
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
# Collect the Control 1 hr data points
gam_hr1_df <- by_gene_df %>%
        filter(ENSEMBL_RAT == gam_gene) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data) %>%
        ungroup() %>%
        unnest(data) %>% 
        filter(animal.key.anirandgroup == "Control - IPE") %>%
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
# Raw counts
gam_pred_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt_jit, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(gam_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Hours Post Acute Exercise (Transformed)") +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
# Raw counts with Line
gam_pred_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt_jit, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(gam_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Hours Post Acute Exercise (Transformed)") +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
# Raw counts with Line with control 7
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
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(gam_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Hours Post Acute Exercise (Transformed)") +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
# Raw counts with Line with control 7 and control 0
gam_pred_df %>%
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt_jit, count), 
               color = ENSEMBL_RAT) +
        geom_point(color = "red", alpha = 1) +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        geom_point(data = gam_hr7_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt_jit, count),
                   color = "blue") +
        geom_point(data = gam_hr1_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt_jit, count),
                   color = "blue") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(gam_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Hours Post Acute Exercise (Transformed)") +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
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
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
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

# Collect the Control 1 hr data points
sin_hr1_df <- by_gene_df %>%
        filter(ENSEMBL_RAT == sin_gene) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data) %>%
        ungroup() %>%
        unnest(data) %>%
        filter(animal.key.anirandgroup == "Control - IPE") %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, specimen.collection.t_death_hour, count)
# Visualize the raw counts, model predictions, and control 7 counts (x=Hours Post Exercise)
# Raw Counts
sin_pred_df %>%
        ggplot(aes(specimen.collection.t_death_hour, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(sin_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Time of Death (Hour)")
# Raw Counts + model
sin_pred_df %>%
        ggplot(aes(specimen.collection.t_death_hour, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = sin_pred_df, 
                  aes(grid_t_death_hour, pred), 
                  size = 1, alpha = 0.8, color = "orange") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(sin_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Time of Death (Hour)")
# Counts, Model, C7 Counts
sin_pred_df %>%
        ggplot(aes(specimen.collection.t_death_hour, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = sin_pred_df, 
                  aes(grid_t_death_hour, pred), 
                  size = 1, alpha = 0.8, color = "orange") +
        geom_point(data = sin_hr7_df, 
                   mapping = aes(specimen.collection.t_death_hour, count),
                   color = "red") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(sin_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Time of Death (Hour)")
# Counts, Model, C7 Counts, C1 Counts
sin_pred_df %>%
        ggplot(aes(specimen.collection.t_death_hour, count), 
               color = ENSEMBL_RAT) +
        geom_point(color = "red") +
        geom_line(data = sin_pred_df, 
                  aes(grid_t_death_hour, pred), 
                  size = 1, alpha = 0.8, color = "orange") +
        geom_point(data = sin_hr7_df, 
                   mapping = aes(specimen.collection.t_death_hour, count),
                   color = "blue") +
        geom_point(data = sin_hr1_df, 
                   mapping = aes(specimen.collection.t_death_hour, count),
                   color = "blue") +
        theme(legend.position = "none") +
        ggtitle(
                paste0("Expression of ",
                       unique(sin_pred_df$SYMBOL_RAT),
                       ":\nExercise Groups & Control IPE (",TISSUE,")")) +
        ylab("Expression (Transformed/Normalized)") + 
        xlab("Time of Death (Hour)")


# Visualize the gam model metrics
(by_gene_df %>%
                dplyr::slice(gene_n) %>%
                select(sin_model) %>%
                ungroup %>% 
                select(sin_model))[[1]][[1]] %>%
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
        select(ENSEMBL_RAT, SYMBOL_RAT, CIRC, r.squared.sin, p.value.sin,
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
        geom_point(alpha = 0.4) +
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
        ggtitle("R2 Comparisons Between Models:\nDifferentially Expressed Genes (C0 -> C7)") +
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
        ggtitle("R2 Comparisons Between Models:\nDifferentially Expressed Genes (C0 -> C7)") +
        geom_point(data = top_metrics,
                   mapping = aes(r.squared.gam, r.squared.sin, color = Expression), 
                   alpha = 1)+
        geom_label_repel(data = top20_metrics,
                         mapping = aes(label=SYMBOL_RAT), alpha = 0.8, 
                         hjust=0, vjust=0)















