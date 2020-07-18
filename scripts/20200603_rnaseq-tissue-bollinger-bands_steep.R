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

# Global options
options(dplyr.print_max = 100)

# Source the functions
source(paste0(WD,'/functions/not_in.R'))
source(paste0(WD,'/functions/rat_mouse_ortho.R'))
source(paste0(WD,'/functions/mouse2rat_ortho.R'))
source(paste0(WD,'/functions/lmp.R'))
source(paste0(WD,'/functions/cor_PC_1_6.R'))
source(paste0(WD,'/functions/elbow_finder.R'))
source(paste0(WD,'/functions/cor_outlier2.R'))
source(paste0(WD,'/functions/sin.R'))
source(paste0(WD,'/functions/cos.R'))

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

models_df <- data.frame()

#TISSUE <- "Kidney"
for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
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
  # Collect the Formula
  FORMULA <- df_tbl %>%
    filter(Tissue == TISSUE) %>%
    select(Formula) %>% unique() %>% 
    unlist() %>% as.character() %>% as.formula()


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
#design = ~1 # Primary variable needs to be last.
design = FORMULA
( title = paste0('Design: ',as.character(design)) )
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

#' ### Adjust Variance

#+ Adjust Variance
################################################################################
########### Adjust Variance  #######################################
################################################################################
if(ADJ_VAR != 'None'){
  for(adj_var in ADJ_VAR){
    # Duplicate the rld object
    rld_final <- rld
    # The batch effect can only be removed with limma
    # https://support.bioconductor.org/p/76099/ (See Michael Love's Comment)
    assay(rld_final) <- limma::removeBatchEffect(assay(rld), rld[[adj_var]])
    
    # Examine the primary variable of interest to see if we've solved our issue
    # Before:
    p <- DESeq2::plotPCA(rld, intgroup =adj_var) +
      guides(color=guide_legend(title=adj_var))
    plot(p)
    # After
    p <- DESeq2::plotPCA(rld_final, intgroup = adj_var, ntop = 500) +
      guides(color=guide_legend(title=adj_var))
    plot(p)
    rld <- rld_final
  }
}

#### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 20000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf(paste0(WD,"/plots/20200426_rnaseq-",TIS,"-PCA-sexmod-modeling_steep.pdf"),
# width = 6, height = 4)
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
  geom_point(size=3) +
  #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  ggtitle(paste0("PCA of ",TISSUE)) +
  guides(color=guide_legend(title="animal.key.anirandgroup")) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())
#dev.off()

#' #### Annotate Data for Modeling By Cluster

#+ Annotate Data for Modeling By Cluster
################################################################################
################ Annotate Data for Modeling By Cluster #########################
################################################################################
time_cols <- tod_cols %>%
        filter(sample_key %in% nona_sams)

# Select the normalized counts
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
# by_gene_df <- by_gene_df %>%
#         mutate(CIRC = ifelse(ENSEMBL_RAT %in% circ_df$ENSEMBL_RAT, 'CIRC', 'NON-CIRC'))
# Add the gene symbol
by_gene_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")

# Join the dataframes and nest
by_gene_df7 <- tod_cols %>%
        left_join(t_counts, by = "sample_key") %>%
        filter(animal.key.anirandgroup %in% c('Control - 7 hr')) %>%
        group_by(ENSEMBL_RAT) %>%
        arrange(sample_key) %>%
        nest()
# Add the gene symbol
#by_gene_df7$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df7$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")

# Cleanup
#rm(count_data, col_data, counts_centered, dds ,dds1,dds2,F_centered,F_counts,
#   M_centered,M_counts,pcaData,t_counts)
# TODO: The entire dataframe cannot be unnested. It may be a memory issue or it may be a corrupted row -- come back and troubleshoot. For now, randomly select 2K rows to continue analysis.
#by_gene_df2 <- by_gene_df
#by_gene_df72 <- by_gene_df7
#by_gene_df <- by_gene_df2
#by_gene_df7 <- by_gene_df72

# randomRows = sample(1:nrow(by_gene_df[,1]), 10, replace=F)
# by_gene_df <- by_gene_df[randomRows, ]
# by_gene_df7 <- by_gene_df7[randomRows, ]

# Must be true
all(by_gene_df$ENSEMBL_RAT == by_gene_df7$ENSEMBL_RAT)

# Generate model functions for the dataframes
gam_mod <- function(df) {
        lm(count ~ ns(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
}
# Generate model functions for the residual
gam_mod2 <- function(df) {
  lm(count ~ ns(resid, df = 4), data = df)
}
# Generate a model function for the dataframes
sin_mod <- function(df) {
        lm(count ~ SIN(specimen.collection.t_death_hour) + 
                   COS(specimen.collection.t_death_hour),
           data = df)
}
# Generate a model function for the residual
sin_mod2 <- function(df) {
  lm(count ~ SIN(resid) + 
       COS(resid),
     data = df)
}

# Generate a model that combines circadian with exercise (circadian first)
ce_mod <- function(df) {
  lm(count ~ SIN(specimen.collection.t_death_hour) + 
       COS(specimen.collection.t_death_hour) +
       ns(specimen.collection.t_exercise_hour_sqrt, df = 4) + 1, data = df)
}
# Generate a model that combines circadian with exercise (exercise first)
ec_mod <- function(df) {
  lm(count ~ SIN(specimen.collection.t_death_hour) + 
       COS(specimen.collection.t_death_hour) +
       ns(specimen.collection.t_exercise_hour_sqrt, df = 4) + 1, data = df)
}

# Add the gene symbol
by_gene_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df7$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")

# In case you'd like to subset the data
#by_gene_df_bk <- by_gene_df
by_gene_df <- by_gene_df_bk

by_gene_df <- by_gene_df %>%
  filter(SYMBOL_RAT == 'Arntl')

# Perform a series of analyses for each model
for(mdl in c('gam', 'sin', 'ce', 'ec')){
  # circadian with exercise Model (circadian first)
  ################################################################################
  # Collect variables
  ( model_col <- as.symbol(paste0(mdl,'_model')) )
  anova_col <- as.symbol(paste0(mdl,'_anova'))
  resid_col <- as.symbol(paste0(mdl,'_resid'))
  metrics_col <- as.symbol(paste0(mdl,'_metrics'))
  summary_col <- as.symbol(paste0(mdl,'_summary'))
  MODEL <- match.fun(paste0(mdl,'_mod'))
  
  # Run models and save as a column
  by_gene_df <- by_gene_df %>%
    mutate(!!model_col := map(data, MODEL))
  # Examine the ANOVA report on models
  by_gene_df <- by_gene_df %>%
    mutate(!!anova_col := map(!!model_col, anova))
  # Add the residuals
  by_gene_df <- by_gene_df %>%
    mutate(!!resid_col := map2(data, !!model_col, modelr::add_residuals))
  # Examine the model metrics
  by_gene_df <- by_gene_df %>%
    mutate(!!metrics_col := map(!!model_col, broom::glance))
  # Examine some model summaries
  by_gene_df <- by_gene_df %>%
    mutate(!!summary_col := map(!!model_col, summary))
}

by_gene_df %>%
  ungroup() %>%
  filter(SYMBOL_RAT == 'Arntl') %>%
  select(ce_model) %>%
  mutate(anova_rep = map(ce_model, .f = anova)) %>%
  unnest(anova_rep)

by_gene_df %>%
  ungroup() %>%
  filter(SYMBOL_RAT == 'Arntl') %>%
  select(ec_model) %>%
  mutate(anova_rep = map(ec_model, .f = anova)) %>%
  unnest(anova_rep)

by_gene_df %>%
  ungroup() %>%
  filter(SYMBOL_RAT == 'Arntl') %>%
  select(gam_model) %>%
  mutate(anova_rep = map(gam_model, .f = anova)) %>%
  unnest(anova_rep)

by_gene_df %>%
  ungroup() %>%
  filter(SYMBOL_RAT == 'Arntl') %>%
  select(sin_model) %>%
  mutate(anova_rep = map(sin_model, .f = anova)) %>%
  unnest(anova_rep)

models_df %>%
  filter(TISSUE == 'Kidney') %>%
  filter(ENSEMBL_RAT == 'ENSRNOG00000014448')

models_df$TISSUE %>% table()




# circadian with exercise Model (circadian first)
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
  mutate(ce_model = map(data, ce_mod))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
  mutate(ce_ANOVA = map(ce_model, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
  mutate(ce_resid = map2(data, ce_model, modelr::add_residuals))
# Examine the model metrics
by_gene_df <- by_gene_df %>%
  mutate(ce_metrics = map(ce_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
  mutate(ce_summary = map(ce_model, summary))
# Save the model metrics
ce_metrics <- by_gene_df %>%
  unnest(ce_metrics)
# Object with the anova metrics
ce_anova <- by_gene_df %>%
  unnest(ce_ANOVA)

# exercise with circadian Model (exercise first)
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
  mutate(ec_model = map(data, ec_mod))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
  mutate(ec_ANOVA = map(ec_model, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
  mutate(ec_resid = map2(data, ec_model, modelr::add_residuals))
# Examine the model metrics
by_gene_df <- by_gene_df %>%
  mutate(ec_metrics = map(ec_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
  mutate(ec_summary = map(ec_model, summary))
# Save the model metrics
ec_metrics <- by_gene_df %>%
  unnest(ec_metrics)
# Object with the anova metrics
ec_anova <- by_gene_df %>%
  unnest(ec_ANOVA)

# Generalized Additive Models (pass1)
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
        mutate(gam_model = map(data, gam_mod))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
        mutate(gam_ANOVA = map(gam_model, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
        mutate(gam_resid = map2(data, gam_model, modelr::add_residuals))
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
        mutate(sin_resid = map2(data, sin_model, modelr::add_residuals))
# Examine the model metrics
by_gene_df <- by_gene_df %>%
        mutate(sin_metrics = map(sin_model, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
        mutate(sin_summary = map(sin_model, summary))
# Examine the model metrics
sin_metrics <- by_gene_df %>%
  mutate(sin_metrics = map(sin_model, broom::glance)) %>%
  unnest(sin_metrics)

# Generalized Additive Model (on SIN Model residuals)
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
  mutate(gam_model2 = map(sin_resid, gam_mod2))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
  mutate(gam_ANOVA2 = map(gam_model2, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
  mutate(gam_resid2 = map2(sin_resid, gam_model2, modelr::add_residuals))
# Examine the model metrics
by_gene_df <- by_gene_df %>%
  mutate(gam_metrics2 = map(gam_model2, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
  mutate(gam_summary2 = map(gam_model2, summary))
# Save the model metrics
gam_metrics2 <- by_gene_df %>%
  mutate(gam_metrics2 = map(gam_model2, broom::glance)) %>%
  unnest(gam_metrics2)

# SIN/COS Model (on SIN Model residuals)
################################################################################
# Run models and save as a column
by_gene_df <- by_gene_df %>%
  mutate(sin_model2 = map(gam_resid, sin_mod2))
# Examine the ANOVA report on models
by_gene_df <- by_gene_df %>%
  mutate(sin_ANOVA2 = map(sin_model2, anova))
# Add the residuals
by_gene_df <- by_gene_df %>%
  mutate(sin_resid2 = map2(gam_resid, sin_model2, modelr::add_residuals))
# Examine the model metrics
by_gene_df <- by_gene_df %>%
  mutate(sin_metrics2 = map(sin_model2, broom::glance))
# Examine some model summaries
by_gene_df <- by_gene_df %>%
  mutate(sin_summary2 = map(sin_model2, summary))
# Examine the model metrics
sin_metrics2 <- by_gene_df %>%
  mutate(sin_metrics2 = map(sin_model2, broom::glance)) %>%
  unnest(sin_metrics2)

# Generate the data frame
modelsr2_df <- data.frame(TISSUE = TISSUE,
           ENSEMBL_RAT = (sin_metrics %>% 
                             arrange(factor(ENSEMBL_RAT, levels = by_gene_df$ENSEMBL_RAT)))$ENSEMBL_RAT,
           SIN1_R2 = (sin_metrics %>% 
                            arrange(factor(ENSEMBL_RAT, levels = by_gene_df$ENSEMBL_RAT)))$adj.r.squared,
           SIN2_R2 = (sin_metrics2 %>% 
                        arrange(factor(ENSEMBL_RAT, levels = by_gene_df$ENSEMBL_RAT)))$adj.r.squared,
           GAM1_R2 = (gam_metrics %>% 
                        arrange(factor(ENSEMBL_RAT, levels = by_gene_df$ENSEMBL_RAT)))$adj.r.squared,
           GAM2_R2 = (gam_metrics2 %>% 
                        arrange(factor(ENSEMBL_RAT, levels = by_gene_df$ENSEMBL_RAT)))$adj.r.squared) %>%
  as_tibble()

#####################

# Concatenate the dataframes
models_df <- rbind(models_df, modelsr2_df)
# Save the final output table
models_file <- paste0(WD,'/data/20200603_rnaseq-tissue-models-residuals-r2-table_steep.txt')
#write.table(models_df, file = models_file,sep = '\t',row.names = F,quote = F)
#}
# Load the file
models_df <- read.table(file = models_file ,sep = '\t', header = T, check.names = F) %>% as_tibble()




# Join this table with DE table
C0_C7_E7_file <- paste0(WD,'/data/20200624_rnaseq-tissue-C0-C7-E7-table_steep.txt')
C0_C7_E7_df <- read.table(file = C0_C7_E7_file ,sep = '\t', header = T, check.names = F) %>% as_tibble()

# Join the tables
bin_df <- full_join(models_df, C0_C7_E7_df, by = c("TISSUE","ENSEMBL_RAT"))
bin_df <- bin_df %>%
  mutate(BIN_TYPE = case_when(SIN1_R2 > GAM1_R2 + 0.15 ~ 'C_STRONG',
                               GAM1_R2 > SIN1_R2 + 0.15 ~ 'E_STRONG',
                               ((SIN1_R2 <= GAM1_R2 + 0.15) & (GAM1_R2 <= SIN1_R2 + 0.15)) ~ 'AMB'))

library("PerformanceAnalytics")
bin_df %>% 
  sample_n(5000) %>%
  select(SIN1_R2, GAM1_R2, SIN2_R2, GAM2_R2) %>%
  chart.Correlation(histogram=TRUE, pch=19)

bin_df %>% 
  sample_n(5000) %>%
  select(GAM1_R2, SIN1_R2, 
         C0C7_padj, C0E7_padj, C7E7_padj) %>%
  chart.Correlation(histogram=TRUE, pch=19)


bin_df <- bin_df %>%
  mutate(C0C7_SIG = 
           ifelse(C0C7_padj <= 0.05 & abs(C0C7_log2FoldChange) >= 0.25,
         'SIG','NOT_SIG')) %>%
  mutate(C0E7_SIG = 
           ifelse(C7E7_padj <= 0.05 & abs(C7E7_log2FoldChange) >= 0.25,
                  'SIG','NOT_SIG')) %>%
  mutate(C7E7_SIG = 
           ifelse(C7E7_padj <= 0.05 & abs(C7E7_log2FoldChange) >= 0.25,
                  'SIG','NOT_SIG'))

bin_df %>% 
  sample_n(5000) %>%
  ggplot(aes(y = SIN1_R2, x = GAM1_R2, color = C7E7_SIG)) +
  geom_point(alpha = 0.5) +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
  xlab("R^2 Natural Spline Model (Exercise)") +
  ylab("R^2 SIN/COS Model (Circadian)") +
  ggtitle("R2 Comparisons Between Models") +
  scale_color_manual(values=c("grey", "red", "grey"))






names(bin_df)
set.seed(100)
p <- bin_df %>% 
  filter(TISSUE == 'Kidney') %>%
  #sample_n(5000) %>%
  ggplot(aes(y = SIN1_R2, x = GAM1_R2, color = BIN_TYPE)) +
  geom_point(alpha = 0.5) +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
  xlab("R^2 Natural Spline Model (Exercise)") +
  ylab("R^2 SIN/COS Model (Circadian)") +
  ggtitle("R2 Comparisons Between Models")
plot(p)
p <- bin_df %>% 
  filter(TISSUE == 'Kidney') %>%
  #sample_n(5000) %>%
  ggplot(aes(y = SIN1_R2, x = GAM1_R2, color = abs(C0C7_padj))) +
  geom_point(alpha = 0.5) +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
  xlab("R^2 Natural Spline Model (Exercise)") +
  ylab("R^2 SIN/COS Model (Circadian)") +
  ggtitle("R2 Comparisons Between Models") +
  scale_color_continuous(high = "red", low = "black")
plot(p)
bin_df %>% 
  filter(TISSUE == 'Kidney') %>%
  sample_n(10000) %>%
  ggplot(aes(y = SIN1_R2, x = GAM1_R2, color = abs(C7E7_log2FoldChange))) +
  geom_point(alpha = 0.5) +
  xlim(0,1) + ylim(0,1) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
  geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
  xlab("R^2 Natural Spline Model (Exercise)") +
  ylab("R^2 SIN/COS Model (Circadian)") +
  ggtitle("R2 Comparisons Between Models") +
  scale_color_continuous(high = "red", low = "grey")


#

 # # Load dependencies
 # pacs...man <- c("ggplot2","vcd","reshape2","ggthemes","PerformanceAnalytics","GGally", "DESeq2","org.Gg.eg.db","stringr","VGAM","tidyverse","janitor","elasticnet","earth","gbm","rpart","randomForest","Cubist","party","doMC","lattice","AppliedPredictiveModeling","ipred","caret")
 # lapply(pacs...man, FUN = function(X) {
 #   do.call("library", list(X)) 
 # })
# pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","BSgenome.Ggallus.UCSC.galGal4","BSgenome.Ggallus.UCSC.galGal5","rtracklayer", "caret","pROC","modelr","ggplot2","e1071","doMC","glmnet", "lattice","MASS","pamr","pls","sparseLDA","lubridate","reshape2","kernlab", "klaR","latticeExtra","earth","partykit","gtools")
# lapply(pacs...man, FUN = function(X) {
#   do.call("library", list(X)) 
# })


################################################################################
############## Create Dummy Variables for Categorical Dependents ###############
################################################################################

#TISSUE
dmy <- dummyVars(" ~ TISSUE ", data = bin_df)
dum_df <- as_tibble(predict(dmy, newdata = bin_df))
# Bind back into main dataframe
bin_df <- cbind(bin_df, dum_df)

################################################################################
############## Choose the Dependent and Response Variables #####################
################################################################################

dim(bin_df)
names(bin_df)

bin_df <- bin_df %>%
  filter(TISSUE == 'Kidney')
# Determine which columns should be considered in the model (response and dependents)
response <- 'BIN_TYPE'
dependents_ungrouped <- c('SIN1_R2','GAM1_R2','SIN2_R2','GAM2_R2',
                          'C0C7_padj','C0C7_log2FoldChange',
                          'C0E7_padj','C0E7_log2FoldChange',
                          'C7E7_padj','C7E7_log2FoldChange')

################################################################################
############## Visualize and Transform Certain Variables #######################
################################################################################

# Visualize before transformations
#####################
hist(bin_df$SIN1_R2, breaks = 100)
hist(bin_df$GAM1_R2, breaks = 100)
hist(bin_df$C0C7_log2FoldChange, breaks = 100)
hist(bin_df$C0C7_padj, breaks = 100)

# # Grouped Predictors
# #####################
# ## Use caret's preProcess function to transform for skewness
# # Center and Scale, BoxCox, Remove zero variance predictors
# pp_grouped <- preProcess(bin_df[,dependents_grouped], 
#                          method = c("center", "scale", "BoxCox", "nzv"))
# 
# ## Apply the transformations
# bin_df <- predict(pp_grouped, newdata = bin_df) %>% as_tibble()
# bin_df$BIN_TYPE <- as.factor(bin_df$BIN_TYPE)

# Ungrouped Predictors
#####################
df_ungrouped <- bin_df[,dependents_ungrouped]
## Use caret's preProcess function to transform for skewness
# Center and Scale, BoxCox, Remove zero variance predictors
# Consider nzv
pp_ungrouped <- preProcess(df_ungrouped[,dependents_ungrouped], 
                           method = c("center", "scale", "BoxCox"))
#pp_ungrouped <- preProcess(df_ungrouped[,dependents_ungrouped], 
#                           method = c("nzv"))

## Apply the transformations
df_ungrouped <- predict(pp_ungrouped, newdata = df_ungrouped) %>% as_tibble()
df_ungrouped$BIN_TYPE <- as.factor(bin_df$BIN_TYPE)

# Visualize after transformations
#####################
hist(df_ungrouped$SIN2_R2, breaks = 100)
hist(df_ungrouped$GAM2_R2, breaks = 100)
hist(df_ungrouped$C0C7_padj, breaks = 100)

# Update the dependents
#dependents_grouped <- names(bin_df)[names(bin_df) %!in% c("SNP_MA")]
dependents_ungrouped <- names(df_ungrouped)[names(df_ungrouped) %!in% c("BIN_TYPE")]

################################################################################
############## Test Predictors for High Correlation ############################
################################################################################

# # Grouped
# #####################
# predCorr <- cor(dfval)
# highCorr <- findCorrelation(predCorr, .99)
# # TODO: Clean up this command
# #bin_df <- bin_df[-highCorr]

# Ungrouped
#####################
predCorr <- cor(df_ungrouped[,dependents_ungrouped])

library("PerformanceAnalytics")

chart.Correlation(df_ungrouped[,dependents_ungrouped], histogram=TRUE, pch=19)

(highCorr <- findCorrelation(predCorr, .99))
# TODO: Clean up this command
#df_ungrouped <- df_ungrouped[-highCorr]

################################################################################
#####     Split the Data (Train/Test/Holdouts/Validation)      #################
################################################################################

# caret provides the createDataPartition() function that creates partitions based on stratified holdout sampling

# # Grouped
# ########################
# set.seed(123)
# train_rows <- createDataPartition(bin_df$SNP_MA, p = 0.70, list = FALSE)
# train_grouped <- bin_df[train_rows, ]
# test_grouped <- bin_df[-train_rows, ]

# Ungrouped
#########################
set.seed(123)
train_rows <- createDataPartition(df_ungrouped$BIN_TYPE, p = 0.70, list = FALSE) %>% as.vector()
train_ungrouped <- df_ungrouped[train_rows, ]
test_ungrouped <- df_ungrouped[-train_rows, ]

dim(train_ungrouped)
dim(test_ungrouped)

########################################################################
######### Determine Model Algorithm to Use (Automate with Caret) #######
########################################################################

#' ## Build Model Formula(s)

# Create the appropriate formula: (only important features)
#formula <- as.formula(paste(response, paste(dependents, collapse=" + "), sep=" ~ "))

#' Model Formula
#print(formula)

#' ## 
################################################################################
########### Feature Engineering Strategy 2 #####################################
################################################################################

# Random sample
train_ungrouped_backup <- train_ungrouped
sub_rows <- sample(row.names(train_ungrouped_backup), 5000)
train_ungrouped <- train_ungrouped_backup[sub_rows,]

train_mat <- train_ungrouped  %>% 
  dplyr::select(-BIN_TYPE) %>% 
  sapply(as.numeric) %>%
  as.matrix()
train_df <- train_ungrouped  %>% 
  dplyr::select(-BIN_TYPE) %>% 
  sapply(as.numeric) %>%
  as.data.frame()

response_vec <- train_ungrouped  %>% 
  dplyr::select(BIN_TYPE) %>% 
  unlist() %>% as.double()

train_mat[is.na(train_mat)] <- 0
train_df[is.na(train_df)] <- 0

################################################################################
########### Lasso ##############################################################
################################################################################
# Generate a grid for lasso
enetGrid <- expand.grid(.lambda = c(0.00, 0.01, 0.10),
                        .fraction = seq(0.05, 1, length = 20))
enetGrid <- expand.grid(.lambda = c(0.10),
                        .fraction = 0.1)
# Tune a lasso model (elastic net)
set.seed(123)
start_time <- Sys.time()
enetTune <- train(x = train_mat, 
                  y = response_vec,
                  method = "enet",
                  tuneGrid = enetGrid,
                  trControl = trainControl(method = "cv"),
                  preProc = c("center", "scale", "nzv"))
end_time <- Sys.time()
start_time - end_time
# Examine the best model
enetTune

# Visualize feature importance
lassoImp <- varImp(enetTune, scale = FALSE)
plot(lassoImp, top = 25, scales = list(y = list(cex = .95)))

################################################################################
########### MARS ###############################################################
################################################################################

marsFit <- earth(train_df, response_vec)
summary(marsFit)

# Generate a grid for MARS
marsGrid <- expand.grid(.degree = 1:2, .nprune = 1:10)
set.seed(123)
marsTuned <- train(train_df, response_vec,
                   method = "earth",
                   tuneGrid = marsGrid,
                   trControl = trainControl(method = "cv"))
x <- response_vec
y <- predict(marsTuned, train_df) %>% as.vector()

df_plot <- data.frame(x,y)

ggplot(df_plot, aes(x = x, y=y)) +
  geom_point()

# Visualize feature importance
marsImp <- varImp(marsTuned, scale = FALSE)
plot(marsImp, top = 25, scales = list(y = list(cex = .95)))

################################################################################
########### Tree-based Methods #################################################
################################################################################

################################################################################
###################### CART Model Development ##################################
################################################################################

# Generate the train control
ctrl <- trainControl(method = "repeatedcv",
                     number = 10, repeats = 10,
                     selectionFunction = "best",
                     savePredictions = TRUE,
                     classProbs = TRUE)

# Ungrouped
#############################
# Train the model
train_ungrouped <- as.data.frame(train_ungrouped)
set.seed(123)
start_time <- Sys.time()
rpart_fit_ungrouped <- train(x = train_ungrouped[,dependents_ungrouped], 
                             y = train_ungrouped[,response],
                             method = "rpart",
                             trControl = ctrl)
end_time <- Sys.time()
end_time - start_time
train_ungrouped <- as_tibble(train_ungrouped)

# # Grouped
# #############################
# # Train the model
# set.seed(123)
# start_time <- Sys.time()
# rpart_fit_grouped <- train(x = train_grouped[,dependents_grouped], 
#                            y = train_grouped[,response],
#                            method = "rpart",
#                            tuneLength = 30,
#                            trControl = ctrl)
# end_time <- Sys.time()
# end_time - start_time

# Proceed with the ungrouped analysis: model performance is identical but easier to interpret
rpart_fit <- rpart_fit_ungrouped

#Examine model
rpart_fit

# Examine the best complexity parameter
rpart_fit$bestTune

################################################################################
######################## Model Evaluation ######################################
################################################################################

# The trainging statistics of the model were manually saved to: ./data/20200203_cart-model-train-stats_steep.txt

# Predict how well the model performed on the test data
test_ungrouped$pred <- predict(rpart_fit, test_ungrouped)
test_ungrouped$pred_prob <- predict(rpart_fit, test_ungrouped, type = "prob")

# Determine how model performs best on test set
test_ungrouped$pred <- as.factor(test_ungrouped$pred)
test_ungrouped$BIN_TYPE <- as.factor(test_ungrouped$BIN_TYPE)
rpart_fit_cm <- confusionMatrix(test_ungrouped$pred,
                                test_ungrouped$BIN_TYPE, positive = "Y")

# Save the confusion matrix to a file
#cm_out <- data.frame(cbind(t(rpart_fit_cm$overall),t(rpart_fit_cm$byClass)))
#write.table(cm_out,
#            file=paste0("./data/20200203_test-confusion-matrix-cart-germline-snps_steep.txt"),
#            sep = '\t', quote = FALSE, row.names = FALSE)

# Visualize the tree (also save to file)
#pdf(paste0('./plots/',date,'_cart-model-vis-transformed_steep.pdf'),width=26,height=14)
plot(as.party(rpart_fit$finalModel))
#dev.off()

# Precision and recall
prec <- posPredValue(test_ungrouped$pred,
                     test_ungrouped$BIN_TYPE, positive = "Y")
rec <- sensitivity(test_ungrouped$pred,
                   test_ungrouped$BIN_TYPE, positive = "Y")
# The F-score (harmonic mean)
f <- (2 * prec * rec) / (prec + rec)
f

# Plot a roc curve
roc_rpart_train <- roc(rpart_fit$pred$obs, rpart_fit$pred$Y)
roc_rpart_test <- roc(test_ungrouped$BIN_TYPE, test_ungrouped$pred_prob$Y)

#pdf(paste0('./plots/',date,'_cart-roc-test_steep.pdf'),width=6,height=6)
## pty sets the aspect ratio of the plot region. Two options:
##                "s" - creates a square plotting region
##                "m" - (the default) creates a maximal plotting region
par(pty = "s") 
#plot(roc_rpart_train, col = 'blue', legacy.axes = TRUE, print.auc=TRUE)
plot(roc_rpart_test, col = 'red', legacy.axes = TRUE, print.auc = TRUE)
#auc(roc_rpart)
#dev.off()

#' Prediction Accuracy (AKA success rate):
#' accuracy = (TP + TN)/(TP+TN+FP+FN)
#' 
#' The Error Rate (proportion of incorrectly classified examples):
#' error rate = (FP + FN)/(TP+TN+FP+FN) = 1 - accuracy
#' 
#' Kappa:
#' The kappa statistic adjusts accuracy by accounting for the possibility of a correct prediction by chance alone
#'Sensitivity: 
#'"The sensitivity of a model (... true positive rate), measures
#'the proportion of positive examples that were correctly classified."
#'
#'sensitivity = TP/(TP+FN)
#'
#'Specificity: 
#'"The specificity of a model (... true negative rate), measures
#'the proportion of negative examples that were correctly classified."
#'
#'specificity = TN/(TN+FP)
#'
#'Precision:
#'The precision (also known as the positive predictive value) is defines as the proportion of positive examples that are truly positive.
#'
#'precision = TP/(TP+FP)
#'
#'Recall:
#'A measure of how complete the results are. A model with high recall captures a large portion of the positive examples (wide breadth)
#'
#'recall = sensitiveity = TP/(TP+FN)
#'
#'~Lantz, B. Machine Learning with R. (Packt, 2019)

# Examine confusion matrix of unfiltered (EDA)
# Pos Pred Value == precision

#'The F-measure (F1-score, F-score):
#'A measure of model performance that combines precision and recall into a single number (harmonic mean).
#'
#'F-score = (2 x precision x recall) / (recall + precision)
#'
#'~Lantz, B. Machine Learning with R. (Packt, 2019)















# Random Forest
######################
# Set up the random forest tuning grid
grid_rf <- expand.grid(mtry = c(1:5))
set.seed(123)
rfTune <- train(train_df, response_vec,
                method = "rf",
                ntree = 1000,
                trControl = trainControl(method = "cv"),
                tuneGrid = grid_rf)

# Examine model
rfTune

# Visualize feature importance
rfImp <- varImp(rfTune, scale = FALSE)
plot(rfImp, top = 25, scales = list(y = list(cex = .95)))

# Boosted Trees
######################
# Set up the tuning grid
gbmGrid <- expand.grid(.interaction.depth = seq(3, 7, by = 2),
                       .n.trees = seq(1000, 2000, by = 500),
                       .shrinkage = c(0.01, 0.1),
                       .n.minobsinnode = seq(1,11, by = 2))
# Tune the model
set.seed(123)
gbmTune <- train(train_df, response_vec,
                 method = "gbm",
                 tuneGrid = gbmGrid,
                 ##distribution = "gaussian",
                 trControl = trainControl(method = "cv"),
                 ## The gbm function produces too much output with verbose
                 verbose = FALSE)
# Exmaine model
gbmTune

# Visualize the tuning
plot(gbmTune, auto.key = list(columns = 4, lines = TRUE))

# Visualize feature importance
gbmImp <- varImp(gbmTune, scale = FALSE)
plot(gbmImp, top = 25, scales = list(y = list(cex = .95)))

# Cubist
######################
cbGrid <- expand.grid(committees = c(1:20, 30, 50, 75, 100), 
                      neighbors = c(0, 1, 5, 9))

set.seed(123)
cubistTune <- train(train_df, response_vec,
                    "cubist",
                    tuneGrid = cbGrid,
                    trControl = trainControl(method = "cv"))
# Examine the model
cubistTune

# Plot the tuning parameters by model performace
plot(cubistTune, auto.key = list(columns = 4, lines = TRUE))

# Exmaine variable performace
cbImp <- varImp(cubistTune, scale = FALSE)
plot(cbImp, top = 25, scales = list(y = list(cex = .95)))

################################################################################
########### Compare Models #####################################################
################################################################################

### Collect the resampling statistics across all the models
rs <- resamples(list("Elastic Net" = enetTune, 
                     "MARS" = marsTuned,
                     "CART" = rpartTune, 
                     "Boosted Tree" = gbmTune,
                     "Random Forest" = rfTune,
                     "Cubist" = cubistTune))

parallelplot(rs, metric = "RMSE")
parallelplot(rs, metric = "Rsquared")












































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
        grid <- data.frame(specimen.collection.t_exercise_hour_sqrt = 
                                   modelr::seq_range(
                                           sub_data$specimen.collection.t_exercise_hour_sqrt, n = length(sub_data$specimen.collection.t_exercise_hour_sqrt)))
        #grid$ENSEMBL_RAT <- g
        mod <- (by_gene_df %>% 
                        filter(ENSEMBL_RAT == g) %>% ungroup() %>% 
                        select(gam_model))[[1]][[1]]
        summary(mod)
        grid <- modelr::add_predictions(grid, mod, "pred") %>% as_tibble()
        names(grid)[1] <- "grid_t_exercise_hour_sqrt_jit"
        grid$specimen.collection.t_exercise_hour_sqrt <- 
                sub_data$specimen.collection.t_exercise_hour_sqrt
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
                                   modelr::seq_range(
                                           sub_data$specimen.collection.t_death_hour, 
                                           n = length(sub_data$specimen.collection.t_exercise_hour_sqrt)))
        #grid$ENSEMBL_RAT <- g
        mod <- (by_gene_df %>% 
                        filter(ENSEMBL_RAT == g) %>% ungroup() %>% 
                        select(sin_model))[[1]][[1]]
        summary(mod)
        grid <- modelr::add_predictions(grid, mod, "pred") %>% as_tibble()
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
        select(ENSEMBL_RAT, SYMBOL_RAT, specimen.collection.t_exercise_hour_sqrt, count)
# Collect the Control 1 hr data points
gam_hr1_df <- by_gene_df %>%
        filter(ENSEMBL_RAT == gam_gene) %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, data) %>%
        ungroup() %>%
        unnest(data) %>% 
        filter(animal.key.anirandgroup == "Control - IPE") %>%
        select(ENSEMBL_RAT, SYMBOL_RAT, specimen.collection.t_exercise_hour_sqrt, count)

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
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
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
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
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
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        geom_point(data = gam_hr7_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
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
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
               color = ENSEMBL_RAT) +
        geom_point(color = "red", alpha = 1) +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        geom_point(data = gam_hr7_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
                   color = "blue") +
        geom_point(data = gam_hr1_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
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
        ggplot(aes(specimen.collection.t_exercise_hour_sqrt, count), 
               color = ENSEMBL_RAT) +
        geom_point() +
        geom_line(data = gam_pred_df, 
                  aes(grid_t_exercise_hour_sqrt_jit, pred), 
                  size = 1, alpha = 0.8, color = "blue") +
        geom_point(data = gam_hr7_df, 
                   mapping = aes(specimen.collection.t_exercise_hour_sqrt, count),
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















