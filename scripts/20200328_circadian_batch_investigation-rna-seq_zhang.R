#'---
#' title: "PASS1A (Rat) RNA-Seq Data: Batch Effects"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
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

#' ## Goals of Analysis
#' * TODO: Add Step-by-step goals of analysis
#' * TODO: Vet "suspect samples" in file: GET_release1_qc_report.pdf and decide if these samples should be removed or not.
#' * Examine if data demonstrate technical batch effects:
#'     * Sequencing batch (sequence date and location)
#' * Examine if data demonstrate biological batch effects:
#'     * Influence from time of day and circadian rhythm
#'     * Influence from feeding schedules
#'     
#' 
#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/SanDisk64GB/MoTrPAC/'
setwd(WD)

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2",
               "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr",
               "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather",
               "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel",
               "AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("PROPER")
BiocManager::install("org.Rn.eg.db")
mypackages = c("tidyverse","GenomicRanges", "DESeq2","devtools","ggplot2","dplyr","PROPER", "lubridate", "hms",
               "org.Rn.eg.db", "AnnotationDbi")
lapply(mypackages, FUN = function(X) {
  do.call("library", list(X)) })

############################################################
##### Functions ############################################
############################################################

# Set select
select <- dplyr::select

# Capture the Date and Author
################################################################################
date <- format.Date( Sys.Date(), '%Y%m%d' )
auth <- "Jiayu"
################################################################################

#' ## Load & Clean Data
#' ##### Data files to load:
#' * RNA-Seq from Mt. Sanai 
#'     * 1 batches & metadata (gastro limited)
#' * RNA-Seq from Stanford
#'     * 2 sequencing batches & metadata (gastro limited)

#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Mt Sanai RNASeq
###################################

# Sanai batch 1 read counts (320 samples)
sanai_1 <- read.table(file = paste0(WD,'Data/rna-seq/sinai/batch_1/summary/rsem_genes_count.txt'), 
                      header = T, row.names = 1, na.strings = "NA", check.names = F)

# Sanai batch 1 meta data
sanai_meta_1 <-  read.csv(file = paste0(WD,'Data/rna-seq/sinai/batch_1/metadata/sample_metadata_20191010.csv'), 
                            header = T, na.strings = "NA", check.names = F)


# Stanford RNASeq
###################################

# Stanford batch 1 read counts (320 samples)
sford_1 <- read.table(file = paste0(WD,'Data/rna-seq/stanford/batch_1/summary/rsem_genes_count.txt'), 
                      header = T, row.names = 1, na.strings = "NA", check.names = F)
# Stanford batch 2 read counts (320 samples)
sford_2 <- read.table(file = paste0(WD,'Data/rna-seq/stanford/batch_2/summary/rsem_genes_count.txt'), 
                      header = T, row.names = 1, na.strings = "NA", check.names = F)

# Stanford batch 1 meta data
sford_meta_1 <-  read.csv(file = paste0(WD,'Data/rna-seq/stanford/batch_1/metadata/sample_metadata_20191010.csv'), 
                          header = T, na.strings = "NA", check.names = F)
# Stanford batch 2 meta data
sford_meta_2 <-  read.csv(file = paste0(WD,'Data/rna-seq/stanford/batch_2/metadata/sample_metadata_20191010.csv'), 
                          header = T, na.strings = "NA", check.names = F)

# Merged Phenotype
###################################

load(paste0(WD,"Data/phenotype/merged/merged_dmaqc_data2019-10-13.Rdata"))

# Combine the Stanford and Mt. Sanai Data
##########################################
# Count Matrixes

# Add custom annotation to column names
colnames(sanai_1) <- paste0(colnames(sanai_1),'_SN1')
colnames(sford_1) <- paste0(colnames(sford_1),'_SF1')
colnames(sford_2) <- paste0(colnames(sford_2),'_SF2')
sanai_meta_1$sample_key <- paste0(sanai_meta_1$vial_label, '_SN1')
sford_meta_1$sample_key <- paste0(sford_meta_1$vial_label, '_SF1')
sford_meta_2$sample_key <- paste0(sford_meta_2$vial_label, '_SF2')

# Combine counts
all.data.m <- cbind(sanai_1, sford_1, sford_2) %>% as.matrix()

# Convert counts to integers
mode(all.data.m) <- "integer"

# Combine metadata (convert factors to characters to prevent error)
sanai_meta_1 <- data.frame(lapply(sanai_meta_1, as.character), stringsAsFactors=FALSE)
sford_meta_1 <- data.frame(lapply(sford_meta_1, as.character), stringsAsFactors=FALSE)
sford_meta_2 <- data.frame(lapply(sford_meta_2, as.character), stringsAsFactors=FALSE)
meta <- bind_rows(sanai_meta_1, sford_meta_1, sford_meta_2)

# Adjust column objects
########################
# To factors

meta$vial_label = as.factor(meta$vial_label)
factor_cols <- c("labelid",
                 "viallabel",
                 "animal.registration.sex")
for(fc in factor_cols){
  merged_dmaqc_data[[fc]] <- as.factor(merged_dmaqc_data[[fc]])
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
               "specimen.collection.d_visit")
for(dc in date_cols){
  merged_dmaqc_data[[dc]] <- dmy(merged_dmaqc_data[[dc]])
}
# To Dates: Jan 2018 
merged_dmaqc_data <- merged_dmaqc_data %>%
  mutate(animal.registration.d_birth = dmy(paste0('1 ',animal.registration.d_birth)))

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
               "specimen.processing.t_freeze")
for(tc in time_cols){
  merged_dmaqc_data[[tc]] <- merged_dmaqc_data[[tc]] %>% as.character() %>% parse_time()
}

# To Times: 01:21
merged_dmaqc_data <- merged_dmaqc_data %>%
  mutate(acute.test.howlongshock = 
           ifelse((merged_dmaqc_data$acute.test.howlongshock != ''), paste0('00:',acute.test.howlongshock), '00:00:00'))
merged_dmaqc_data$acute.test.howlongshock <- parse_time(merged_dmaqc_data$acute.test.howlongshock)

# Adjust labels within columns (e.g. 1 and 2 to female and male)
# see ../phenotype/merged/20191013_merged-column-dictionary_steep-edits.xlsx
########################
merged_dmaqc_data <- merged_dmaqc_data %>%
  mutate(animal.registration.sex = 
           case_when(animal.registration.sex == '1' ~ 'Female',
                     animal.registration.sex == '2' ~ 'Male'))
merged_dmaqc_data$animal.registration.sex <- as.factor(merged_dmaqc_data$animal.registration.sex)

#' ### Incorporate Annotation from Phenotypic Data
#' 

################################################################################
######### Incorporate Phenotypic Annotations ###################################
################################################################################

# Replace the name of the vial label column (for sake of join)
names(merged_dmaqc_data) <- str_replace(names(merged_dmaqc_data), 'viallabel', 'vial_label')

# To investigate if the "vial_label" key is indeed unique
#meta$vial_label %>% length()
#meta$vial_label %>% unique %>% length()
#merged_dmaqc_data$vial_label %>% length()
#merged_dmaqc_data$vial_label %>% unique %>% length()

# Perform a left join to incorporate phenotypic information
status <- left_join(meta, merged_dmaqc_data, by = 'vial_label') %>% as.data.frame()
#dim(status)

# Adjust objects in columns as needed

# Create new columns
status <- status %>% 
  mutate(Seq_batch = 
           case_when((GET_site == 'Stanford' & Seq_date == '190426') ~ "Stanford_1",
                     (GET_site == 'Stanford' & Seq_date == '190703') ~ "Stanford_2",
                     (GET_site == 'MSSM' & Seq_date == '190409') ~ "MSSM_1",
                     (GET_site == 'MSSM' & Seq_date == '190626') ~ "MSSM_2",
                     (GET_site == 'MSSM' & Seq_date == '190723') ~ "MSSM_3")
  )
#table(status$Seq_batch)

# To factors
factor_cols <- c("sample_key",
                 "Seq_batch")
for(fc in factor_cols){
  status[[fc]] <- as.factor(status[[fc]])
}

# "It is absolutely critical that the columns of the count matrix and the rows of the 
# column data (information about samples) are in the same order. DESeq2 will not make 
# guesses as to which column of the count matrix belongs to which row of the column data, 
# these must be provided to DESeq2 already in consistent order." ~ Mike Love

# TODO: Come back and figure out about missing samples and annotation

# Make sure that all values in status "sample_key" are unique
#status$sample_key %>% length()
#status$sample_key %>% unique %>% length()
# Make sure that all columns in counts matrix are unique
#colnames(all.data.m) %>% length()
#colnames(all.data.m) %>% unique %>% length()

rownames(status) <- status$sample_key
#all(rownames(status) %in% colnames(all.data.m))
#all(colnames(all.data.m) %in% rownames(status))
#all(rownames(status) == colnames(all.data.m))
all.data.m <- all.data.m[, rownames(status)]


#' ### PCA Visualization of Sequencing Batches (Unsupervised)
#' 

##########################################################################
############ PCA Sequencing Batches (Unsupervised) #######################
##########################################################################

# Perform unsupervised clustering
# Build model
count_data <- all.data.m
col_data <- status # To look at all samples 
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowMeans(counts(dds)) >= 10
dds <- dds[keep,]


#' ** Summary of counts and annotation data in a DESeqDataSet
#'
dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")

#' **Here we examine all samples by batch and tissue
#' 
# first PCA with intGroup of tissue
DESeq2::plotPCA(rld, intgroup ="Tissue") +
  guides(color=guide_legend(title="Tissues"))

# subset the dataset to only Gastro, which is my tissue of interest
rld.sub <- rld[ , (rld$Tissue == "Gastrocnemius") ]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
  geom_point(size=3) +
  #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #coord_fixed() +
  ggtitle("Gastrocnemius Samples Sequenced Across Batches")

#' can see a batch variance between Stf_1 and other batches
#' 

# here I subset the genes (select by circadian rythm) and divided into 2 groups: Stf_1 vs MtSN_1+Stf_2

# I referred to KEGG pathway id to identify circadian rythm genes 
# TODO: better way to select?

count_data_sub = count_data[which(counts$path == "04710"),]


# Build model
col_data <- status # To look at all samples 
design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds_sub <- DESeqDataSetFromMatrix(countData = count_data_sub,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowMeans(counts(dds_sub)) >= 10
dds_sub <- dds_sub[keep,]


#' ** Summary of counts and annotation data in a DESeqDataSet
#'
dds_sub <- estimateSizeFactors(dds_sub)
rs <- rowSums(counts(dds_sub))
# Normalize the counts
rld <- varianceStabilizingTransformation(dds_sub) #vst

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts_sub <- as.data.frame(norm_counts)

# Annotate normalized counts
counts_sub$ensembl <- rownames(counts_sub)
counts_sub$symbol <- mapIds(org.Rn.eg.db, counts_sub$ensembl, "SYMBOL", "ENSEMBL")
counts_sub$entrez <- mapIds(org.Rn.eg.db, counts_sub$ensembl, "ENTREZID", "ENSEMBL")
counts_sub$genename <- mapIds(org.Rn.eg.db, counts_sub$ensembl, "GENENAME", "ENSEMBL")
counts_sub$go <- mapIds(org.Rn.eg.db, counts_sub$ensembl, "GO", "ENSEMBL")
counts_sub$path <- mapIds(org.Rn.eg.db, counts_sub$ensembl, "PATH", "ENSEMBL")

rld_gastro = rld[, rld$Tissue == "Gastrocnemius"]
DESeq2::plotPCA(rld_gastro, intgroup ="Seq_batch") +
  guides(color=guide_legend(title="batch"))

# Variables of Interest
#rld_stf1 = rld.sub[, rld.sub$Seq_batch == "Stanford_1"]
#rld_ref = rld.sub[, rld.sub$Seq_batch != "Stanford_1"]

