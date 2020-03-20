#'---
#' title: "Examination of Rat RNA-Seq Data: Examine Batch Effects"
#' author: "Alec Steep and Jiayu Zhang"
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     github_document: default
#'     pdf_document: default
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' # Goals of Analysis
#' * TODO: Add Step-by-step goals of analysis
#' 
#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("org.Rn.eg.db")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2",
                "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr",
                "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather",
                "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################

# Set select
select <- dplyr::select

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and AUthor
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

#' ## Load & Clean Data
#' ##### Data Files to Load:
#' * RNA-Seq from Mt. Sanai
#'     * 3 batches
#'     * Metadata
#' * RNA-Seq from Stanford
#'     * 2 Batches
#'     * Metadata

#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Mt Sanai RNASeq
###################################

# Sanai batch 1 read counts (320 samples)
sanai_1 <- read.table(file = paste0(WD,'/sinai/batch_1/summary/featureCounts.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)
# Sanai batch 2 read counts (320 samples)
sanai_2 <- read.table(file = paste0(WD,'/sinai/batch_2/summary/featureCounts.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE) %>% as_tibble()
# Sanai batch 3 read counts (80 samples)
sanai_3 <- read.table(file = paste0(WD,'/sinai/batch_3/summary/featureCounts.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE) %>% as_tibble()

# Sanai batch 1 meta data
sanai_meta_1 <-  read.table(file = paste0(WD,'/sinai/batch_1/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)
# Sanai batch 2 meta data
sanai_meta_2 <-  read.table(file = paste0(WD,'/sinai/batch_2/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE) %>% as_tibble()
# Sanai batch 3 meta data
sanai_meta_3 <-  read.table(file = paste0(WD,'/sinai/batch_3/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE) %>% as_tibble()

# Stanford RNASeq
###################################

# Sanford batch 1 read counts (320 samples)
sford_1 <- read.table(file = paste0(WD,'/stanford/batch_1/summary/rsem_genes_count.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE) %>% as_tibble()
# Sanford batch 2 read counts (320 samples)
sford_2 <- read.table(file = paste0(WD,'/stanford/batch_2/summary/featureCounts.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE) %>% as_tibble()

# Sanford batch 1 meta data
sford_meta_1 <-  read.table(file = paste0(WD,'/stanford/batch_1/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE) %>% as_tibble()
# Sanford batch 2 meta data
sford_meta_2 <-  read.table(file = paste0(WD,'/stanford/batch_2/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE) %>% as_tibble()

# General Phenotype Data
##################################
pheno_file <- paste0(WD,'/../phenotype/merged/merged_dmaqc_data2019-10-13.txt')
gen_pheno <- read.table(file = pheno_file, header = TRUE, sep = '\t', check.names = FALSE) %>% 
        as_tibble()

# Adjust column objects
########################
sford_meta_1$vial_label <- as.factor(sford_meta_1$vial_label)

# Adjust column objects
########################
# To factors
factor_cols <- c("labelid",
                 "viallabel",
                 "animal.registration.sex")
for(fc in factor_cols){
        gen_pheno[[fc]] <- as.factor(gen_pheno[[fc]])
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
        gen_pheno[[dc]] <- dmy(gen_pheno[[dc]])
}
# To Dates: Jan 2018 
gen_pheno <- gen_pheno %>%
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
        gen_pheno[[tc]] <- gen_pheno[[tc]] %>% as.character() %>% parse_time()
}

# To Times: 01:21
gen_pheno <- gen_pheno %>%
        mutate(acute.test.howlongshock = 
                       ifelse((gen_pheno$acute.test.howlongshock != ''), paste0('00:',acute.test.howlongshock), '00:00:00'))
gen_pheno$acute.test.howlongshock <- parse_time(gen_pheno$acute.test.howlongshock)

# Adjust labels within columns (e.g. 1 and 2 to female and male)
# see ../phenotype/merged/20191013_merged-column-dictionary_steep-edits.xlsx
########################
gen_pheno <- gen_pheno %>%
        mutate(animal.registration.sex = 
                       case_when(animal.registration.sex == '1' ~ 'Female',
                                 animal.registration.sex == '2' ~ 'Male'))
gen_pheno$animal.registration.sex <- as.factor(gen_pheno$animal.registration.sex)

#' # Incorporate Annotation from the phenotypic data

#+ Incorporate Pheno
################################################################################
######### Incorporate Phenotypic Annotations ###################################
################################################################################

# Replace the name of the vial label column (for sake of join)
names(gen_pheno) <- str_replace(names(gen_pheno), 'viallabel', 'vial_label')

# Perform a left join to incorporate phenotypic information
status <- left_join(sford_meta_1, gen_pheno, by = 'vial_label') %>% as.data.frame()

# Create a counts matrix
sford_1 <- sford_1 %>% as.data.frame()
row.names(sford_1) <- sford_1$gene_id
all.data <- sford_1 %>% select(-gene_id)
# Convert counts to integers
all.data.m <- as.matrix(all.data)
mode(all.data.m) <- "integer"

# "It is absolutely critical that the columns of the count matrix and the rows of the 
# column data (information about samples) are in the same order. DESeq2 will not make 
# guesses as to which column of the count matrix belongs to which row of the column data, 
# these must be provided to DESeq2 already in consistent order." ~ Mike Love
rownames(status) <- status$vial_label
all(rownames(status) %in% colnames(all.data.m))
all(rownames(status) == colnames(all.data.m))
all.data.m <- all.data.m[, rownames(status)]
all(rownames(status) == colnames(all.data.m))

##########################################################################
############ Annotate Data for BioJupies #################################
##########################################################################
# Annotate Data for BioJupies
# Annotate normalized counts
#biojupes <- all.data
#biojupes$ensembl <- rownames(biojupes)
#biojupes$symbol <- mapIds(org.Gg.eg.db, biojupes$ensembl, "SYMBOL", "ENSEMBL")
#biojupes <- biojupes[!duplicated(biojupes$symbol),]
#biojupes <- biojupes[!is.na(biojupes$symbol),]
#rownames(biojupes) <- biojupes$symbol
#biojupes <- biojupes %>% dplyr::select(symbol, everything())
#biojupes <- biojupes %>% dplyr::select(-ensembl)
#write.table(biojupes, file = paste0('./data/',date,'_biojupiesmatrix_steep.txt'), quote = FALSE, row.names = FALSE, sep = '\t')
##########################################################################

#' ## PCA Analysis
#' TODO: Generate a summary of the inferences from PCA.
#' 

#+ PCA Analysis (Unsupervised model)

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
dds
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
dds

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

#' ### Primary variable of interest, unsupervised
#' Let's examine our primary dichotomous relationship of interest--exercise vs control. How well does this relationship explain the variance in our data? Are there other factors influencing our data?
names(status)

# PCA Analysis (tissue)
# Examine the Distribution of male and female samples in Liver
rld.sub <- rld[ , (rld$Tissue %in% c("Liver") & 
                           rld$animal.registration.sex %in% c("Female")) ]

# Examine groups by exercise
DESeq2::plotPCA(rld.sub, intgroup ="animal.key.anirandgroup") +
        guides(color=guide_legend(title="Group"))
