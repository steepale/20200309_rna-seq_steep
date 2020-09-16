#'---
#' title: "PASS1A (Rat) RNA-Seq: Build the Gene Count Matrix (TPM)"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200309"
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
#' * TODO: Investigate samples for a correlation of 1 (duplicate samples) and remove them
#' * TODO: Consider removing control probes from the analysis
#' * TODO: Vet "suspect samples" in file: GET_release1_qc_report.pdf and decide if these samples should be removed or not.
#' * Examine if data demonstrate technical batch effects:
#'     * Sequencing batch (sequence date and location)
#' * Examine if data demonstrate biological batch effects:
#'     * Influence from time of day and circadian rhythm
#'     * Influence from feeding schedules
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
#BiocManager::install("org.Rn.eg.db")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2",
                "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr",
                "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt","feather",
                "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","reshape2")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

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
#' ##### Data files to load:
#' * RNA-Seq from Mt. Sanai
#'     * 3 sequencing batches & metadata
#' * RNA-Seq from Stanford
#'     * 2 sequencing batches & metadata

# Stops evaluation of code
knitr::opts_chunk$set(eval = FALSE)

#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# Mt Sanai RNASeq
###################################

# Sanai batch 1 read counts (320 samples)
sanai_1 <- read.table(file = paste0(WD,'/sinai/batch_1/summary/rsem_genes_tpm.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)
# Sanai batch 2 read counts (320 samples)
sanai_2 <- read.table(file = paste0(WD,'/sinai/batch_2/summary/rsem_genes_tpm.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)
# Sanai batch 3 read counts (80 samples)
sanai_3 <- read.table(file = paste0(WD,'/sinai/batch_3/summary/rsem_genes_tpm.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)

# Sanai batch 1 meta data
sanai_meta_1 <-  read.table(file = paste0(WD,'/sinai/batch_1/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)
# Sanai batch 2 meta data
sanai_meta_2 <-  read.table(file = paste0(WD,'/sinai/batch_2/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)
# Sanai batch 3 meta data
sanai_meta_3 <-  read.table(file = paste0(WD,'/sinai/batch_3/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)

# Stanford RNASeq
###################################

# Stanford batch 1 read counts (320 samples)
sford_1 <- read.table(file = paste0(WD,'/stanford/batch_1/summary/rsem_genes_tpm.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)
# Stanford batch 2 read counts (320 samples)
sford_2 <- read.table(file = paste0(WD,'/stanford/batch_2/summary/rsem_genes_tpm.txt'), 
                      header = TRUE, sep = '\t', check.names = FALSE)

# Stanford batch 1 meta data
sford_meta_1 <-  read.table(file = paste0(WD,'/stanford/batch_1/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)
# Stanford batch 2 meta data
sford_meta_2 <-  read.table(file = paste0(WD,'/stanford/batch_2/metadata/sample_metadata_20191010.csv'), 
                            header = TRUE, sep = ',', check.names = FALSE)
# Determine sequencing location and batch (batch by date)
#str(sford_meta_1)
#str(sford_meta_2)
#str(sanai_meta_3)

# Location under GET_site variable
#summary(sford_meta_2$GET_site)
#summary(sford_meta_2$GET_site)
#summary(sanai_meta_1$GET_site)
#summary(sanai_meta_2$GET_site)
#summary(sanai_meta_3$GET_site)

# Date (batch) under Seq_date variable
#unique(sford_meta_1$Seq_date)
#unique(sford_meta_2$Seq_date)
#unique(sanai_meta_1$Seq_date)
#unique(sanai_meta_2$Seq_date)
#unique(sanai_meta_3$Seq_date)

# Are the combination of these values unique?
#sford_meta_1 %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()
#sford_meta_2 %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()
#sanai_meta_1 %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()
#sanai_meta_2 %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()
#sanai_meta_3 %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()

# Combine the Stanford and Mt. Sanai Data
##########################################
# Count Matrixes

# Add custom annotation to column names
colnames(sanai_1) <- paste0(colnames(sanai_1),'_SN1')
colnames(sanai_2) <- paste0(colnames(sanai_2),'_SN2')
colnames(sanai_3) <- paste0(colnames(sanai_3),'_SN3')
colnames(sford_1) <- paste0(colnames(sford_1),'_SF1')
colnames(sford_2) <- paste0(colnames(sford_2),'_SF2')
sanai_meta_1$sample_key <- paste0(sanai_meta_1$vial_label, '_SN1')
sanai_meta_2$sample_key <- paste0(sanai_meta_2$vial_label, '_SN2')
sanai_meta_3$sample_key <- paste0(sanai_meta_3$vial_label, '_SN3')
sford_meta_1$sample_key <- paste0(sford_meta_1$vial_label, '_SF1')
sford_meta_2$sample_key <- paste0(sford_meta_2$vial_label, '_SF2')

# Combine counts
all.data.m <- cbind(sanai_1[,-1], sanai_2[,-1], sanai_3[,-1], sford_1[,-1], sford_2[,-1] ) %>% as.matrix()
rownames(all.data.m) <- sanai_1[,1]
# Convert counts to integers
#mode(all.data.m) <- "integer"

# Combine metadata (convert factors to characters to prevent error)
sanai_meta_1 <- data.frame(lapply(sanai_meta_1, as.character), stringsAsFactors=FALSE)
sanai_meta_2 <- data.frame(lapply(sanai_meta_2, as.character), stringsAsFactors=FALSE)
sanai_meta_3 <- data.frame(lapply(sanai_meta_3, as.character), stringsAsFactors=FALSE)
sford_meta_1 <- data.frame(lapply(sford_meta_1, as.character), stringsAsFactors=FALSE)
sford_meta_2 <- data.frame(lapply(sford_meta_2, as.character), stringsAsFactors=FALSE)
meta <- bind_rows(sanai_meta_1, sanai_meta_2, sanai_meta_3, sford_meta_1, sford_meta_2)

# Double check unique values
#dim(meta)
#meta %>% 
#        select(vial_label, GET_site, Seq_date) %>%
#        unique() %>% dim()

# General Phenotype Data
##################################
pheno_file <- paste0(WD,'/../phenotype/merged/merged_dmaqc_data2019-10-13.txt')
gen_pheno <- read.table(file = pheno_file, header = TRUE, sep = '\t', check.names = FALSE) %>% 
        as_tibble()

# Adjust column objects
########################
meta$vial_label <- as.factor(meta$vial_label)

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
        gen_pheno[[dc]] <- lubridate::dmy(gen_pheno[[dc]])
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

################################################################################
###################### Cleanup Memory ##########################################
################################################################################
rm("pacs...man","sanai_1","sanai_2","sanai_3","sanai_meta_1","sanai_meta_2","sanai_meta_3","sford_1","sford_2","sford_meta_1","sford_meta_2","pheno_file","factor_cols","date_cols","time_cols")
################################################################################

#' ## Incorporate Annotation from Phenotypic Data

#+ Incorporate Pheno
################################################################################
######### Incorporate Phenotypic Annotations ###################################
################################################################################

# Replace the name of the vial label column (for sake of join)
names(gen_pheno) <- str_replace(names(gen_pheno), 'viallabel', 'vial_label')

# To investigate if the "vial_label" key is indeed unique
#meta$vial_label %>% length()
#meta$vial_label %>% unique %>% length()
#gen_pheno$vial_label %>% length()
#gen_pheno$vial_label %>% unique %>% length()

# Perform a left join to incorporate phenotypic information
status <- left_join(meta, gen_pheno, by = 'vial_label') %>% as.data.frame()
#dim(status)

#' ## Generate New Annotations
#' * Sequencing batch
#' * TODO: Time since last fed
#' * Exercise and control (redundent -- for visualizations)
#' * Exercise group <= 4 hrs

#+ New Annotations
################################################################################
####### Build New Annotations and Adjust Objects in Columns ####################
################################################################################

# Sequencing Batch
status <- status %>% 
        mutate(Seq_batch = 
                       case_when((GET_site == 'Stanford' & Seq_date == '190426') ~ "Stanford_1",
                                 (GET_site == 'Stanford' & Seq_date == '190703') ~ "Stanford_2",
                                 (GET_site == 'MSSM' & Seq_date == '190409') ~ "MSSM_1",
                                 (GET_site == 'MSSM' & Seq_date == '190626') ~ "MSSM_2",
                                 (GET_site == 'MSSM' & Seq_date == '190723') ~ "MSSM_3")
        )
#table(status$Seq_batch)

# Exercise and Control
status <- status %>% 
        mutate(animal.key.exvsctrl = 
                       factor(case_when( grepl('Control',animal.key.anirandgroup) ~ "Control",
                                         grepl('Exercise',animal.key.anirandgroup) ~ "Exercise" ))
        )

# Exercise group <= 4 hrs
status <- status %>% 
        mutate(animal.key.exlt4 = 
                       factor(ifelse(
                               animal.key.anirandgroup %in% 
                                       c('Exercise - IPE',
                                         'Exercise - 0.5 hr',
                                         'Exercise - 1 hr',
                                         'Exercise - 4 hr'), 
                               '1','0'))
        )

# Create bins for exercise and control groups
# Thse bins represent the major clusters we see from RNASeq in Liver (likely other tissues)
status <- status %>% 
        mutate(animal.key.anirandgroup.bins.1 = 
                       factor(case_when(
                               animal.key.anirandgroup %in% 
                                       c('Exercise - IPE',
                                         'Exercise - 0.5 hr',
                                         'Exercise - 1 hr',
                                         'Exercise - 4 hr') ~ "Exercise - 0-4 hr",
                               animal.key.anirandgroup %in% 
                                       c('Exercise - 7 hr',
                                         'Exercise - 24 hr',
                                         'Exercise - 48 hr') ~ "Exercise - 7-48 hr",
                               animal.key.anirandgroup %in% 
                                       c('Control - IPE',
                                         'Control - 7 hr') ~ "Control - 0-7 hr")))

# Bin time of death of PASS1A rats and visualize time of death by major experimental groups
# Adjust the time of death column to be on an hourly scale
status$specimen.collection.t_death_hour = 
        hour(status$specimen.collection.t_death) + 
        minute(status$specimen.collection.t_death)/60 + 
        second(status$specimen.collection.t_death)/3600
#Create Bins
bins=c(paste0(rep(c(paste0(0,0:9),10:23), each=4),".", c("00",25,50,75))[-1],"24:00")
# Divide Time of Death into Subbins (ever quarter hour)
status$specimen.collection.t_death_bins = cut(status$specimen.collection.t_death_hour, breaks=seq(0, 24, 0.25), labels=bins)
#Reformat to Numeric
status$specimen.collection.t_death_bins <- as.numeric(as.character(status$specimen.collection.t_death_bins))
# Place the bins in even more bins (bins were choosen in subjectively/arbitrarily TODO: Articulate decision making)
status <- status %>% 
        mutate(specimen.collection.t_death_bins.type = 
                       factor(case_when(
                               (specimen.collection.t_death_bins >= 10 & 
                                        specimen.collection.t_death_bins <= 11.25) ~ "10.00-11.25",
                               (specimen.collection.t_death_bins >= 11.75 & 
                                        specimen.collection.t_death_bins <= 12.5) ~ "11.75-12.50",
                               (specimen.collection.t_death_bins >= 13.25 & 
                                        specimen.collection.t_death_bins <= 14.00) ~ "13.25-14.00",
                               (specimen.collection.t_death_bins >=  14.25 & 
                                        specimen.collection.t_death_bins <= 15.00) ~ "14.25-15.00",
                               (specimen.collection.t_death_bins >=17.00 & 
                                        specimen.collection.t_death_bins <= 18.00) ~ "17.00-18.00")))

# To factors
factor_cols <- c("sample_key",
                 "Seq_batch")
for(fc in factor_cols){
        status[[fc]] <- as.factor(status[[fc]])
}

# Releveling factors
status$animal.key.anirandgroup <- as.character(status$animal.key.anirandgroup)
status$animal.key.anirandgroup <- factor(status$animal.key.anirandgroup, 
                                         levels = c("Control - IPE",
                                                    "Control - 7 hr",
                                                    "Exercise - IPE",
                                                    "Exercise - 0.5 hr",
                                                    "Exercise - 1 hr",
                                                    "Exercise - 4 hr",
                                                    "Exercise - 7 hr",
                                                    "Exercise - 24 hr",
                                                    "Exercise - 48 hr"))

#' #### Data Save: Count Matrix & Metadata
#' Count matrix: motrpac/20200309_rna-seq_steep/data/20200309_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv
#' Metadata: motrpac/20200309_rna-seq_steep/data/20200309_rnaseq-meta-pass1a-stanford-sinai_steep.txt

#+ Data Save: Count Matrix & Metadata
################################################################################
############## Data Save: Count Matrix & Metadata ##############################
################################################################################

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
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(status) == colnames(all.data.m))

#' #### TODO: Data can easily be formatted for BioJupies if desired
#' #+ BioJupies
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

# Count matrix
out_file <- paste0(WD,'/data/20200309_rnaseq-countmatrix-pass1a-stanford-sinai-tpm_steep.csv')
#write.table(all.data.m, file=out_file, row.names = TRUE, quote = FALSE, sep = ',')

# Meatdata table
out_file <- paste0(WD,'/data/20200309_rnaseq-meta-pass1a-stanford-sinai-tpm_steep.txt')
#write.table(status, file=out_file, row.names = FALSE, quote = FALSE, sep = '\t')

#' #### Histograms and density plots demonstrate the frequency of sampling based on time of day as well as our binning strategies for incorporating time of day into analyses. 
#' #### Tissue type collectd by exercise/control group
#' TODO: Articulate decision making: bins were choosen subjectively/arbitrarily. Instead, we should come back once we have a nice cohort of circadian genes in certain tissues and cluster these dates into bins that better represent the (cosine) models of circadian rhythm. 

#+ Design Summary: Heatmaps, Histograms, Density Plots, and Tables
################################################################################
####### Design Summary: Heatmaps, Histograms, Density Plots, and Tables ########
################################################################################

# Generate a heatmap of major annotations of interest
# Annotations of interest:
# Samples: sample_key
# Sequencing Batch: Seq_batch
# Sex: animal.registration.sex
# Exercise/Control Group: animal.key.anirandgroup
# Time of Death (binned): specimen.collection.t_death_bins.type
# Tissue: Tissue

# Select annotations of interest
ann_df <- status %>%
        select(sample_key,
               Tissue,
               Seq_batch,
               animal.registration.sex,
               animal.key.anirandgroup,
               specimen.collection.t_death_bins.type)
annotations <- c("sample_key",
                 "Tissue",
                 "Seq_batch",
                 "animal.registration.sex",
                 "animal.key.anirandgroup",
                 "specimen.collection.t_death_bins.type")
for(anno in annotations){
        ann_df[[anno]] <- as.character(ann_df[[anno]])
}
# Adjust names
names(ann_df) <- c('Sample','Tissue','Sequencing Batch',
                   'Sex','Exercise/Control Group','Time of Death (binned)')
ann_df <- ann_df %>%
        select('Sample','Exercise/Control Group','Time of Death (binned)',
               'Sex','Sequencing Batch','Tissue')
# Melt dataframe
ann_df_melt <- melt(ann_df, id.var = 'Sample')
# Adjust colors
n <- length(unique(ann_df_melt$value))
colors <- colorRampPalette(c("blue", "green", "yellow", "red"))(n)

#' #### The distribution of Samples by Annotations of Primary Interest
# Plot heatmap
ggplot(ann_df_melt, 
       aes(variable, Sample)) + 
        geom_tile(aes(fill = value),
                  colour = "white") +
        scale_fill_manual(values=colors) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        theme(legend.position="none") +
        theme(axis.title.x=element_blank()) +
        theme(text = element_text(size=20))

# Tissue type collectd by exercise/control group
table(status$Tissue, status$animal.key.anirandgroup)
# Tissue type collectd by time of death
table(status$Tissue, status$specimen.collection.t_death_bins.type)
# Demonstrate the confounding
table(status$animal.key.anirandgroup, status$specimen.collection.t_death_bins.type)

#Histogram
ggplot(status, aes(x=specimen.collection.t_death_bins,
                   color=specimen.collection.t_death_bins.type, 
                   fill=specimen.collection.t_death_bins.type)) +
        geom_histogram(aes(y=..density..), position="identity", alpha=0.5, bins = 96) +
        labs(title="Histogram: \nAdjusted Bins for Time of Death",
             x="Time of Death (Hour scale)", 
             y = "Density") +
        scale_x_continuous(breaks=seq(0, 24, 0.5))
# Density plot of sampling based on experimental group
ggplot(status, aes(x=specimen.collection.t_death, 
                   fill=animal.key.anirandgroup)) +
        geom_density(alpha=0.4) +
        labs(title="Density Plot \nTime of Death for Experimental Groups",
             x="Time of Death (Hour scale)",
             y = "Density") +
        guides(fill=guide_legend(title="Exercise/Control Groups"))

#' ## PCA Visualization of Sequencing Batches (Unsupervised)
#' TODO: Generate a summary of major inferences
#' 

#+ PCA Sequencing Batches (Unsupervised)
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
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet
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
#counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
#counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
#counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
#counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")

#' #### Let's examine the stucture in the data to see if there is significant sample correlation
#' It looks like much of the sample correlation is driven by sequencing batch
s <- svd(norm_counts)
dim(s$v)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
image ( cor(norm_counts) ,col=cols,zlim=c(-1,1))

#' #### Variance-explained plot: explained variance for PCs
#' ##### This is what independent data would look like:
norm_counts0 <- matrix( rnorm( nrow(norm_counts)*ncol(norm_counts) ) , nrow(norm_counts), ncol(norm_counts) )
d0 <- svd(norm_counts0)$d
plot(d0^2/sum(d0^2),ylim=c(0,.25))

#' ##### This is what these data look like: Which shows that a majority of variance is explained by just 2 principle components
plot(s$d^2/sum(s$d^2))

#' #### Here We examine all samples by batch and tissue. In our experience, if there is no batch effect, then samples should cluster by tissue type. However, if samples cluster by "batch," which we define as sequencing date and location, then a batch effect might be occurring. 
#' 
#' #### When we examine groups by tissue, it's difficult to distinguish tissues colors, but it seems like mutliple tissues are represented in different clusters.
DESeq2::plotPCA(rld, intgroup ="Tissue") +
        guides(color=guide_legend(title="Tissues"))
#' 
#' #### When we visualize sequencing batches see that Stanford batch 1 is significantly isolated.
DESeq2::plotPCA(rld, intgroup ="Seq_batch") +
        guides(color=guide_legend(title="Batches"))

#' #### Let's examine the stucture in the data to see if there is significant sample correlation
#' It looks like much of the sample correlation is driven by sequencing batch
s <- svd(norm_counts)
dim(s$v)
cols=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100)
image ( cor(norm_counts) ,col=cols,zlim=c(-1,1))


#' #### Variance-explained plot: explained variance for PCs
#' ##### This is what independent data would look like:
norm_counts0 <- matrix( rnorm( nrow(norm_counts)*ncol(norm_counts) ) , nrow(norm_counts), ncol(norm_counts) )
d0 <- svd(norm_counts0)$d
plot(d0^2/sum(d0^2),ylim=c(0,.25))
# Cleanup
rm(d0,norm_counts0)

#' ##### This is what these data look like: Which shows that a majority of variance is explained by just 2 principle components
plot(s$d^2/sum(s$d^2))

#' #### Boxplot of PCs
#' ##### If we examine the first 4 PC's and startify the data by batch, we can see that the majority of variance in PCs 1 and two are associated with Stanford Batch 1.

# length(unique(col_data$Seq_batch))
variable <- col_data$Seq_batch
mypar(2,2)
for(i in 1:4){
        boxplot(split(s$v[,i],variable),las=2,range=0)
        stripchart(split(s$v[,i],variable),add=TRUE,vertical=TRUE,pch=1,cex=.5,col=1)
}

#' #### We can use alalysis of variance to see which PCs correlate with batch
#' ##### As expected, PCs 1 and 2 correlate strongly with sequencing batch. In total, the first 7 PCs or so show a relatively strong correlation to batch.

corr <- sapply(1:ncol(s$v),function(i){
        fit <- lm(s$v[,i]~as.factor(variable))
        return( summary(fit)$adj.r.squared  )
})
mypar()
plot(seq(along=corr), corr, xlab="PC")

#' #### F-statistics comparing within sequencing batch to across sequencing batch variability show how extreme the F-statistic is for PCs 1 and 2.
i <- 1
Fstats<- sapply(1:ncol(s$v),function(i){
        fit <- lm(s$v[,i]~as.factor(variable))
        Fstat <- summary(aov(fit))[[1]][1,4]
        return(Fstat)
})
mypar()
plot(seq(along=Fstats),sqrt(Fstats))
p <- length(unique(variable))
abline(h=sqrt(qf(0.995,p-1,ncol(s$v)-1)))

#' #### It's possible that certain tissues were isolated in certain batches. Let's examine what tissues were sequenced in respective batches.

# Create a contingency table of Tissues sequenced by batch
ctable <- status %>%
        select(Seq_batch,Tissue) %>% 
        table() %>% as.data.frame()

#' #### This balloon plot demonstrates that sample tissue types were not evenly distributed across batches of RNASeq data. The overlap between batches and tissues is far from ideal. Size of dots and color of dots demonstrate the same attribute: frequency.
#' 
# Create a balloon plot
theme_set(theme_pubr())
ggballoonplot(ctable, fill = "value") +
        scale_fill_viridis_c(option = "C") +
        ggtitle("Tissue Samples Distributed Across Sequencing Batches")

#' #### Tissue types sequenced across multiple batches include:
#' * Gastrocnemius
#' * Heart
#' * Kidney
#' * Lung

ctable %>% 
        filter(Tissue %in% c("Gastrocnemius","Heart","Kidney","Lung")) %>%
        ggballoonplot(fill = "value") + 
        scale_fill_viridis_c(option = "C") +
        ggtitle("Tissue Samples Distributed Across Sequencing Batches")

#' #### Some of these above samples are "reference samples" sequenced across batches
#' Note: Notice how the maximum frequency has dramatically decreased but the circle size has not.
ref_table <- status %>%
        filter(Sample_category == 'ref') %>%
        select(Seq_batch,Tissue) %>% 
        table() %>% as.data.frame()
ref_table %>% 
        ggballoonplot(fill = "value") + 
        scale_fill_viridis_c(option = "C") +
        ggtitle("Reference Sample Sequence Distribution Across Batches")

#' #### We examine a PCA plot of **only** samples with tissues represented across mutliple batches.
#' 
#' ##### The tissue types in which samples cluster by batch:
#' * Gastrocnemius
#' 
#' ##### The tissue types in which samples cluster by tissue:
#' * Heart
#' * Kidney
#' * Lung
rld.sub <- rld[ , (rld$Tissue %in% c("Gastrocnemius","Heart","Kidney","Lung")) ]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("Tissue","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Seq_batch, shape=Tissue)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle("Tissues Sequenced Across More than One Batch")

#' #### Gastrocnemius samples cluster by sequencing batch.
#' 
#' Note: variance is extreme.
#' 
#' Also shown: variance is not driven by sex.
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

#' #### Heart samples cluster by sex. 
#' 
#' Grey samples represent reference samples (some distributed across sequencing sites).
rld.sub <- rld[ , (rld$Tissue == "Heart") ]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Heart Samples Sequenced Across Batches")

#' #### Kidney samples cluster by sex, but were all primarily sequenced in one batch (except for a pair of reference samples)
rld.sub <- rld[ , (rld$Tissue == "Kidney") ]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle("Kidney Samples Sequenced Across Batches")

#' #### Lung samples cluster by sex. 
rld.sub <- rld[ , (rld$Tissue == "Lung") ]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle("Lung Samples Sequenced Across Batches")

#' #### Now let's examine the distribution of **only** reference samples across tissues and batches.
rld.sub <- rld[ , (rld$Sample_category == 'ref') ]
pcaData <- DESeq2::plotPCA(rld.sub, 
                           intgroup=c("animal.registration.sex","Seq_batch","Tissue"), 
                           returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Seq_batch)) +
        geom_point(size=3) +
        geom_text_repel(aes(label=Tissue)) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed() +
        ggtitle("Reference Samples Sequenced Across Batches")

#' #### Inference:
#' These plots demonstrate agreement across most tissues sequenced across multiple batches. However, Gastocnemius was the only tissue and reference sample tissue to be sequenced across Stanford batch 1 and other batches (Mt. Sanai batch 1 & Stanford batch 2). These samples help suggest that samples from Stanford batch 1 are experiencing a heavy batch effect. All Stanford batch 1 samples cluster together, while other batches are more spread out by tissue. One might argue that the tissues in Stanford batch 1 are similar to one another (e.g. that they are fatty): White Adipose, Brown Adipose, Liver & Gastrocnemius. However, Gastrocnemius is not a fatty tissue whatsoever and the observation that a large amount of variance drives this muscle to cluster with fatty tissues **only in Stanford batch 1** suggests a strong batch effect. Finally, most of the variance in these data are driven by Stanford batch 1 samples.

#' #### Systematic Search for Categorical Metadata to Correct for Batch
#' As we systematically searched the meta data and phenotypic annotation for possible culprits of batch. The categorical variables most likely associated with batch (i.e. Seq_batch) are:
#' * Seq_date (not suprising)
#' * Seq_flowcell_ID 
for(clm in colnames(status)[colnames(status) != 'Seq_batch'] ){
        if(is.character(status[[clm]]) | is.factor(status[[clm]])){
                col_num <- status %>% select("Seq_batch", clm) %>% table() %>% ncol()
                if(col_num < 10 && clm %in% c("Seq_date","Seq_flowcell_ID")){
                        print(clm)
                        status %>% select("Seq_batch", clm) %>% table() %>% print()
                        print("###########")
                        print("")
                }
                
        }
        
}

#' TODO: Compare differential gene expression analysis between Gastroc samples from different batches
#' * Show a correlation matrix of svd data and expect uneven correlatation -- structure in data
#'     * Look for certain samples being correlated more than other samples (uninform structure)
#' * Look at "d's" for randomized data of singular value composition -- should be straight line around 0
#' * Look at d's from svd -- see if it shows structure
#' * Show batch with PCA or MDS
#' * Use box-plots and split them by factor we think is related to batch
#'     * Could support visualization with Anova by showing that PC correlated with batch category
#' * Show correlation matrix and correlation table with factor analysis
#'     * Show a fitted model and see how much of the variance can be explained
#'     * Correlation matrix before and after SVA
#' * Histogram of p-values next to volcano plot (before SVA)
#' * Histogram of p-values next to volcano plot (after SVA)
#' * Heatmap of statified genes (before SVA)
#' * Heatmap of statified genes (after SVA)
#' 

#' There's no categorical variable (alone) that can isolate Stanford batch 1 from the rest of the group. Here we experiment with a new binary variable, "Seq_batch_effect", and run a supervised model to determine if we can correct for the suspected batch effect from Stanford batch 1.

#+ Supervised PCA: Correction for Stanford Batch 1
##########################################################################
############ Supervised PCA: Correction for Stanford Batch 1 #############
##########################################################################

#' #### The Seq_batch_effect varibale:
#status <- status %>%
#        mutate(Seq_batch_effect = as.factor(ifelse(Seq_batch == 'Stanford_1', '1', '0')))
#status %>%
#        select(Seq_batch, Seq_batch_effect) %>% table()

# Perform supervised clustering
#col_data <- status # To look at all samples 
#design = ~ pcaData$PC1 + 1 # Primary variable needs to be last
#title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
#dds <- DESeqDataSetFromMatrix(countData = count_data,
#                              colData = col_data,
#                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
#keep <- rowSums(counts(dds))/ncol(dds) >= 10
#dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet
#dds

#dds <- estimateSizeFactors(dds)
#rs <- rowSums(counts(dds))
# Normalize the counts
#rld <- vst(dds) #vst and rlog comparable with all samples


#' #### Now let's examine the distribution of **only** reference samples across tissues and batches.
#rld.sub <- rld[ , (rld$Sample_category == 'ref') ]
#pcaData <- DESeq2::plotPCA(rld.sub, 
#                           intgroup=c("animal.registration.sex","Seq_batch","Tissue"), 
#                           returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, color=Seq_batch)) +
#        geom_point(size=3) +
#        geom_text_repel(aes(label=Tissue)) +
#        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#        coord_fixed() +
#        ggtitle("Reference Samples Sequenced Across Batches")

#' #### The TODOs:
#' * Demonstrate a volcano plot showing how many genes 
#' * Examine if the batch from Stanford 1 can be corrected for possible culprits:
#'     * Sequencing machine
#'     * Sequencing library
#'     * Time from tissue collection to sequencing
#'     * Experimental condition

# #### Sesh:
session_info()

# Stops evaluation of code
#knitr::opts_chunk$set(eval = FALSE)



