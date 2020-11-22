
#' ## Goals of Analysis:
#' 
#' ## Setup the Environment

#+ Setup Environment, message=FALSE, results='hide', warning = FALSE
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/home/acsteep/motrpac/20200309_rna-seq_steep'
#WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("gam")
#install.packages("XML")
# for(p in packages){
#         if(!require(p, character.only = T)){
#                 install.packages(p)
#         }
#         library(p, character.only = T)
# }

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","reshape2","xtable","kohonen","caret","enrichR","gplots","tiff","splines","gam","DESeq2","car","KEGGREST")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
slice <- dplyr::slice
filter <- dplyr::filter
counts <- DESeq2::counts
map <- purrr::map
seq_range <- modelr::seq_range

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
source(paste0(WD,'/functions/circleFun.R'))

# A function to collect the non-parametric p-value
################################################################################
p.test.cos <- function(data,TOD,iter,every=10) {
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

# A function to collect the non-parametric p-value
################################################################################
p.test.ce <- function(data,TOD,iter,every=10) {
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


# Collect the percent variance explained 
################################################################################
PVE_e <- function(df){
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output pve exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained 
################################################################################
PVE_e2 <- function(df){
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output pve exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained 
################################################################################
PVE_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_res <- ss_res/ss_total
        # Output all 3 PVE (circ, exer, res)
        pve_circ
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb_e <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb2_e <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_exer
}
################################################################################

# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'ns(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_circ
}
################################################################################


# Collect the percent variance explained from a sin or cosine aspect of a model
################################################################################
PVE_comb2_c <- function(df){
        # Collect the sum of squares from each model term
        ss_sin <- df %>% 
                filter(term == 'SIN(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_cos <- df %>% 
                filter(term == 'COS(specimen.collection.t_death_hour_mc)') %>%
                select(sumsq) %>% as.numeric()
        ss_ns <- df %>% 
                filter(term == 'poly(specimen.collection.t_exercise_hour_sqrt, df = 4)') %>%
                select(sumsq) %>% as.numeric()
        ss_res <- df %>% 
                filter(term == 'Residuals') %>%
                select(sumsq) %>% as.numeric()
        ss_total <- df %>% select(sumsq) %>% sum()
        # Collect the PVEs
        pve_circ <- (ss_sin + ss_cos)/ss_total
        pve_exer <- ss_ns/ss_total
        pve_res <- ss_res/ss_total
        # Output exercise
        pve_circ
}
################################################################################

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

by_gene_by_tissue_df <- tibble()
models_df <- data.frame()
#TISSUE <- "Kidney"
# for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
for(TISSUE in c('White Adipose')){
        print(TISSUE)
        TISSUE1 <- TISSUE
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
        
        # TODO: Take this if statement from this script and incorporate it into the bollinger script
        # Files last saved in: 20200603_rnaseq-tissue-data-assembly_steep.R
        
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
                
                # Determine the median values by which to center
                median_c0 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Control - IPE') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e24 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 24 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e48 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 48 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e0.5 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 0.5 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e0 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - IPE') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e1 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 1 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e4 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 4 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_c7 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Control - 7 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                median_e7 <- col_data %>%
                        filter(Tissue == 'Kidney') %>%
                        filter(animal.key.anirandgroup == 'Exercise - 7 hr') %>%
                        select(specimen.collection.t_death_hour) %>%
                        unlist() %>% as.numeric() %>% median()
                
                # Apply the median centered groups to time of death (modeling purposes)
                col_data <- col_data %>%
                        mutate(specimen.collection.t_death_hour_mc = 
                                       case_when(animal.key.anirandgroup == 'Control - IPE' ~
                                                         median_c0,
                                                 animal.key.anirandgroup == 'Control - 7 hr' ~
                                                         median_c7,
                                                 animal.key.anirandgroup == 'Exercise - IPE' ~
                                                         median_e0,
                                                 animal.key.anirandgroup == 'Exercise - 0.5 hr' ~
                                                         median_e0.5,
                                                 animal.key.anirandgroup == 'Exercise - 1 hr' ~
                                                         median_e1,
                                                 animal.key.anirandgroup == 'Exercise - 4 hr' ~
                                                         median_e4,
                                                 animal.key.anirandgroup == 'Exercise - 7 hr' ~
                                                         median_e7,
                                                 animal.key.anirandgroup == 'Exercise - 24 hr' ~
                                                         median_e24,
                                                 animal.key.anirandgroup == 'Exercise - 48 hr' ~
                                                         median_e48))
                
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
        # gf_file <- paste0(WD,'/data/20200603_Rnor-6.0.96-GRanges_steep.sqlite')
        # Rn_TxDb <- loadDb(gf_file)
        # # Define Female specific sex genes (X chromosome)
        # # To examine chromosome names
        # seqlevels(Rn_TxDb)[1:23]
        # # Extract genes as GRanges object, then names
        # X_genes_gr <- genes(Rn_TxDb, columns = 'TXCHROM', filter = list(tx_chrom=c('X')))
        # # Collect ensembl gene ids for female specific genes
        # X_ens_id <- names(X_genes_gr)
        # # Examine the gene symbols
        # X_sym <- mapIds(org.Rn.eg.db, names(X_genes_gr), 'SYMBOL', 'ENSEMBL')
        # # Extract genes as GRanges object, then names
        # Y_genes_gr <- genes(Rn_TxDb, columns = 'TXCHROM', filter = list(tx_chrom=c('Y')))
        # # Collect ensembl gene ids for female specific genes
        # Y_ens_id <- names(Y_genes_gr)
        # sex_ens_id <- c(X_ens_id,Y_ens_id)
        # # Examine the gene symbols
        # Y_sym <- mapIds(org.Rn.eg.db, names(Y_genes_gr), 'SYMBOL', 'ENSEMBL')
        
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
        # tod_cols <- tod_cols %>%
        #         mutate(calculated.variables.deathtime_after_acute =
        #                        ifelse(animal.key.anirandgroup == 'Control - IPE', 
        #                               calculated.variables.deathtime_after_acute - 3600,
        #                               calculated.variables.deathtime_after_acute))
        # tod_cols <- tod_cols %>%
        #         mutate(specimen.collection.t_exercise_hour_sqrt = ifelse(
        #                 calculated.variables.deathtime_after_acute < 0, 
        #                 (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)*(-1), 
        #                 (sqrt(abs(calculated.variables.deathtime_after_acute))/60/60)))
        
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
                filter(sample_key %!in% OUTLIERS) %>%
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
        # zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% 
        #         nrow() %>% as.character()
        # print(paste0('genes removed with < 1 count on average across samples: ',zero_n))
        # reads_n <- 1
        # keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
        # dds2 <- dds1[keep,]
        #' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
        #' TODO: Critic from Jun: Here we are removing features that have a low average expression. This may be removing important features that might have zero counts in some samples and higher counts in specific groups. Consider developing an algorithm that will account for features with expression in n or more samples.
        #' dds2
        #' filter_n <- nrow(dds1) - nrow(dds2) - as.numeric(zero_n)
        #' filter_p <- filter_n/(nrow(dds1) - as.numeric(zero_n))
        #' total_n <- nrow(dds1) - nrow(dds2)
        #' #' ##### Note: Number of genes with average counts between zero and 1 is `r zero_n` but removing reads with less than or equal to `r reads_n` removes an additional `r filter_n` features or removes `r filter_p*100`% of the non-zero reads (total of `r total_n` features removed).
        # dds <- dds2
        dds <- dds1
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
        
        #' ### Adjust Variance
        
        #+ Adjust Variance
        ################################################################################
        ########### Adjust Variance  #######################################
        ################################################################################
        #adj_var <- 'animal.key.batch'
        #adj_var <- 'animal.registration.sex'
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
        p <- ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
                geom_point(size=3) +
                #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
                xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
                #coord_fixed() +
                ggtitle(paste0("PCA of ",TISSUE)) +
                guides(color=guide_legend(title="animal.key.anirandgroup")) +
                scale_color_manual(values=ec_colors) +
                theme(legend.title=element_blank())
        plot(p)
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
        t_counts <- setNames(reshape2::melt(tod_counts), 
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
        poly4_mod <- function(df) {
                lm(count ~ poly(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        
        # Generate a model function for the dataframes
        sin_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc),
                   data = df)
        }
        
        # Generate a model that combines circadian with exercise (circadian first)
        ce_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc) +
                           ns(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        ce2_mod <- function(df) {
                lm(count ~ SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc) +
                           poly(specimen.collection.t_exercise_hour_sqrt, df = 4), data = df)
        }
        # Generate a model that combines circadian with exercise (exercise first)
        ec_mod <- function(df) {
                lm(count ~ ns(specimen.collection.t_exercise_hour_sqrt, df = 4) +
                           SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc), data = df)
        }
        ec2_mod <- function(df) {
                lm(count ~ poly(specimen.collection.t_exercise_hour_sqrt, df = 4) +
                           SIN(specimen.collection.t_death_hour_mc) + 
                           COS(specimen.collection.t_death_hour_mc), data = df)
        }
        
        # Add the gene symbol
        by_gene_df$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")
        by_gene_df7$SYMBOL_RAT = mapIds(org.Rn.eg.db, as.character(by_gene_df7$ENSEMBL_RAT), "SYMBOL", "ENSEMBL")
        
        # In case you'd like to subset the data
        #by_gene_df_bk <- by_gene_df
        #by_gene_df <- by_gene_df_bk
        
        # # Load in the MANOVA files (TODO add to top)
        # manova_file <- paste0(WD,'/data/20200413_rnaseq-tissue-manova-table_steep.txt')
        # manova_df <- read.table(file = manova_file ,sep = '\t', header = T, check.names = F) %>%
        #         as_tibble()
        # 
        # 
        # keepers <- manova_df %>% 
        #         filter(TISSUE == 'Kidney') %>%
        #         filter(MANOVA_PVAL <= 0.05) %>%
        #         select(ENSEMBL_RAT) %>% unlist() %>% as.character()
        by_gene_df <- by_gene_df %>%
                ungroup()
                # filter(ENSEMBL_RAT %in% keepers) #%>%
                # dplyr::sample_n(size = 50)
        # Arntl, Tbx10
        # by_gene_df
        # bin_df %>%
        #         filter(TISSUE == 'Kidney') %>%
        #         arrange(desc(pve_ce2_e))
        
        # by_gene_df <- by_gene_df %>% filter(SYMBOL_RAT == 'Nr4a2')
        # by_gene_df <- by_gene_df %>% filter(SYMBOL_RAT == 'Hsp90aa1')
        # by_gene_df <- by_gene_df %>% filter(SYMBOL_RAT == 'Arntl')
        # by_gene_df <- by_gene_df %>% filter(SYMBOL_RAT == 'Pdk4')
        # by_gene_df <- by_gene_df %>% filter(SYMBOL_RAT %in% kps)
        
        # kps <- d %>%
        #          filter(GNAME != '') %>%
        #          select(SYMBOL_RAT, BIN_TYPE, GROUP_ce_n,pve_ce2_e,pve_ce2_c, pve_e2, pve_c) %>%
        #          arrange(desc(pve_ce2_e)) %>% select(SYMBOL_RAT) %>% head(n = 20) %>%
        #         unlist() %>% as.character()
        
        # Capture the rejects mdl <- 'ce2'
        rejects_df <- data.frame()
        for(mdl in c('sin','poly4','ce2')){
                # circadian with exercise Model (circadian first)
                ################################################################################
                # Collect variables
                ( model_col <- as.symbol(paste0(mdl,'_model')) )
                anova_int <- as.symbol(paste0(mdl,'_anova_int'))
                anova_col <- as.symbol(paste0(mdl,'_anova'))
                resid_col <- as.symbol(paste0(mdl,'_resid'))
                metrics_col <- as.symbol(paste0(mdl,'_metrics'))
                summary_col <- as.symbol(paste0(mdl,'_summary'))
                deviance_col <- as.symbol(paste0(mdl,'_deviance'))
                MODEL <- match.fun(paste0(mdl,'_mod'))
                # Run models and save as a column
                by_gene_df <- by_gene_df %>%
                        mutate(!!model_col := map(data, MODEL))
                # Examine the ANOVA report on models
                #car::Anova defaults to type 2
                # Filter out genes with deviance of zero
                
                by_gene_df <- by_gene_df %>%
                        mutate(!!deviance_col := map_dbl(!!model_col, deviance))
                assign(paste0(mdl,'_rejects_df'), (by_gene_df %>%
                                                           filter(!!deviance_col <= sqrt(.Machine$double.eps))))
                
                rejects_df <- get(paste0(mdl,'_rejects_df'))
        }
        
        # Perform a series of analyses for each model
        for(mdl in c('sin','poly4','ce2')){
                # circadian with exercise Model (circadian first)
                ################################################################################
                # Collect variables
                ( model_col <- as.symbol(paste0(mdl,'_model')) )
                anova_int <- as.symbol(paste0(mdl,'_anova_int'))
                anova_col <- as.symbol(paste0(mdl,'_anova'))
                resid_col <- as.symbol(paste0(mdl,'_resid'))
                metrics_col <- as.symbol(paste0(mdl,'_metrics'))
                summary_col <- as.symbol(paste0(mdl,'_summary'))
                deviance_col <- as.symbol(paste0(mdl,'_deviance'))
                MODEL <- match.fun(paste0(mdl,'_mod'))
                # Run models and save as a column
                by_gene_df <- by_gene_df %>%
                        mutate(!!model_col := map(data, MODEL))
                # Examine the ANOVA report on models
                #car::Anova defaults to type 2
                # Filter out genes with deviance of zero
                
                by_gene_df <- by_gene_df %>%
                        mutate(!!deviance_col := map_dbl(!!model_col, deviance))
                by_gene_df <- by_gene_df %>%
                        filter(!!deviance_col > sqrt(.Machine$double.eps))
                by_gene_df <- by_gene_df %>%
                        mutate(!!anova_int := map(!!model_col, car::Anova)) %>%
                        mutate(!!anova_col := map(!!anova_int, broom::tidy)) %>%
                        select(-all_of(anova_int))
                # Add the residuals
                # by_gene_df <- by_gene_df %>%
                #   mutate(!!resid_col := map2(data, !!model_col, modelr::add_residuals))
                # # Examine the model metrics
                by_gene_df <- by_gene_df %>%
                   mutate(!!metrics_col := map(!!model_col, broom::glance))
                # # Examine some model summaries
                by_gene_df <- by_gene_df %>%
                   mutate(!!summary_col := map(!!model_col, summary))
        }
        
        # Collect all the PVEs
        pve_df <- by_gene_df %>%
                # select(ENSEMBL_RAT, SYMBOL_RAT,gam_anova,sin_anova,ce_anova,ec_anova,poly4_anova,ec2_anova,ce2_anova)
                select(ENSEMBL_RAT, SYMBOL_RAT,sin_anova,poly4_anova,ce2_anova)
        # Collect the PVE from models and different aspects of combined models
        names(by_gene_df$data[[1]])
        by_gene_df$data[[1]][['specimen.collection.t_death_hour']] %>%
                as.data.frame() %>% print(row.names = F)
        pve_df <- pve_df %>%
                # mutate(pve_e = map_dbl(gam_anova, PVE_e)) %>%
                mutate(pve_e2 = map_dbl(poly4_anova, PVE_e2)) %>%
                mutate(pve_c = map_dbl(sin_anova, PVE_c)) %>%
                # mutate(pve_ec_c = map_dbl(ec_anova, PVE_comb_c)) %>%
                # mutate(pve_ec_e = map_dbl(ec_anova, PVE_comb_e)) %>%
                # mutate(pve_ce_c = map_dbl(ce_anova, PVE_comb_c)) %>%
                # mutate(pve_ce_e = map_dbl(ce_anova, PVE_comb_e)) %>%
                # mutate(pve_ec2_c = map_dbl(ec2_anova, PVE_comb2_c)) %>%
                # mutate(pve_ec2_e = map_dbl(ec2_anova, PVE_comb2_e)) %>%
                mutate(pve_ce2_c = map_dbl(ce2_anova, PVE_comb2_c)) %>%
                mutate(pve_ce2_e = map_dbl(ce2_anova, PVE_comb2_e)) %>%
                # select(-gam_anova,-sin_anova,-ce_anova,-ec_anova,
                #        -poly4_anova,-ce2_anova,-ec2_anova) %>%
                select(-sin_anova,-poly4_anova,-ce2_anova) %>%
                mutate(TISSUE = TISSUE1)
        pve_rejects_df <- data.frame(ENSEMBL_RAT = rejects_df$ENSEMBL_RAT,
                                     SYMBOL_RAT = rejects_df$SYMBOL_RAT,
                                     pve_e2 = 0,
                                     pve_c = 0,
                                     pve_ce2_c = 0,
                                     pve_ce2_e = 0,
                                     TISSUE = TISSUE1)
        pve_df <- rbind(pve_df, pve_rejects_df)
        
        # Concatenate the dataframes
        #models_df <- pve_df
        models_df <- rbind(models_df, pve_df)
        # Save dfs
        by_gene_df$TISSUE <- TISSUE1
        # by_gene_by_tissue_df <- rbind(by_gene_by_tissue_df, by_gene_df)
        by_gene_by_tissue_df_file <- paste0(WD,'/data/20200727_rnaseq-',TISSUE1,'-models-data_steep.rds')
        # saveRDS(by_gene_df, file = by_gene_by_tissue_df_file)
        # Save the final output table
        models_file <- paste0(WD,'/data/20200603_rnaseq-',TISSUE1,'-models-pve-table2_steep.txt')
        write.table(models_df, file = models_file,sep = '\t',row.names = F,quote = F)
}
