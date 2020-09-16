#'---
#' title: "PASS1A Rat Tissue: -- Data Assembly"
#' author: "Alec Steep"
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

#' ## Setup the Environment

#+ Setup Environment, message=FALSE, results='hide', warning = FALSE
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source('https://bioconductor.org/biocLite.R')
BiocManager::install('bigstatsr')
#install.packages('tidyverse')

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr')
lapply(pacs...man, FUN = function(X) {
        do.call('library', list(X)) })

#' ## Load Custom Functions

#+ Functions, message=FALSE, results='hide', warning = FALSE
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

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
#' #### Section to generate a table of major decisions

# TISSUE: Hypothalamus, Liver, Kidney, Aorta, Adrenal, Brown Adipose, Cortex, Gastrocnemius, Heart, Hippocampus,Lung,Ovaries,PaxGene,Spleen,Testes, White Adipose

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

# Create the
df_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
names(df_tbl) <- c('Tissue','Tis','Seq_Batch','Misidentified',
                   'Adjusted_Variables','Outliers','Formula')
for(TISSUE in c('Kidney')){
        # for(TISSUE in c('Hypothalamus', 'Kidney', 'Aorta', 'Adrenal', 'Brown Adipose', 'Cortex', 'Gastrocnemius', 'Heart', 'Hippocampus','Lung','Ovaries','Spleen', 'White Adipose','Liver','Testes')){
        # Declare Outliers
        if(TISSUE == 'Kidney'){
                TIS <- 'KID'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Liver'){
                TIS <- 'LIV'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                OUTLIERS <- c('None')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
        }else if(TISSUE == 'Hypothalamus'){
                TIS <- 'SCN'
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('90010015402_SF2','90011015402_SF2')
                MIS_ID <- c('90010015402_SF2','90011015402_SF2')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Aorta'){
                TIS <- 'AOR'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90159016502_SN2')
        }else if(TISSUE == 'Gastrocnemius'){
                TIS <- 'SKM'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Heart'){
                TIS <- 'HAT'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('90052015802_SN1')
        }else if(TISSUE == 'Adrenal'){
                TIS <- 'ADG'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Brown Adipose'){
                TIS <- 'BAT'
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('90121016905_SF1','90014016905_SF1','90128016905_SF1','90047016905_SF1','90001016905_SF1',
                             '90018016905_SF1')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('90023016905_SF1','90047016905_SF1','90128016905_SF1')
        }else if(TISSUE == 'White Adipose'){
                TIS <- 'WAT'
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Cortex'){
                TIS <- c('COR')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Hippocampus'){
                TIS <- c('HIP')
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('90115015202_SF2','90041015202_SF2')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90040015202_SF2','90115015202_SF2','90041015202_SF2',
                              '90009015202_SF2','90005015202_SF2','90139015202_SF2')
        }else if(TISSUE == 'Lung'){
                TIS <- 'LUNG'
                OUTLIERS <- c('None')
                MIS_ID_VAR <- c('animal.key.batch')
                PRI_VAR <- c('animal.key.batch','animal.registration.sex')
                MIS_SUS <- c('90120016604_SN3', '90146016604_SN3', '90112016604_SN3',
                             '90129016604_SN3', '90124016604_SN3','90114016604_SN3',
                             '90112016604_SN3','90126016604_SN3')
                MIS_ID <- c('None')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],' + ',PRI_VAR[2],
                                                ' + animal.key.anirandgroup'))
        }else if(TISSUE == 'Ovaries'){
                TIS <- c('OVR')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(TISSUE == 'Testes'){
                TIS <- c('TES')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90017016302_SN2','90031016302_SN2','90127016302_SN2','90015016302_SN2','90129016302_SN2')
        }else if(TISSUE == 'Spleen'){
                TIS <- 'SPL'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }
        # Create the row for the dataframe
        df_r <- data.frame(Tissue = TISSUE,
                           Tis = TIS,
                           Seq_Batch = 'None',
                           Misidentified = MIS_ID,
                           Adjusted_Variance = PRI_VAR,
                           Outliers = OUTLIERS,
                           Formula = as.character(FINAL_FORMULA))
        df_tbl <- rbind(df_tbl,df_r)
}

# Save the decision table
table_file <- paste0(WD,'/data/20200603_rnaseq-tissue-data-assambly-table_steep.txt')
#write.table(df_tbl, file = table_file,sep = '\t',row.names = F,quote = F)
#' #### Final table saved as:
#' `r table_file`
#'
#'
#'
#' ## Take Homes
#' #### Tissue:
#' `r TISSUE`
#'
#' #### Mis-identified Samples:
#' `r MIS_ID`
#'
#' #### Covariables to Investigate:
#' `r PRI_VAR`
#'
#' #### Final Modeling Formula:
#' `r FINAL_FORMULA`
#'
#' #### Outlier Samples:
#' `r OUTLIERS`
#'
#'
#'
#'
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
                col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% parse_time() %>% as.numeric()
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
        
        # Generate a time-last-fed variable
        col_data <- col_data %>%
                mutate(animal_time_last_fed = case_when(
                        animal.key.anirandgroup %!in% c('Control - 7 hr', 'Exercise - 7 hr') ~ parse_time('8:00'),
                        (animal.key.anirandgroup %in% c('Control - 7 hr') &
                                 animal.registration.sex == 'Male') ~ parse_time('11:50'),
                        (animal.key.anirandgroup %in% c('Control - 7 hr') &
                                 animal.registration.sex == 'Female') ~ parse_time('12:30'),
                        (animal.key.anirandgroup %in% c('Exercise - 7 hr') &
                                 animal.registration.sex == 'Male') ~ parse_time('12:30'),
                        (animal.key.anirandgroup %in% c('Exercise - 7 hr') &
                                 animal.registration.sex == 'Female') ~ parse_time('12:50')) %>% as.numeric())
        
        # Generate a time-fasted variable
        col_data$calculated.variables.deathtime_after_fed <- (col_data$specimen.collection.t_death - col_data$animal_time_last_fed) %>% as.numeric()
        
        # Generate a non-zero variable
        # The fraction of non-zero counts (function)
        nonzero <- function(x) sum(x != 0)
        # Non zero counts
        nzc <- apply(count_data,2,nonzero)
        ginic <- apply(count_data,2,reldist::gini)
        # Non zero fraction
        nzf <- nzc/nrow(count_data)
        # Mean gene expresison
        m <- apply(count_data,2,mean)
        # SD gene expression
        sd <- apply(count_data,2,sd)
        # Coefficient of variation
        cv <- sd/m
        # Join the dataframe into meta data
        nz_tis_df <- data.frame('sample_key' = colnames(count_data),
                                'NZF' = nzf,
                                'GINI' = ginic,
                                'CV' = cv)
        col_data <- left_join(col_data, nz_tis_df, by = c('sample_key'))
        
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

meta_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
#' #### Polished metadata saved as:
#' `r meta_file`
#'
count_file <- paste0(WD, '/data/20200603_rnaseq-counts-pass1a-stanford-sinai-processed_steep.rds')
#' #### Polished read counts saved as:
#' `r count_file`

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
#' Data from Supplementary Table 2 from:
#' Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#'
#' ## Place Genes in Genomic Ranges
#'
#' #### Reference Genome and Annotation:
#' Rnor_6.0 (GCA_000001895.4) assembly from Ensembl database (Release 96)
#' Found at: http://uswest.ensembl.org/Rattus_norvegicus/Info/Index.
#'
#' FASTA: Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
#' ftp://ftp.ensembl.org/pub/release-96/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
#'
#' GTF: Rattus_norvegicus.Rnor_6.0.96.gtf.gz
#' ftp://ftp.ensembl.org/pub/release-96/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.96.gtf.gz
#'
#' ## Annotate Genes by Chromosome

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
#seqlevels(Rn_TxDb)[1:23]
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
if(TISSUE == c('Gastrocnemius_MSSM_1')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex))
}else if(TISSUE == c('Gastrocnemius_Stanford_1')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'Stanford_1') %>%
                filter(!is.na(animal.registration.sex))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex))
}
rownames(tod_cols) <- tod_cols$sample_key

# Collect samples without NA values in TOD
nona_sams <- tod_cols %>%
        filter(!is.na(specimen.collection.t_death_hour)) %>%
        filter(!is.na(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Collect tissue specific counts
tod_counts <- count_data[,nona_sams]

#' ##### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
design = ~ 1 # Primary variable needs to be last.
title = paste0('Design: ',as.character(design))
dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = design)
# Reasoning from:
#citation('PROPER')
#dds
#' #### We remove genes with an average sequencing depth (across samples) of less than 1
#' Summary of DESeqDataSet Before Filtering:
dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>%
        nrow() %>% as.character()
reads_n <- 1
keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
dds2 <- dds1[keep,]
#' Summary of counts and annotation data in a DESeqDataSet after filtering
# TODO: Critic from Jun: Here we are removing features that have a low average expression. This may be removing important features that might have zero counts in some samples and higher counts in specific groups. Consider developing an algorithm that will account for features with expression in n or more samples.
dds2
filter_n <- nrow(dds1) - nrow(dds2) - as.numeric(zero_n)
filter_p <- filter_n/(nrow(dds1) - as.numeric(zero_n))
total_n <- nrow(dds1) - nrow(dds2)

#' ##### Number of Genes Removed:
#' Number of genes with average counts (across samples) less than 1 is `r zero_n`
#'
#' ##### If more stringent filtering thresholds employed:
#' Removing reads with less than or equal to `r reads_n` removes an additional `r filter_n` features or removes `r filter_p*100`% of the non-zero reads (total of `r total_n` features removed).
dds <- dds2

#' #### Reads per million for each sample
sort(colSums(assay(dds)))/1e6

#' #### Size Factor Estimates (summary across samples):
#' Size facotrs are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
# estimateSizeFactors() gives us a robust estimate in sequencing depth.
dds <- estimateSizeFactors(dds)
summary(sizeFactors(dds))

#' #### Transform:
rld <- DESeq2::vst(dds, blind = FALSE)
#' Variance Stabalizing Transform (vst)
# for(n in 1){
#         start_time <- Sys.time()
#         #rld <- DESeq2::rlog(dds)
#         end_time <- Sys.time()
#         print(end_time - start_time)
# }
# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

#' ## Early Identification/Removal of Outliers
#' Certain Tissues (e.g. Hippocampus & Testes) demonstrate obvious outlier samples that disrupt covariant identification and correction (could not be "batch-corrected"); therefore, we exercise an early identification and removal of outliers.

#+ Early Identification/Removal of Outliers, results="asis",echo=FALSE,message=FALSE
################################################################################
################### Early Identification/Removal of Outliers ###################
################################################################################

#' #### Outlier Sample Detection Strategy (PCA-Approach)
#' ###### Samples are plotted across principle components (PCs) with PC values used to calculate distances across multiple dimensions. Distance metrics include: Infinite (median and MAD), Mahalanobis, and Local Outlier Factor (LOF). Samples are classified as outliers by their distances across one or more "high variance" PCs. Often, a threshold of greater than or equal to 6 standard deviations (Tukey) from the mean distances in one or more PC(s) justifies outlier status, however, classification is at the discretion of the analysis. Thresholds are adjusted in consideration of skewed distance distributions and multiple tests. This strategy was adopted from Florian Privé, PhD (https://www.r-bloggers.com/detecting-outlier-samples-in-pca/).
#'
#' #### Description of Plots
#' ###### Variance by PCs:  A plot of the variance explained by PCs. Dashed vertical lines represent the elbow.
#' ###### Samples by Multidimensional Distance (2 plots):  Samples are plotted by their Local Outlier Factor (log10) and Mahalanobis Distance (log10). Multivariant distances apply across PCs and the PC (above the elbow in the previous plot) with the maximum infinite distance--unique to each PC--is plotted. Thresholds in red represent a potential outlier detection threshold accounting for skewness and multiple testing. In the PC-annotated plot, samples are labeled if they appear above either threshold.
#' ###### PCA of Testes (PCs 1 & 2):  Samples plotted across major PCs 1 & 2. Suspected outliers--samples above either threshold in the Multidimesnional Distance Plots--are labeled by a unique sample identifier. All samples are colored by their maximum infinite distance within major PCs above the "elbow" from plot, "Variance by PCs."
#'

# A few Hippocampus samples demonstrate considerably low Lib_molarity that cannot be compensated for. These samples will be treated as outliers and removed.
if(TISSUE %in% c('Hippocampus','Testes')){
        if(TISSUE == 'Hippocampus'){
                pri_var <- 'Lib_molarity'
        }else if(TISSUE == 'Testes'){
                pri_var <- 'Lib_frag_size..bp.'
                #pri_var <- 'Lib_molarity..nM.'
        }
        pca <- prcomp(t(assay(rld)), scale. = F)
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        PC <- pca$x %>% as.data.frame()
        PC$sample_key <- row.names(pca$x)
        df_plot <- data.frame(PC = seq_along(pca$sdev^2/sum(pca$sdev^2)),
                              VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
                mutate(ELBOW = ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE'))
        p <- df_plot %>%
                ggplot(aes(PC,VAR, color = ELBOW)) +
                geom_point() +
                geom_vline(xintercept = elbow_finder(df_plot$PC,df_plot$VAR)[1], linetype='dashed') +
                ggtitle('Variance by PCs') +
                ylab('Variance')
        plot(p)
        PC_meta <- t(pca$x) %>%
                as_tibble() %>%
                mutate(PC = seq_along(pca$sdev^2/sum(pca$sdev^2))) %>%
                mutate(VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
                mutate(ELBOW = ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE')) %>%
                filter(ELBOW == 'UPPER') %>%
                mutate(PC = as.character(sprintf('PC%d',PC))) %>%
                filter(PC %in% c('PC1','PC2','PC3','PC4'))
        # Reset the pca object
        pca <- prcomp(t(assay(rld)), scale. = F,
                      rank. = elbow_finder(df_plot$PC,df_plot$VAR)[1])
        # Collect Standard Deviations
        pca$sd <- apply(pca$x, 2, function(x) abs(x - median(x)) / mad(x))
        join_meta <- PC_meta[,c('PC','VAR')]
        names(join_meta) <- c('SD_DIST_PC','PC_VAR')
        PC <- PC %>%
                select(
                        c(sprintf('PC%d',1:(elbow_finder(df_plot$PC,df_plot$VAR)[1])),'sample_key')) %>%
                mutate(sample_n = seq_along(PC1)) %>%
                mutate(SD_DIST = apply(pca$x, 2, function(x) abs(x - median(x)) / mad(x)) %>%
                               apply(1, max)) %>%
                mutate(SD_DIST_PC = colnames(pca$sd)[apply(pca$sd,1,which.max)]) %>%
                mutate(TUKEY_OUTLIER = ifelse(SD_DIST > 6, 'Outlier', 'Normal')) %>%
                mutate(MAHALANOBIS_DIST = dist_ogk(pca$x)) %>%
                # Bonferroni correction
                mutate(MAHALANOBIS_FDR = pchisq(MAHALANOBIS_DIST, df = elbow_finder(df_plot$PC,df_plot$VAR)[1], lower.tail = FALSE)*length(MAHALANOBIS_DIST)) %>%
                mutate(MAHALANOBIS_OUTLIER = ifelse(MAHALANOBIS_FDR <= 0.05,
                                                    'Outlier', 'Normal')) %>%
                mutate(LOF = LOF(pca$x, log = F)) %>%
                left_join(y = col_data, by = 'sample_key') %>%
                left_join(y = join_meta, by = 'SD_DIST_PC') %>%
                arrange(dplyr::desc(SD_DIST))
        p <- PC %>%
                mutate(PC_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) |
                                                 log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), SD_DIST_PC, '')) %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = animal.key.anirandgroup)) +
                geom_point(size = 2) +
                scale_color_manual(values=ec_colors) +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle('Samples by Multidimensional Distance') +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)')
        # Samples plotted against 2 multidimension distance metrics--local
        plot(p)
        # Examine the variance from of the PC that is associated with the outlier
        p <- PC %>%
                mutate(PC_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) |
                                                 log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), SD_DIST_PC, '')) %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = PC_VAR)) +
                geom_point(size = 2) +
                geom_label_repel(aes(label=PC_LABEL)) +
                #geom_label_repel(aes(label=sample_key)) +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle('Samples by Multidimensional Distance') +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)') +
                labs(color='Variance per PC')
        plot(p)
        # Visualize the gini index
        p <- PC %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = GINI)) +
                geom_point(size = 2) +
                #geom_label_repel(aes(label=sample_key)) +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle('Samples by Multidimensional Distance') +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)') +
                labs(color='Gini Index') +
                scale_color_continuous(high = "#132B43", low = "#56B1F7",
                                       name = "Gini Index")
        plot(p)
        
        PC <- PC %>%
                mutate(OUT_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) |
                                                 log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, ''))
        # Plot the Gini Index for these samples
        p <- ggplot() +
                geom_violin(data = PC, aes(x = TISSUE, y= GINI)) +
                geom_jitter(data = PC, aes(x = TISSUE, y= GINI,
                                           color = log(MAHALANOBIS_DIST)),
                            height = NULL, width = 0.15,
                            alpha = 0.9, size = 2) +
                geom_label_repel(data = PC, aes(x = TISSUE, y= GINI,label=OUT_LABEL)) +
                xlab("") +
                ylab("Gini Index") +
                theme_light() +
                theme(axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      axis.title.y = element_text(size = 20),
                      legend.title = element_text(),
                      legend.text = element_text(size = 16)) +
                scale_color_continuous(high = "#132B43", low = "#56B1F7",
                                       name = "Mahalanobis Distance (log)")
        plot(p)
        
        p <- PC %>%
                mutate(OUTLIER = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) |
                                                log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, '')) %>%
                ggplot(aes(x = PC1, y = PC2, color = log(MAHALANOBIS_DIST))) +
                geom_point(size = 3) +
                coord_equal() +
                geom_label_repel(aes(label=OUTLIER)) +
                xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
                ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) +
                ggtitle(paste0(TISSUE)) +
                scale_color_continuous(high = "#132B43",low = "#56B1F7",
                                       name = "Mahalanobis Distance (log)")
        plot(p)
}

#' ## Qualitative Assessment of Variance in Tissue

# Qualitative Assessment of Variance in Tissue
################################################################################
#####     Qualitative Assessment of Variance in Tissue      ####################
################################################################################

# Examine the Annotations most strongly associated with PCs
################################################################################

# Find variables associated with PCs
pc_cor_df <- cor_PC_1_6(rld = rld, ntop = 20000, intgroups = names(tod_cols)) %>%
  filter(Adjusted_R_Sq > 0.2) %>%
  arrange(PC)
# Examine Variables of interest for Pcs 1 - 4
pc_cor_df %>%
  filter(PC %in% c(1:4))
# # Quick Test
# animal.key.batch %in% (pc_cor_df %>% filter(PC %in% c(1,2)) %>% select(Condition) %>% unlist() %>% as.character() %>% unique())

# Examine strength of PC
pca <- prcomp(t(assay(rld)), scale. = F)
df_plot <- data.frame(PC = seq_along(pca$sdev^2/sum(pca$sdev^2)),
                      VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
  mutate(ELBOW = ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE'))
p <- df_plot %>%  
  ggplot(aes(PC,VAR, color = ELBOW)) +
  geom_point() +
  geom_vline(xintercept = elbow_finder(df_plot$PC,df_plot$VAR)[1], linetype='dashed')
plot(p)

# colData(rld) %>% 
#   as.data.frame() %>% select(animal.key.anirandgroup, 
#                         animal.registration.sex,
#                         animal_time_last_fed)

if(T){
# Visualize variables in PC's 1 - 2
voi12 <- pc_cor_df %>%
  filter(PC %in% c(1:2)) %>%
  select(Condition) %>% unlist() %>% as.character %>% unique()
# Plot the top variables associated with PCs 1 or 2
for( v in voi12){
  # Quick investigation of variables
  title_parts <- pc_cor_df %>%
    filter(Condition == v)
  title <- paste0(title_parts$Condition,'\nPC = ',
                  title_parts$PC,'; R2 = ',
                  round(title_parts$Adjusted_R_Sq,3))
  p <- DESeq2::plotPCA(rld, intgroup =v, ntop = 20000) +
    ggtitle(title) +
    coord_equal()
  plot(p)
}
}

if(T){
  # Visualize variables in PC's 3 - 4
  voi34 <- pc_cor_df %>%
    filter(PC %in% c(3:4)) %>%
    select(Condition) %>% unlist() %>% as.character %>% unique()
  
  # Prepare the data for PC3 and PC4 visualization
  pca <- prcomp(t(assay(rld)), scale. = F)
  df_plot <- pca$x %>% as.data.frame()
  df_plot$sample_key <- row.names(df_plot)
  df_plot <- df_plot %>% 
    left_join(y = col_data, by = 'sample_key')
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  # Plot the top variables associated with PCs 1 or 2
  for( v in voi34){
    # Quick investigation of variables
    title_parts <- pc_cor_df %>%
      filter(Condition == v)
    title <- paste0(title_parts$Condition,'\nPC = ',
                    title_parts$PC,'; R2 = ',
                    round(title_parts$Adjusted_R_Sq,3))
    p <- df_plot %>%
      ggplot(aes(x = PC3, y = PC4, color = !!as.symbol(v))) +
      geom_point(size = 2) +
      coord_equal() +
      xlab(paste0('PC3: ',round(percentVar[3]*100,2),'% variance')) +
      ylab(paste0('PC4: ',round(percentVar[4]*100,2),'% variance')) + 
      ggtitle(title)
    plot(p)
  }
}

################################################################################
###### Manual Annotation Opportunity: Fill in MIS_ID_VAR variable ##############
################################################################################

# Examine a labled PCA Plot
################################################################################
if(TISSUE %!in% c('Ovaries','Testes')){
  for(mis_id_var in MIS_ID_VAR){
    if(MIS_ID_VAR == 'None'){
      mis_id_var <- 'animal.registration.sex'
      }
  #mis_id_var <- MIS_ID_VAR[1]
  print(mis_id_var)
# Examine metadata varibales most correlated with suspected misidentified samples to determine if there is an effect are genuine
cor_out_df <- cor_outlier2(rld = rld, ntop = 20000, intgroups = names(colData(rld)))
cor_out_df <- cor_out_df %>% 
  filter(Condition != 'outlier')
cor_out_df

# Show a PCA plot labeled with suspected misidentified samples
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c(mis_id_var, 'sample_key','animal.key.anirandgroup',
                                      'animal.registration.sex','animal.registration.cagenumber'), 
                           returnData=TRUE, ntop = 20000)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
#pdf(paste0(WD,'/plots/20200505_rnaseq-',TIS,'-PCA-naive-modeling_steep.pdf'),
#    width = 6, height = 4)
p <- pcaData %>%
  mutate(sample_key = as.character(sample_key)) %>%
  mutate(sample_misid_suspect = ifelse(sample_key %in% MIS_SUS, sample_key, '')) %>%
  ggplot(aes_string('PC1', 'PC2', color=mis_id_var)) +
  geom_point(size=3) +
  #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
  geom_label_repel(aes(label=sample_misid_suspect),hjust=0, vjust=0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) + 
  #coord_fixed() +
  ggtitle(paste0('PCA of ',TISSUE,' Gene Expression:\nNaive Model (~ 1)')) +
  guides(color=guide_legend(title=mis_id_var)) +
  theme(legend.title=element_blank())
plot(p)
#dev.off()
# Examine the PCA annotated with variables of interest
p <- pcaData %>%
  mutate(sample_key = as.character(sample_key)) %>%
  mutate(sample_misid_suspect = ifelse(sample_key %in% MIS_SUS, sample_key, '')) %>%
  ggplot(aes(PC1, PC2, color=animal.key.anirandgroup)) +
  geom_point(size=3) +
  #geom_label_repel(aes(label=!!as.symbol(mis_id_var)),hjust=0, vjust=0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  #coord_fixed() +
  ggtitle(paste0('PCA of ',TISSUE,' Gene Expression')) +
  guides(color=guide_legend(title='animal.key.anirandgroup')) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())
plot(p)

# Plot Sex-Specific Plots
################################################################################
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

#### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.
chryexp <- colMeans(assay(rld)[Y_genes,])
chryexp_df <- data.frame('counts' = chryexp)
chryexp_df$sample <- row.names(chryexp_df) %>% as.factor()
chryexp_df <- chryexp_df %>%
        mutate(sex = ifelse(sample %in% male_tis, 'Male', 'Female'))

# Visualize genes on the Y Chromosome
##############################################
##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution.
mypar()
p <- chryexp_df %>%
  mutate(plot_labs = ifelse(sample %in% MIS_ID, as.character(sample), '')) %>%
  ggplot(aes(counts, colour = sex)) +
  geom_freqpoly(alpha = 0.6) +
  xlab('Median counts on Y genes (normalized)') +
  ylab('Frequency') +
  ggtitle('Frequency Polyplot & Sample Dotplot: \nSamples plotted by Y gene coverage (median)') +
  geom_dotplot(aes(fill = sex), alpha = 0.5) +
  #geom_text_repel(aes(label = sample), y = 2) 
  geom_text_repel(aes(label = plot_labs), box.padding = unit(2, 'lines'), y = 0)
plot(p)

################################################################################
###### Manual Annotation Opportunity: Fill in MIS_SUS variable #################
################################################################################

# ## Sample Classification: High Dimension (Gene) Outlier Prediction + Supervised

# Sample Classification: High Dimension Outlier Prediction + Supervised 
################################################################################
##     Sample Classification: High Dimension Outlier Prediction + Supervised #
################################################################################

# (iii) standardization of the data so that the expression measures for each array have mean 0 and variance 1 across genes.
edata <- scale(assay(rld))

# check that we get mean of 0 and sd of 1
apply(edata , 2, mean)
apply(edata , 2, sd)

 ## Reassign Annotation Groups to MisIdentified Samples: Outlier Detection (Per Gene)

# Reassign Annotation Groups to MisIdentified Samples
################################################################################
#####     Reassign Annotation Groups to MisIdentified Samples      #############
################################################################################
# TODO: Add a conditional statement if samples were misidentieid or not

# identify multivariate outliers
x.out <- pcout(edata, makeplot = T)

# Create a dataframe with 3 columns: gene, weight_combined, gene_index
df_out <- data.frame(names(x.out$wfinal), x.out$wfinal) %>% as_tibble()
names(df_out) <- c('gene', 'weight_combined')
df_out$index <- row.names(df_out)

# Arrange the dataframe by weight
df_out <- df_out %>%
  arrange(weight_combined)

# Plot the Outliers by combined weight
p <- df_out %>%
  sample_n(500) %>%
  ggplot(aes(index, weight_combined)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0.1) +
  ggtitle('Weights of 500 randomly selected genes')
plot(p)

# Number of strong outlier genes
genes_out <- df_out %>%
  filter(weight_combined < 0.10) %>%
  select(gene) %>% unlist() %>% as.character()

# Unsupervised Model
################################################################################
set.seed(100)
rld_out <- rld[names(rld) %in% genes_out,]
df <- t(assay(rld_out)) %>% as.data.frame()
df$sample_key <- row.names(df)
join_df <- col_data %>%
  dplyr::select(animal.registration.sex, sample_key, animal.key.anirandgroup)
df <- left_join(df,join_df, by = 'sample_key')
row.names(df) <- make.names(df$animal.registration.sex, unique = TRUE)
#df <- df %>%
#  filter(animal.key.anirandgroup %!in% c('Control - IPE', 'Control - 7 #hr',
#                                         'Exercise - 24 hr', 'Exercise #- 48 hr',
#                                         'Exercise - 7 hr'))
d <- dist( df )
hc <- hclust(d)
sex <- as.character(df$animal.registration.sex)
sample <- as.character(df$sample_key)
ex_grp <- as.character(df$animal.key.anirandgroup)
myplclust(hc, labels=sex, lab.col=as.fumeric(sex), cex=1)
myplclust(hc, labels=sample, lab.col=as.fumeric(sex), cex=1)

# Supervised Model
################################################################################
# Classification Model, Predict Class with outlier gene set
set.seed(100)
rld_out <- rld[names(rld) %in% genes_out,]
# Remove the Sample in question
rld_mis <- rld_out[,rld_out$sample_key %in% MIS_SUS]
rld_out <- rld_out[,rld_out$sample_key %!in% MIS_SUS]

# Create a train and test set
rld_out[[mis_id_var]] <- as.character(rld_out[[mis_id_var]])
rld_out[[mis_id_var]] <- as.factor(rld_out[[mis_id_var]])
index = createDataPartition(y=rld_out[[mis_id_var]], p=0.7, list=FALSE)
sample_train <- as.character(rld_out$sample_key[index])
sample_test <- as.character(rld_out$sample_key[-index])
rld_train = rld_out[,rld_out$sample_key %in% sample_train]
rld_test = rld_out[,rld_out$sample_key %in% sample_test]
# Reorganize the train and test DFs
join_df <- tod_cols %>%
  select('sample_key', all_of(mis_id_var))
# Train
train_df <- data.frame(t(assay(rld_train)))
train_df$sample_key <- row.names(train_df)
train_df <- train_df %>%
  left_join(join_df, by = c('sample_key')) %>%
  as_tibble() %>%
  select(-sample_key)
# Test
test_df <- data.frame(t(assay(rld_test)))
test_df$sample_key <- row.names(test_df)
test_df <- test_df %>%
  left_join(join_df, by = c('sample_key')) %>%
  as_tibble() %>%
  select(-sample_key)
# Misidentified
mis_df <- data.frame(t(assay(rld_mis)))
mis_df$sample_key <- row.names(mis_df)
mis_df <- mis_df %>%
  left_join(join_df, by = c('sample_key')) %>%
  as_tibble() %>%
  select(-sample_key)
# Double check (must be true)
nrow(colData(rld)) == nrow(train_df) + nrow(test_df) + nrow(mis_df)

# genes to remove from model
x = nearZeroVar(train_df, saveMetrics = TRUE)
str(x, vec.len=2)
# The Zero variance predictors
x[x[,'zeroVar'] > 0, ]
# The near Zero variance predictors
x[x[,'zeroVar'] + x[,'nzv'] > 0, ]

# prepare training scheme
control <- trainControl(method='repeatedcv', number = 10, repeats = 10)
control <- trainControl(method='cv')
#control <- trainControl(method='LOOCV')
design <- formula(paste0(mis_id_var,' ~ .'))

# Random Forest
######################
# Set up the random forest tuning grid
# TODO: Automate mtry selection
grid_rf <- expand.grid(mtry = c(50,100,250,500,750,1000))
set.seed(100)
design <- formula(paste0(mis_id_var,' ~ .'))
rf.fit <- train(design, data=train_df,
                method = 'rf',
                ntree = 1000,
                trControl = control,
                tuneGrid = grid_rf)

# Examine model
rf.fit

# Visualize feature importance
rfImp <- varImp(rf.fit, scale = FALSE)
plot(rfImp, top = 25, scales = list(y = list(cex = .95)))
# Annotate gene ids with symbols
row.names(rfImp$importance) <- mapIds(org.Rn.eg.db, row.names(rfImp$importance), 
                                      'SYMBOL', 'ENSEMBL') %>% 
  make.names(unique=TRUE)
plot(rfImp, top = 25, scales = list(y = list(cex = .95)))

pred_sex = predict(rf.fit, test_df)
table(pred_sex, test_df[[mis_id_var]])
pred.accuracy = round(mean(pred_sex == test_df[[mis_id_var]])*100,2)
pred.accuracy

# For sake of conditions below
correct_df <- mis_df
if(nrow(mis_df) >= 1){
  true_var = predict(rf.fit, mis_df)
  true_df <- col_data %>%
    filter(sample_key %in% MIS_SUS) %>%
    select('sample_key', all_of(mis_id_var))
  names(true_df) <- c('Recorded_sample_key', 'Recorded_Var')
  true_df$Corrected_Var <- true_var
  if(TISSUE %in% c('Hypothalamus')){
    true_df$Corrected_sample_key <- c(as.character(true_df$Recorded_sample_key[2]),
                                      as.character(true_df$Recorded_sample_key[1]))
  }
  # Remove samples that were correctly identified
  if(nrow(true_df %>%
          filter(Recorded_Var != Corrected_Var)) >= 1){
    correct_df <- true_df %>%
      filter(Recorded_Var != Corrected_Var) %>%
      mutate(Corrected_sample_key = as.character(Corrected_sample_key))
  }
}

################################################################################
###### Manual Annotation Opportunity: Fill in MIS_ID variable ##################
################################################################################

# TODO: Improve this adjustment of metadata to be independent of mis_id_var class
# TODO: Adjust the vial label
if(TISSUE %in% c('Hypothalamus')){
    for(i in 1:nrow(correct_df)){
      col_data <- col_data %>%
        mutate(sample_key = 
                 ifelse(
                   vial_label == (correct_df[i,]$Recorded_sample_key %>% 
                                    as.character() %>% 
                                    strsplit(split = '_') %>% 
                                    unlist())[1],
                   as.character(correct_df[i,]$Corrected_sample_key),
                   as.character(sample_key)
                 ))
  }
}

# Save the updated metadata
################################################################################
if(nrow(true_df %>%
        filter(Recorded_Var != Corrected_Var)) >= 1){
  # # Meta Data
  meta_file <- paste0(WD,'/data/20200603_rnaseq-meta-pass1a-stanford-sinai-proc_steep.rds')
  saveRDS(col_data, file = meta_file)
  
  # Update the rld object
  ################################################################################
  tod_cols <- col_data %>%
    filter(Tissue == TISSUE) %>%
    filter(!is.na(animal.registration.sex))
  rownames(tod_cols) <- tod_cols$sample_key
  nona_sams <- tod_cols %>%
    filter(!is.na(specimen.collection.t_death_hour)) %>%
    filter(!is.na(animal.registration.sex)) %>%
    select(sample_key) %>% unlist() %>% as.character()
  tod_counts <- count_data[,nona_sams]
  all(rownames(tod_cols) == colnames(tod_counts))
  dds <- DESeqDataSetFromMatrix(countData = tod_counts,
                                      colData = tod_cols,
                                      design = formula('~ 1'))
  dds <- dds[rowSums(counts(dds))/ncol(dds) >= 1,]
  dds <- estimateSizeFactors(dds)
  rld <- DESeq2::vst(dds, blind = F)
}
}
}

################################################################################
#########     Manual Annotation Opportunity: Assign PRI_VAR      ###############
################################################################################

#' ## Correct for Variance: Examine Genes Driving Primary Variance of Interest

#+ Examine Genes Driving Primary Variance of Interest
################################################################################
#########     Examine Genes Driving Primary Variance of Interest      ##########
################################################################################

for(pri_var in PRI_VAR){
  if(PRI_VAR == 'None'){
    pri_var <- 'animal.registration.sex'
  } 
  #pri_var <- PRI_VAR[1]
  print(pri_var)
  # Find variables associated with PCs
  pc_cor_df <- cor_PC_1_6(rld = rld, ntop = 20000, intgroups = names(colData(rld))) %>%
    filter(Adjusted_R_Sq > 0.2) %>%
    arrange(PC)
  # Examine Variables of interest for Pcs 1 and 2
  p <- pc_cor_df %>%
    filter(PC %in% c(1:4))
  print(p)
  
  if(T){
  # Visualize variables in PC's 1 and 2
  voi <- pc_cor_df %>%
    filter(PC %in% c(1,2)) %>%
    select(Condition) %>% unlist() %>% as.character %>% unique()
  # Plot the top variables associated with PCs 1 or 2
  for( v in voi){
    #v <- 'RNA_extr_conc'
    # Quick investigation of variables
    title_parts <- pc_cor_df %>%
      filter(Condition == v)
    title <- paste0(title_parts$Condition,'\nPC = ',
                    title_parts$PC,'; R2 = ',
                    round(title_parts$Adjusted_R_Sq,3))
    p <- DESeq2::plotPCA(rld, intgroup =v, ntop = 20000) +
      ggtitle(title)
    plot(p)
  }
  }
  if(T){
  # Visualize variables in PC's 3 - 4
  voi34 <- pc_cor_df %>%
    filter(PC %in% c(3:4)) %>%
    select(Condition) %>% unlist() %>% as.character %>% unique()
  
  # Prepare the data for PC3 and PC4 visualization
  pca <- prcomp(t(assay(rld)), scale. = F)
  df_plot <- pca$x %>% as.data.frame()
  df_plot$sample_key <- row.names(df_plot)
  df_plot <- df_plot %>% 
    left_join(y = col_data, by = 'sample_key')
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  # Plot the top variables associated with PCs 1 or 2
  for( v in voi34){
    # Quick investigation of variables
    title_parts <- pc_cor_df %>%
      filter(Condition == v)
    title <- paste0(title_parts$Condition,'\nPC = ',
                    title_parts$PC,'; R2 = ',
                    round(title_parts$Adjusted_R_Sq,3))
    p <- df_plot %>%
      ggplot(aes(x = PC3, y = PC4, color = !!as.symbol(v))) +
      geom_point(size = 2) +
      coord_equal() +
      xlab(paste0('PC3: ',round(percentVar[3]*100,2),'% variance')) +
      ylab(paste0('PC4: ',round(percentVar[4]*100,2),'% variance')) + 
      ggtitle(title)
    plot(p)
  }
  }
  
  # If PC1 Unknown, then create a pho PC1 variable (binary), perform DE + Pathway analysis to investigate
  ######################################
  # Make this if statement true if you'd like to examine a PC of interest that is not obviously correlated with a meta data variable of interest
  if(TISSUE %in% c('Aorta','Ovaries','Testes')){
    # Fill in manually
    PC_THRESH <- 0
    pc_int <- 'PC1'
    #all(tod_cols$sample_key == row.names((prcomp(t(assay(rld)), scale. = F))$x))
    tod_cols$PC_INT <- (prcomp(t(assay(rld)), scale. = F))$x %>%
      as.data.frame() %>% select(all_of(pc_int)) %>% unlist() %>% as.numeric()
    tod_cols <- tod_cols %>%
      mutate(PC_UNKN = factor(ifelse(PC_INT < PC_THRESH, '0', '1'), levels = c('0','1')))
    pri_var <- 'PC_UNKN'
    # Examine PCA
    df_plot <- (prcomp(t(assay(rld)), scale. = F))$x %>%
      as.data.frame()
    df_plot$sample_key <- row.names(df_plot)
    df_plot <- df_plot %>% 
      left_join(y = col_data, by = 'sample_key')
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    p <- df_plot %>%
      ggplot(aes_string(x = pc_int, y = 'PC2', color = 'animal.registration.sex')) +
      geom_point(size = 2) +
      coord_equal() +
      xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
      ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance'))
    plot(p)
    pri_var <- PRI_VAR[1]
    if(pri_var == 'None'){pri_var <- 'animal.registration.sex'}
  }
  
  # DE Test for Variable (Ensure factors are properly ordered)
  ( design = formula(paste0('~ 1')) ) # Primary variable needs to be last.
  dds_vp <- DESeqDataSetFromMatrix(countData = tod_counts,
                                   colData = tod_cols,
                                   design = design)
  dds_vp <- dds_vp[rowSums(counts(dds_vp))/ncol(dds_vp) >= 1,]
  dds_vp <- estimateSizeFactors(dds_vp)
  res_vp <- DESeq(dds_vp) %>%
    results(alpha = 0.05,lfcThreshold=0.5)
  
  if(pri_var == 'animal.registration.sex') {
    # Adjust the result dataframe
    res_vp <- res_vp %>%
      as.data.frame() %>%
      mutate(ENSEMBL_RAT = row.names(res_vp)) %>%
      mutate(SIG = ifelse((abs(log2FoldChange) >= 0.5 & padj <= 0.05), 'SIG','NS')) %>%
      mutate(CHROM = factor(case_when(ENSEMBL_RAT %in% Y_ens_id ~ 'Y',
                                      ENSEMBL_RAT %in% X_ens_id ~ 'X',
                                      ENSEMBL_RAT %!in% c(X_ens_id,Y_ens_id) ~ 'AUTO'),
                            levels = c('X', 'Y', 'AUTO'))) %>%
      arrange(padj) %>%
      as_tibble()
    res_vp 
    top20 <- res_vp %>%
      head(n = 20) %>%
      select(ENSEMBL_RAT) %>% unlist() %>% as.character()
    n_sig <- res_vp %>% filter(SIG == 'SIG') %>% nrow()
    # Generate a volcano plot
    p <- res_vp %>%  
      ggplot(aes(log2FoldChange, -log10(padj), color = CHROM, size = CHROM)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = -log10(0.05), linetype='dashed') +
      geom_vline(xintercept = 0.5, linetype='dashed') +
      geom_vline(xintercept = -0.5, linetype='dashed') +
      ggtitle('DE Genes from Females to Males:
Adjusted p-value <= 0.05; |Log2FC| >= 0.5') +
      ylab(bquote('-'~Log[10]~ 'Adjusted p-value (Bonferroni)')) +
      scale_color_manual(values=c('blue', 'red', 'grey')) +
      scale_size_manual(values=c(3, 3, 2)) +
      theme_bw()
    plot(p)
    # Generate an accompanying pvalue histogram
    p <- res_vp %>%
      ggplot(aes(x = padj, fill = CHROM)) +
      #geom_freqpoly() +
      ggtitle('DE Genes Females to Males (adjusted p-value <= 0.05)') +
      geom_dotplot(aes(x = padj), dotsize = 0.5, alpha = 0.5) +
      geom_vline(xintercept = 0.05, linetype='dashed') +
      xlab('Adjusted p-value (Bonferroni)') +
      ylab('Density') +
      scale_fill_manual(values=c('blue', 'red', 'grey')) +
      theme_bw()
    plot(p)
  } else {
    res_vp <- res_vp %>%
      as.data.frame() %>%
      mutate(ENSEMBL_RAT = row.names(res_vp)) %>%
      mutate(SIG = ifelse((abs(log2FoldChange) >= 0.5 & padj <= 0.05), 'SIG','NS')) %>%
      arrange(padj) %>%
      as_tibble()
    res_vp 
    top20 <- res_vp %>%
      head(n = 20) %>%
      select(ENSEMBL_RAT) %>% unlist() %>% as.character()
    n_sig <- res_vp %>% filter(SIG == 'SIG') %>% nrow()
    # Generate a volcano plot
    p <- res_vp %>%
      mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, 'SYMBOL', 'ENSEMBL')) %>%
      #mutate(plot_label = ifelse(ENSEMBL_RAT %in% top20, SYMBOL_RAT, '')) %>%
      mutate(plot_label = ifelse(ENSEMBL_RAT %in% top20, ENSEMBL_RAT, '')) %>%
      ggplot(aes(log2FoldChange, -log10(padj), color = SIG)) +
      geom_point(alpha = 1) +
      #geom_label_repel(aes(label=plot_label),hjust=0, vjust=0) +
      #geom_text_repel(aes(label = plot_label), box.padding = unit(2, 'lines')) +
      geom_hline(yintercept = -log10(0.05), linetype='dashed') +
      geom_vline(xintercept = 0.5, linetype='dashed') +
      geom_vline(xintercept = -0.5, linetype='dashed') +
      ggtitle(paste0(pri_var,' DE Genes (n=',n_sig,'; dir:',levels(col_data[[pri_var]])[1],' to ',levels(col_data[[pri_var]])[2],')
Adjusted p-value <= 0.05; |Log2FC| >= 0.5')) +
      ylab(bquote('-'~Log[10]~ 'Adjusted p-value (Bonferroni)')) +
      scale_color_manual(values=c('grey','black')) +
      theme_bw()
    plot(p)
    # Generate an accompanying pvalue histogram
    p <- res_vp %>%
      ggplot(aes(x = padj, fill = SIG)) +
      #geom_freqpoly() +
      ggtitle(paste0(pri_var,' DE Genes (n=',n_sig,'; dir:',levels(col_data[[pri_var]])[1],' to ',levels(col_data[[pri_var]])[2],')
Adjusted p-value <= 0.05; |Log2FC| >= 0.5')) +
      geom_dotplot(aes(x = padj), dotsize = 0.5, alpha = 0.5) +
      geom_vline(xintercept = 0.05, linetype='dashed') +
      xlab('Adjusted p-value (Bonferroni)') +
      ylab('Density') +
      scale_fill_manual(values=c('grey','black')) +
      theme_bw()
    plot(p)
  }
  
    # Pathway Analysis on DE genes
    ################################################################################
    de_df <- res_vp %>% filter(SIG == 'SIG')
    # Convert gene symbols to mouse orthologs
    de_df <- rat_mouse_ortho(as.data.frame(de_df), column = 'ENSEMBL_RAT', direction = 'rat2mouse')
    
    # COllect the gene symbols
    de_df$SYMBOL_MOUSE <- mapIds(org.Mm.eg.db, de_df$ENSEMBL_MOUSE, 'SYMBOL', 'ENSEMBL')
    de_df$SYMBOL_MOUSE <- as.character(de_df$SYMBOL_MOUSE)
    
    # Collect gene symbols
    symbol_up <- de_df %>%
      filter(!is.na(SYMBOL_MOUSE)) %>%
      filter(log2FoldChange > 0) %>%
      select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()
    symbol_down <- de_df %>%
      filter(!is.na(SYMBOL_MOUSE)) %>%
      filter(log2FoldChange < 0) %>%
      select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()
    
    # Collect the desired pathways to query for enrichment
    desired <- c('KEGG_2019_Mouse')
    # db <- c('KEGG_2019_Mouse')
    enriched_up <- enrichr(symbol_up, desired)
    enriched_down <- enrichr(symbol_down, desired)
    n<-1
    enrichr_df <- data.frame()
    for(db in desired){
      # Collapse the data into flat format
      if(nrow(enriched_up[[desired]]) > 0){
        enriched_up[[n]]$Data.Base <- db
        enriched_up[[n]]$UPDOWN <- 'Up'
      }
      if(nrow(enriched_down[[desired]]) > 0){
        enriched_down[[n]]$Data.Base <- db
        enriched_down[[n]]$UPDOWN <- 'Down'
      }
      if(nrow(enriched_down[[desired]]) > 0){
        enrichr_df <- rbind(enrichr_df,enriched_up[[n]])
      }else if(nrow(enriched_up[[desired]]) > 0){
        enrichr_df <- rbind(enrichr_df,enriched_down[[n]])
      }else{
        enrichr_df <- rbind(enrichr_df,enriched_up[[n]],enriched_down[[n]])
      }
      n <- n + 1
    }
    # Generate a new column for the number of genes per enriched term
    if(nrow(enrichr_df) > 0){
    enrichr_df  <- enrichr_df %>% separate(Overlap, sep = '/', into = c('Overlap.N','Total.N'))
    enrichr_df$FDR_nlog <- -log(enrichr_df$Adjusted.P.value)
    # Filter the results
    enrichr_plot <- enrichr_df %>%
      filter(Total.N <= 500) %>%
      filter(Total.N >= 10) %>%
      filter(Adjusted.P.value <= 0.1) %>%
      arrange(desc(FDR_nlog)) %>%
      dplyr::slice(1:40) %>%
      group_by(Term, UPDOWN) %>%
      arrange(FDR_nlog) %>%
      ungroup()
    # Adjust levels
    enrichr_plot$UPDOWN <- factor(enrichr_plot$UPDOWN, levels = c('Up','Down'))
    # Create a special column that will allow for proper ordering of Pathways in Plot:
    # Ordering will have 2 priorities in this order: Up/Down regulation & -Log(FDR)
    enrichr_plot <- enrichr_plot %>%
      mutate(FDR_ord = ifelse(UPDOWN == 'Down', -FDR_nlog, FDR_nlog))
    # Visualize results with barplot
    #pdf(paste0(WD,'/plots/20200426_rnaseq-',TIS,'-pathways-sexmod-C0vsC7_steep.pdf'),
    #    width = 12, height = 6)
    p <- ggplot(enrichr_plot, 
                aes(x=reorder(Term, FDR_ord), y=FDR_nlog, fill = UPDOWN)) +
      geom_bar(stat='identity') +
      coord_flip() +
      ylab('-Log(FDR)') +
      xlab('') +
      theme_linedraw() +
      guides(fill=guide_legend(title='Gene Expression')) +
      scale_fill_manual(values = c('Down' = '#00BFC4','Up' = '#F8766D')) +
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(size=12, colour = 'black'),
            axis.ticks = element_line(colour = 'black'),
            axis.title=element_text(size=12,face='bold'),
            strip.text = element_blank()) +
      theme(panel.grid.major = element_blank()) +
      ggtitle(paste0('Pathway Enrichment (',desired,')
',pri_var,' DE Genes (n=',n_sig,'; dir:',levels(col_data[[pri_var]])[1],' to ',levels(col_data[[pri_var]])[2],')
Adjusted p-value <= 0.05; |Log2FC| >= 0.5'))
    plot(p)
    #dev.off()
    print('No DE Genes')
    }
  
  #' ### Adjust Variance
  
  #+ Adjust Variance
  ################################################################################
  ########### Adjust Variance  #######################################
  ################################################################################
  if(PRI_VAR != 'None'){
    #pri_var <- 'animal.registration.sex'
    #pri_var <- 'calculated.variables.deathtime_after_fed'
    #pri_var <- 'animal_time_last_fed'
    # Set primary variable of interest
    if(TISSUE == c('Gastrocnemius_MSSM_1')){
      tod_cols <- col_data %>%
        filter(Tissue == 'Gastrocnemius') %>%
        filter(Seq_batch == 'MSSM_1') %>%
        filter(!is.na(animal.registration.sex))
    }else if(TISSUE == c('Gastrocnemius_Stanford_1')){
      tod_cols <- col_data %>%
        filter(Tissue == 'Gastrocnemius') %>%
        filter(Seq_batch == 'Stanford_1') %>%
        filter(!is.na(animal.registration.sex))
    }else if(TISSUE == c('Hippocampus')){
      tod_cols <- col_data %>%
        filter(Tissue == 'Hippocampus') %>%
        filter(sample_key %!in% OUTLIERS) %>%
        filter(!is.na(animal.registration.sex))
      }else{
      # Filter Samples (meta)
      tod_cols <- col_data %>%
        filter(Tissue == TISSUE) %>%
        filter(!is.na(animal.registration.sex))
    }
    rownames(tod_cols) <- tod_cols$sample_key
    
    # Collect samples without NA values in TOD
    nona_sams <- tod_cols %>%
      filter(!is.na(specimen.collection.t_death_hour)) %>%
      filter(!is.na(animal.registration.sex)) %>%
      select(sample_key) %>% unlist() %>% as.character()
    # Collect tissue specific counts
    tod_counts <- count_data[,nona_sams]
    #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
    all(rownames(tod_cols) == colnames(tod_counts))
    
    # Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
    design = formula(paste0('~ ',pri_var)) # Primary variable needs to be last.
    dds_final <- DESeqDataSetFromMatrix(countData = tod_counts,
                                        colData = tod_cols,
                                        design = design)
    dds_final <- dds_final[rowSums(counts(dds_final))/ncol(dds_final) >= 1,]
    dds_final <- estimateSizeFactors(dds_final)
    rld_final <- DESeq2::vst(dds_final, blind = F)
    # The batch effect can only be removed with limma
    # https://support.bioconductor.org/p/76099/ (See Michael Love's Comment)
    assay(rld_final) <- limma::removeBatchEffect(assay(rld), rld[[pri_var]])
    
    # Examine the primary variable of interest to see if we've solved our issue
    # Before:
    p <- DESeq2::plotPCA(rld, intgroup = pri_var, ntop = 20000) +
      guides(color=guide_legend(title='')) + coord_equal() +
      ggtitle('PCA of Kidney')
    plot(p)
    # After
    p <- DESeq2::plotPCA(rld_final, intgroup = pri_var, ntop = 20000) +
      guides(color=guide_legend(title='')) + coord_equal() +
      ggtitle('PCA of Kidney')
    plot(p)
    rld <- rld_final
  }
}

################################################################################
###### Manual Annotation Opportunity: Take Note of FINAL_FORMULA ###############
################################################################################

################################################################################
# Move back to beginning of for loop and investigate additional batch factors ##
################################################################################

# Find variables associated with PCs
pc_cor_df <- cor_PC_1_4(rld = rld, ntop = 20000, intgroups = names(colData(rld))) %>%
  filter(Adjusted_R_Sq > 0.15) %>%
  arrange(PC)
# Examine Variables of interest for Pcs 1 and 2
pc_cor_df %>%
  filter(PC %in% c(1:4))

# Visualize variables in PC's 1 and 2
voi <- pc_cor_df %>%
  filter(PC %in% c(1,2)) %>%
  select(Condition) %>% unlist() %>% as.character %>% unique()
# Plot the top variables associated with PCs 1 or 2
for( v in voi){
  # Quick investigation of variables
  #v <- 'RIN'
  title_parts <- pc_cor_df %>%
    filter(Condition == v)
  title <- paste0(title_parts$Condition,'\nPC = ',
                  title_parts$PC,'; R2 = ',
                  round(title_parts$Adjusted_R_Sq,3))
  p <- DESeq2::plotPCA(rld, intgroup =v, ntop = 20000) +
    ggtitle(title)
  plot(p)
}
# To examine PC's 3 and 4
if(T){
  # Visualize variables in PC's 3 - 4
  voi34 <- pc_cor_df %>%
    filter(PC %in% c(3:4)) %>%
    select(Condition) %>% unlist() %>% as.character %>% unique()
  
  # Prepare the data for PC3 and PC4 visualization
  pca <- prcomp(t(assay(rld)), scale. = F)
  df_plot <- pca$x %>% as.data.frame()
  df_plot$sample_key <- row.names(df_plot)
  df_plot <- df_plot %>% 
    left_join(y = col_data, by = 'sample_key')
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  # Plot the top variables associated with PCs 1 or 2
  for( v in voi34){
    # Quick investigation of variables
    title_parts <- pc_cor_df %>%
      filter(Condition == v)
    title <- paste0(title_parts$Condition,'\nPC = ',
                    title_parts$PC,'; R2 = ',
                    round(title_parts$Adjusted_R_Sq,3))
    p <- df_plot %>%
      ggplot(aes(x = PC3, y = PC4, color = !!as.symbol(v))) +
      geom_point(size = 2) +
      coord_equal() +
      xlab(paste0('PC3: ',round(percentVar[3]*100,2),'% variance')) +
      ylab(paste0('PC4: ',round(percentVar[4]*100,2),'% variance')) + 
      ggtitle(title)
    plot(p)
  }
}
# Plot PCA with corrected counts
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c(pri_var, 'sample_key','animal.key.anirandgroup',
                                      'animal.registration.sex','animal.registration.cagenumber','calculated.variables.deathtime_after_fed'), 
                           returnData=TRUE, ntop = 20000)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
# Examine the PCA annotated with variables of interest
pcaData %>%
  mutate(sample_key = as.character(sample_key)) %>%
  mutate(sample_misid_suspect = ifelse(sample_key %in% MIS_SUS, sample_key, '')) %>%
  ggplot(aes(PC1, PC2, color=animal.key.anirandgroup)) +
  geom_point(size=4) +
  #geom_label_repel(aes(label=!!as.symbol(pri_var)),hjust=0, vjust=0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  coord_fixed() +
  ggtitle(paste0(TISSUE,' Gene Expression (Post Batch Correction)')) +
  guides(color=guide_legend(title='animal.key.anirandgroup')) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())

#' ## Outlier Detection (By Sample)

#+ Outlier Detection (By Sample)
################################################################################
#####     Outlier Detection (By Sample)     ####################################
################################################################################

# outlier strategy is adapted from Florian Privé (https://www.r-bloggers.com/detecting-outlier-samples-in-pca/)
pca <- prcomp(t(assay(rld)), scale. = F)
PC <- pca$x %>% as.data.frame()
PC$sample_key <- row.names(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
#' #### Variance-explained plot: explained variance for PCs
# This iterates through PCs and identifies 
df_plot <- data.frame(PC = seq_along(pca$sdev^2/sum(pca$sdev^2)),
                      VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
  mutate(ELBOW = ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE'))
p <- df_plot %>%  
  ggplot(aes(PC,VAR, color = ELBOW)) +
  geom_point() +
  geom_vline(xintercept = elbow_finder(df_plot$PC,df_plot$VAR)[1], linetype='dashed')
plot(p)
PC_meta <- t(pca$x) %>%
  as_tibble() %>%
  mutate(PC = seq_along(pca$sdev^2/sum(pca$sdev^2))) %>%
  mutate(VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
  mutate(ELBOW = ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE')) %>%
  filter(ELBOW == 'UPPER') %>%
  mutate(PC = as.character(sprintf('PC%d',PC))) %>%
  filter(PC %in% c('PC1','PC2','PC3','PC4'))
# Reset the pca object
pca <- prcomp(t(assay(rld)), scale. = F, 
              rank. = elbow_finder(df_plot$PC,df_plot$VAR)[1] - 1)

# SUbset the PC dataframe to only include PCs above the elbow of PC importance
# We can also detect the distance from the median in units of standard deviation (SD_DIST); use median() instead of mean() and mad() instead of sd() because they are more robust estimators in outlier context.
# Tukey considered data > 6 standard deviations as outliers
# #Instead of using the infinite distance, Mahalanobis distance is a multivariate distance based on all variables (PCs here) at once. We use a robust version of this distance, which is implemented in packages {robust} and {robustbase} (Gnanadesikan and Kettenring 1972, Yohai and Zamar (1988), Maronna and Zamar (2002), Todorov, Filzmoser, and others (2009)) and that is reexported in {bigutilsr}.
#This new criterion provides similar results for this data. These robust Mahalanobis distances are approximately Chi-square distributed, which enables deriving p-values of outlierness.
#Local Outlier Factor (LOF)
#LOF statistic (Breunig et al. 2000) has been cited more than 4000 times. Instead of computing a distance from the center, it uses some local density of points. We make use of the fast K nearest neighbours implementation of R package {nabor} (Elseberg et al. 2012) to implement this statistic efficiently in {bigutilsr}.
#To solve the problem of skewness, the medcouple (mc) has been introduced (Hubert and Vandervieren 2008) and is implemented in robustbase::adjboxStats().
#Tukey’s rule uses a fixed coefficient (1.5) that does not account for multiple testing, which means that for large samples, you will almost always get some outliers if using 1.5.
#To solve these two issues, we implemented tukey_mc_up() that accounts both for skewness and multiple testing by default.

# Collect Standard Deviations
pca$sd <- apply(pca$x, 2, function(x) abs(x - median(x)) / mad(x))
join_meta <- PC_meta[,c('PC','VAR')]
names(join_meta) <- c('SD_DIST_PC','PC_VAR')

PC <- PC %>%
  select(
    c(sprintf('PC%d',1:(elbow_finder(df_plot$PC,df_plot$VAR)[1]-1)),'sample_key')) %>%
  mutate(sample_n = seq_along(PC1)) %>%
  mutate(SD_DIST = apply(pca$x, 2, function(x) abs(x - median(x)) / mad(x)) %>%
           apply(1, max)) %>%
  mutate(SD_DIST_PC = colnames(pca$sd)[apply(pca$sd,1,which.max)]) %>%
  mutate(TUKEY_OUTLIER = ifelse(SD_DIST > 6, 'Outlier', 'Normal')) %>%
  mutate(MAHALANOBIS_DIST = dist_ogk(pca$x)) %>%
  # Bonferroni correction
  mutate(MAHALANOBIS_FDR = pchisq(MAHALANOBIS_DIST, df = elbow_finder(df_plot$PC,df_plot$VAR)[1], lower.tail = FALSE)*length(MAHALANOBIS_DIST)) %>%
  mutate(MAHALANOBIS_OUTLIER = ifelse(MAHALANOBIS_FDR <= 0.05, 
                                      'Outlier', 'Normal')) %>%
  mutate(LOF = LOF(pca$x, log = F)) %>%
  left_join(y = col_data, by = 'sample_key') %>%
  left_join(y = join_meta, by = 'SD_DIST_PC') %>%
  arrange(desc(SD_DIST))

# Plot 2 distances (MAHALANOBIS & Local Outlier Factor) on opposing axises and apply Tukey's rule for outlier detection
p <- PC %>%
  mutate(PC_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                             log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), SD_DIST_PC, '')) %>%
  ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = animal.key.anirandgroup)) +
  geom_point(size = 2) +
  scale_color_manual(values=ec_colors) +
  geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
  geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red')
plot(p)
# Examine the variance from of the PC that is associated with the outlier
p <- PC %>%
  mutate(PC_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                             log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), SD_DIST_PC, '')) %>%
  ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = PC_VAR)) +
  geom_point(size = 2) +
  geom_label_repel(aes(label=PC_LABEL)) +
  #geom_label_repel(aes(label=sample_key)) +
  geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
  geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red')
plot(p)

# Examine a PCA plot of points
p <- PC %>%
  mutate(OUTLIER = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                             log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, '')) %>%
  ggplot(aes(x = PC1, y = PC2, color = SD_DIST)) +
  geom_point(size = 3) +
  coord_equal() +
  geom_label_repel(aes(label=OUTLIER)) +
  xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
  ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) + 
  ggtitle(paste0('PCA of ',TISSUE,' Gene Expression:\n',FINAL_FORMULA))
plot(p)

hist(PC$SD_DIST, breaks = 78)
hist(PC$SD_DIST, breaks = 78)

################################################################################
###### Manual Annotation Opportunity: Fill in OUTLIER variable #################
################################################################################
# Note: If Samples Labled as Outliers are in Acute exercise groups, then they likely do not meet our criteria for outliers.

# In case you need to make last minute changes
if(TISSUE %in% c('Brown Adipose')){
  # Filter Samples (meta)
  tod_cols <- col_data %>%
    filter(Tissue == TISSUE) %>%
    filter(sample_key %!in% OUTLIERS) %>%
    filter(!is.na(animal.registration.sex))
  rownames(tod_cols) <- tod_cols$sample_key
  # Collect samples without NA values in TOD
nona_sams <- tod_cols %>%
  filter(!is.na(specimen.collection.t_death_hour)) %>%
  filter(!is.na(animal.registration.sex)) %>%
  select(sample_key) %>% unlist() %>% as.character()
# Collect tissue specific counts
tod_counts <- count_data[,nona_sams]
#### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))
# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
design = formula(paste0('~ ',pri_var)) # Primary variable needs to be last.
dds <- DESeqDataSetFromMatrix(countData = tod_counts,
                                    colData = tod_cols,
                                    design = design)
dds <- dds[rowSums(counts(dds))/ncol(dds) >= 1,]
dds <- estimateSizeFactors(dds)
rld <- DESeq2::vst(dds, blind = F)
# The batch effect can only be removed with limma
# https://support.bioconductor.org/p/76099/ (See Michael Love's Comment)
assay(rld) <- limma::removeBatchEffect(assay(rld), rld[[pri_var]])
}

# Plot PCA with corrected counts
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c(pri_var, 'sample_key','animal.key.anirandgroup',
                                      'animal.registration.sex'), 
                           returnData=TRUE, ntop = 20000)
percentVar <- round(100 * attr(pcaData, 'percentVar'))
# Examine the PCA annotated with variables of interest
pcaData %>%
  mutate(sample_key = as.character(sample_key)) %>%
  mutate(sample_misid_suspect = ifelse(sample_key %in% MIS_SUS, sample_key, '')) %>%
  ggplot(aes(PC1, PC2, color=animal.key.anirandgroup)) +
  geom_point(size=3) +
  #geom_label_repel(aes(label=!!as.symbol(pri_var)),hjust=0, vjust=0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  coord_fixed() +
  ggtitle(paste0(TISSUE,' Gene Expression (Post Batch Correction)')) +
  guides(color=guide_legend(title='animal.key.anirandgroup')) +
  scale_color_manual(values=ec_colors) +
  theme(legend.title=element_blank())

################################################################################
################ Additional Code Notes #########################################
################################################################################


# # TODO: Generate a decision tree for which Var-Specific Strategy to utilize
# # Strategy Symbol Key:
# # Remove - Rm
# # Median center - Mc
# # Model - Mo
# # Sex Genes - Sg
# # Autosomal Genes - Au
# # DE Genes (Sex) - DEg
# # All Genes - Ag
# 
# # Collect Samples from dichotomous groups
# L1_samples <- col_data %>%
#   filter(Tissue == TISSUE) %>%
#   filter(!!as.symbol(pri_var) == levels(col_data[[pri_var]])[1]) %>%
#   dplyr::select(sample_key) %>% unlist() %>% as.character()
# L2_samples <- col_data %>%
#   filter(Tissue == TISSUE) %>%
#   filter(!!as.symbol(pri_var) == levels(col_data[[pri_var]])[2]) %>%
#   dplyr::select(sample_key) %>% unlist() %>% as.character()
# SIG_genes <- res_vp %>%
#   filter(SIG == 'SIG') %>%
#   dplyr::select(ENSEMBL_RAT) %>% unlist() %>% as.character()
# NS_genes <- res_vp %>%
#   filter(SIG == 'NS') %>%
#   dplyr::select(ENSEMBL_RAT) %>% unlist() %>% as.character()
# # Select the counts
# L1SIG_counts <- assay(rld_final)[SIG_genes, L1_samples]
# L2SIG_counts <- assay(rld_final[SIG_genes, L2_samples])
# L1NS_counts <- assay(rld_final[NS_genes, L1_samples])
# L2NS_counts <- assay(rld_final[NS_genes, L2_samples])
# 
# 
# if(VAR_STGY == 'MC_DEg'){
#   # Median Center Significant Genes
#   # Collects median of each row, then subtracts by row medians
#   L1_medians <- apply(L1SIG_counts,1,median)
#   L1SIG_counts <- L1SIG_counts - L1_medians
#   L2_medians <- apply(L2SIG_counts,1,median)
#   L2SIG_counts <- L2SIG_counts - L2_medians
#   SIG_counts <- cbind(L1SIG_counts, L2SIG_counts)
#   
#   # Non-Significant Genes
#   NS_counts <- cbind(L1NS_counts, L2NS_counts)
#   MC_DEg_counts <- rbind(SIG_counts, NS_counts)
#   MC_DEg_counts <- MC_DEg_counts[, colnames(assay(rld_final))]
#   assay(rld_final) <- MC_DEg_counts
# }else if(VAR_STGY == 'MC_Ag'){
#   # Median Center Significant Genes
#   # Collects median of each row, then subtracts by row medians
#   M_medians <- apply(MSIG_counts,1,median)
#   MSIG_counts <- MSIG_counts - M_medians
#   F_medians <- apply(FSIG_counts,1,median)
#   FSIG_counts <- FSIG_counts - F_medians
#   SIG_counts <- cbind(MSIG_counts, FSIG_counts)
#   
#   # Non-Significant Genes
#   M_medians <- apply(MNS_counts,1,median)
#   MNS_counts <- MNS_counts - M_medians
#   F_medians <- apply(FNS_counts,1,median)
#   FNS_counts <- FNS_counts - F_medians
#   NS_counts <- cbind(MNS_counts, FNS_counts)
#   MC_Ag_counts <- rbind(SIG_counts, NS_counts)
#   MC_Ag_counts <- MC_Ag_counts[, colnames(assay(rld_final))]
#   #rld <- rld
#   assay(rld_final) <- MC_Ag_counts
# }


# # Fit the LDA
# lda.fit = train(design, data=train_df, method='lda',
#                 trControl = control)
# 
# lda.fit
# summary(lda.fit)
# 
# pred_sex = predict(lda.fit, test_df)
# table(pred_sex, test_df[[mis_id_var]])
# 
# pred.accuracy = round(mean(pred_sex == test_df[[mis_id_var]])*100,2)
# pred.accuracy



# Examine counts of top gene
# assay(rld)['ENSRNOG00000050000',MIS_SUS] %>% as.numeric()
# t <- assay(rld)['ENSRNOG00000050000',] %>% as.data.frame()
# t$sample <- row.names(t)
# t %>%
#   filter(sample %!in% MIS_SUS) %>% as.numeric()
# 
# t <- counts(dds, normalized = TRUE) %>% as.data.frame()
# t$ENSEMBL_RAT <- row.names(t)
# t <- t %>% 
#   filter(ENSEMBL_RAT == 'ENSRNOG00000050000') %>%
#   t() %>% as.data.frame()
# t$SAMPLE <- row.names(t)
# t <- t[-79,]
# b1 <- (col_data %>% filter(animal.key.batch == '1'))$sample_key %>% as.character()
# t <- t %>% 
#   mutate(COHORT = ifelse(SAMPLE %in% MIS_SUS, 'SUS', 'NORM')) %>%
#   mutate(BATCH = ifelse(SAMPLE %in% b1, '1', '2'))
# t <- t %>% as_tibble()
# t <- t %>% mutate(COUNT = ENSRNOG00000050000 %>% unlist() %>% as.character() %>% as.numeric())
# 
# t %>%
#   ggplot(aes(x = BATCH, y = COUNT, color = COHORT)) +
#   geom_point()
#plotCounts(dds, gene='ENSRNOG00000050000', intgroup='animal.key.batch')



