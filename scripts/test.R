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
#BiocManager::install('ImpulseDE2')
#install.packages('tidyverse')

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr','ImpulseDE2')
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
for(TISSUE in c('Hippocampus')){
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





#


# tod_cols <- col_data %>%
#         filter(Tissue == TISSUE) %>%
#         filter(sample_key %!in% OUTLIERS) %>%
#         filter(!is.na(animal.registration.sex))
# rownames(tod_cols) <- tod_cols$sample_key
# nona_sams <- tod_cols %>%
#         filter(!is.na(specimen.collection.t_death_hour)) %>%
#         filter(!is.na(animal.registration.sex)) %>%
#         select(sample_key) %>% unlist() %>% as.character()
# tod_counts <- count_data[,nona_sams]
# all(rownames(tod_cols) == colnames(tod_counts))
# design = ~ 1 # Primary variable needs to be last.
# title = paste0('Design: ',as.character(design))
# dds1 <- DESeqDataSetFromMatrix(countData = tod_counts,
#                                colData = tod_cols,
#                                design = design)
# zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% 
#         nrow() %>% as.character()
# reads_n <- 1
# keep <- rowSums(counts(dds1))/ncol(dds1) >= reads_n
# dds2 <- dds1[keep,]
# dds2
# dds <- dds2
# sort(colSums(assay(dds)))/1e6
# dds <- estimateSizeFactors(dds)
# summary(sizeFactors(dds))
# rld <- DESeq2::vst(dds, blind = FALSE)
# rs <- rowSums(counts(dds))
