#'---
#' title: "PASS1B Rat RNASeq: -- Data Assembly"
#' author: "Alec Steep"
#' date: "20200718"
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
#' * Working directory
#' * Load R packages

#+ Setup Environment, message=FALSE, results='hide', warning = FALSE
################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source('https://bioconductor.org/biocLite.R')
#BiocManager::install('MotrpacBicQC')
#install.packages('prettydoc')
#devtools::install_github("MoTrPAC/MotrpacBicQC", build_vignettes = TRUE)

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr','ImpulseDE2','data.table','prettydoc','MotrpacBicQC','zoo')
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
source(paste0(WD,'/functions/simpleCap.R'))

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
#' #### Section to generate a table of major decisions

# Obtain tissue codes
bic_animal_tissue_code = data.table(MotrpacBicQC::bic_animal_tissue_code)
bic_animal_tissue_code <- bic_animal_tissue_code %>%
        mutate(tissue_key = paste(bic_tissue_code, bic_tissue_name)) %>%
        mutate(tissue_key = str_remove_all(tissue_key, ' Powder')) %>%
        mutate(tissue_key = str_replace_all(tissue_key, ' ', '-')) %>%
        mutate(tissue_key = tolower(tissue_key))
bic_animal_tissue_code$tissue_key

# Create the 
df_tbl <- data.frame(matrix(ncol = 8, nrow = 0))
names(df_tbl) <- c('Tissue','Tis','Seq_Batch','Misidentified',
                   'Adjusted_Variables','Outliers','Formula')

for(tissue in c('t65-aorta','t52-hippocampus',
                't54-hypothalamus','t55-gastrocnemius','t56-vastus-lateralis')){
        # for(TISSUE in c('Hypothalamus', 'Kidney', 'Aorta', 'Adrenal', 'Brown Adipose', 'Cortex', 'Gastrocnemius', 'Heart', 'Hippocampus','Lung','Ovaries','Spleen', 'White Adipose','Liver','Testes')){
        # Declare Outliers
        if(tissue == 'Kidney'){
                TIS <- 'KID'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 'Liver'){
                TIS <- 'LIV'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                OUTLIERS <- c('None')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
        }else if(tissue == 't54-hypothalamus'){
                TISSUE <- 'Hypothalamus'
                TIS <- 'SCN'
                MIS_ID_VAR <- c('')
                MIS_SUS <- c('')
                MIS_ID <- c('')
                PRI_VAR <- c('')
                FINAL_FORMULA <- ('')
                OUTLIERS <- c('')
        }else if(tissue == 't56-vastus-lateralis'){
                TISSUE <- 'Vastus Lateralis'
                TIS <- 'VLS'
                MIS_ID_VAR <- c('')
                MIS_SUS <- c('')
                MIS_ID <- c('')
                PRI_VAR <- c('')
                FINAL_FORMULA <- ('')
                OUTLIERS <- c('')
        }else if(tissue == 't65-aorta'){
                TISSUE <- 'Aorta'
                TIS <- 'AOR'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90159016502_SN2')
        }else if(tissue == 't55-gastrocnemius'){
                TISSUE <- 'Gastrocnemius'
                TIS <- 'SKM'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 'Heart'){
                TIS <- 'HAT'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('90052015802_SN1')
        }else if(tissue == 'Adrenal'){
                TIS <- 'ADG'
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 'Brown Adipose'){
                TIS <- 'BAT'
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('90121016905_SF1','90014016905_SF1','90128016905_SF1','90047016905_SF1','90001016905_SF1',
                             '90018016905_SF1')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('90023016905_SF1','90047016905_SF1','90128016905_SF1')
        }else if(tissue == 'White Adipose'){
                TIS <- 'WAT'
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 'Cortex'){
                TIS <- c('COR')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('animal.registration.sex')
                FINAL_FORMULA <- formula(paste0('~ ',PRI_VAR[1],
                                                ' + animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 't52-hippocampus'){
                TISSUE <- 'Hippocampus'
                TIS <- c('HIP')
                MIS_ID_VAR <- c('animal.registration.sex')
                MIS_SUS <- c('90115015202_SF2','90041015202_SF2')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90040015202_SF2','90115015202_SF2','90041015202_SF2',
                              '90009015202_SF2','90005015202_SF2','90139015202_SF2')
        }else if(tissue == 'Lung'){
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
        }else if(tissue == 'Ovaries'){
                TIS <- c('OVR')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('None')
        }else if(tissue == 'Testes'){
                TIS <- c('TES')
                MIS_ID_VAR <- c('None')
                MIS_SUS <- c('None')
                MIS_ID <- c('None')
                PRI_VAR <- c('None')
                FINAL_FORMULA <- formula(paste0('~ animal.key.anirandgroup'))
                OUTLIERS <- c('90017016302_SN2','90031016302_SN2','90127016302_SN2','90015016302_SN2','90129016302_SN2')
        }else if(tissue == 't62-spleen'){
                TISSUE <- 'Spleen'
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
table_file <- paste0(WD,'/data/20200718_pass1b-rnaseq-tissue-data-assembly-table_steep.txt')
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

# Load the decision table
table_file <- paste0(WD,'/data/20200603_rnaseq-tissue-data-assambly-table_steep.txt')
df_tbl <- read.table(file = table_file,sep = '\t', header = T, check.names = F)

models_df <- data.frame()
#TISSUE <- "Kidney"tissue
for(tissue in c('t65-aorta','t52-hippocampus',
                't54-hypothalamus','t55-gastrocnemius','t56-vastus-lateralis')){
        TISSUE <- str_sub(tissue, start = 5) %>% 
                str_replace_all(pattern = '-', replacement = ' ') %>%
                sapply(simpleCap) %>% as.character()
        print(TISSUE)
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
        # FORMULA <- df_tbl %>%
        #         filter(Tissue == TISSUE) %>%
        #         select(Formula) %>% unique() %>% 
        #         unlist() %>% as.character() %>% as.formula()
        
        # Obtain tissue codes
        bic_animal_tissue_code = data.table(MotrpacBicQC::bic_animal_tissue_code)
        bic_animal_tissue_code <- bic_animal_tissue_code %>%
                mutate(tissue_key = paste(bic_tissue_code, bic_tissue_name)) %>%
                mutate(tissue_key = str_remove_all(tissue_key, ' Powder')) %>%
                mutate(tissue_key = str_replace_all(tissue_key, ' ', '-')) %>%
                mutate(tissue_key = tolower(tissue_key))
        bic_animal_tissue_code$tissue_key
        
        # Tissue read counts (raw)
        count_data <- read.table(file = paste0(WD,'/transcriptomics/',tissue,'/transcript-rna-seq/results/motrpac_pass1b-06_',tissue,'_transcript-rna-seq_rsem-genes-count.txt'), 
                                 header = TRUE, sep = '\t', check.names = FALSE)
        
        # Tissue meta data
        col_data <-  fread(file = paste0(WD,'/phenotype/motrpac_pass1b-06_pheno_viallabel-data.txt'), header = TRUE, sep = '\t') %>% as_tibble()
        cols_1 <- fread(file = paste0(WD,'/phenotype/motrpac_pass1b_pheno_merged-dictionary.txt'), header = TRUE, sep = '\t') %>% as_tibble()
        
        
        # Iterate through the column names of meta table and adjust them accordingly
        new_names <- c()
        for(cn in colnames(col_data)){
                if(grepl('_1', cn)){
                        fn <- cols_1 %>%
                                filter(BICUniqueID == str_remove_all(cn, '_1')) %>%
                                select(FullName) %>% as.character()
                        fn <- paste0(fn,'_1')
                }else if(grepl('_2', cn)){
                        fn <- cols_1 %>%
                                filter(BICUniqueID == str_remove_all(cn, '_2')) %>%
                                select(FullName) %>% as.character()
                        fn <- paste0(fn,'_2')        
                }else{
                        fn <- cols_1%>%
                                filter(BICUniqueID == cn) %>%
                                select(FullName) %>% as.character()
                }
                new_names <- c(new_names, fn)
        }
        colnames(col_data) = new_names
        
        # Adjust variable nomenclature (intuitive)
        ############################################
        col_data <- col_data %>%
                mutate(Key.protocol = case_when(
                        Key.protocol == 1 ~ 'Phase 1A',
                        Key.protocol == 2 ~ 'Phase 1B')) %>%
                mutate(Key.agegroup = case_when(
                        Key.agegroup == 1 ~ '6 months',
                        Key.agegroup == 2 ~ '18 months')) %>%
                mutate(Key.intervention = case_when(
                        Key.intervention == 1 ~ 'Exercise',
                        Key.intervention == 2 ~ 'Training',
                        Key.intervention == 3 ~ 'Control')) %>%
                mutate(Key.sacrificetime = case_when(
                        Key.sacrificetime == 1 ~ 'Immediate Post Exercise (Phase 1A)',
                        Key.sacrificetime == 2 ~ '0.5 hour (Phase 1A)',
                        Key.sacrificetime == 3 ~ '1 hour (Phase 1A)',
                        Key.sacrificetime == 4 ~ '4 hour (Phase 1A)',
                        Key.sacrificetime == 5 ~ '7 hour (Phase 1A)',
                        Key.sacrificetime == 6 ~ '24 hour (Phase 1A)',
                        Key.sacrificetime == 7 ~ '48 hour (Phase 1A)',
                        Key.sacrificetime == 8 ~ '1 week of training or control time (Phase 1B)',
                        Key.sacrificetime == 9 ~ '2 weeks of training (Phase 1B)',
                        Key.sacrificetime == 10 ~ '4 weeks of training (Phase 1B)',
                        Key.sacrificetime == 11 ~ '8 weeks of training or control time (Phase 1B)')) %>%
                mutate(Key.sitename = case_when(
                        Key.sitename == 910 ~ 'Joslin',
                        Key.sitename == 930 ~ 'University of Florida',
                        Key.sitename == 940 ~ 'University of Iowa')) %>%
                mutate(Registration.siteid = case_when(
                        Registration.siteid == 910 ~ 'Joslin',
                        Registration.siteid == 930 ~ 'University of Florida',
                        Registration.siteid == 940 ~ 'University of Iowa')) %>%
                mutate(Registration.sex = case_when(
                        Registration.sex == 1 ~ 'Female',
                        Registration.sex == 2 ~ 'Male')) %>%
                mutate(Familiarization.siteid = case_when(
                        Familiarization.siteid == 910 ~ 'Joslin',
                        Familiarization.siteid == 930 ~ 'University of Florida',
                        Familiarization.siteid == 940 ~ 'University of Iowa')) %>%
                mutate(Familiarization.compliant = case_when(
                        Familiarization.compliant == 1 ~ 'Yes',
                        Familiarization.compliant == 0 ~ 'No')) %>%
                mutate(Training.siteid = case_when(
                        Training.siteid == 910 ~ 'Joslin',
                        Training.siteid == 930 ~ 'University of Florida',
                        Training.siteid == 940 ~ 'University of Iowa')) %>%
                mutate(VO2.Max.Test.siteid_1 = case_when(
                        VO2.Max.Test.siteid_1 == 910 ~ 'Joslin',
                        VO2.Max.Test.siteid_1 == 930 ~ 'University of Florida',
                        VO2.Max.Test.siteid_1 == 940 ~ 'University of Iowa')) %>%
                mutate(VO2.Max.Test.siteid_2 = case_when(
                        VO2.Max.Test.siteid_2 == 910 ~ 'Joslin',
                        VO2.Max.Test.siteid_2 == 930 ~ 'University of Florida',
                        VO2.Max.Test.siteid_2 == 940 ~ 'University of Iowa')) %>%
                mutate(NMR.Testing.siteid_1 = case_when(
                        NMR.Testing.siteid_1 == 910 ~ 'Joslin',
                        NMR.Testing.siteid_1 == 930 ~ 'University of Florida',
                        NMR.Testing.siteid_1 == 940 ~ 'University of Iowa')) %>%
                mutate(NMR.Testing.siteid_2 = case_when(
                        NMR.Testing.siteid_2 == 910 ~ 'Joslin',
                        NMR.Testing.siteid_2 == 930 ~ 'University of Florida',
                        NMR.Testing.siteid_2 == 940 ~ 'University of Iowa')) %>%
                mutate(Specimen.Collection.siteid = case_when(
                        Specimen.Collection.siteid == 910 ~ 'Joslin',
                        Specimen.Collection.siteid == 930 ~ 'University of Florida',
                        Specimen.Collection.siteid == 940 ~ 'University of Iowa')) %>%
                mutate(Specimen.Collection.bloodcomplete = case_when(
                        Specimen.Collection.bloodcomplete == 1 ~ 'full',
                        Specimen.Collection.bloodcomplete == 2 ~ 'partial',
                        Specimen.Collection.bloodcomplete == 3 ~ 'N/A',
                        Specimen.Collection.bloodcomplete == 4 ~ 'Not')) %>%
                mutate(Specimen.Collection.uteruscomplete = case_when(
                        Specimen.Collection.uteruscomplete == 1 ~ 'full',
                        Specimen.Collection.uteruscomplete == 2 ~ 'partial',
                        Specimen.Collection.uteruscomplete == 3 ~ 'N/A',
                        Specimen.Collection.uteruscomplete == 4 ~ 'Not')) %>%
                mutate(Specimen.Collection.deathtype = case_when(
                        Specimen.Collection.deathtype == 1 ~ 'Decapitation',
                        Specimen.Collection.deathtype == 2 ~ 'Heart Removal')) %>%
                mutate(Specimen.Processing.siteid = case_when(
                        Specimen.Processing.siteid == 910 ~ 'Joslin',
                        Specimen.Processing.siteid == 930 ~ 'University of Florida',
                        Specimen.Processing.siteid == 940 ~ 'University of Iowa')) %>%
                mutate(Specimen.Processing.hemolyzed = case_when(
                        Specimen.Processing.hemolyzed == 1 ~ 'Yes',
                        Specimen.Processing.hemolyzed == 2 ~ 'No')) %>%
                mutate(sample_key = viallabel) %>%
                mutate(Omic_batch = case_when(
                        BICLabelData.shiptositeid == 810 ~ 'PNNL',
                        BICLabelData.shiptositeid == 820 ~ 'BROAD',
                        BICLabelData.shiptositeid == 821 ~ 'BROAD_DUKE',
                        BICLabelData.shiptositeid == 822 ~ 'BROAD_CARR',
                        BICLabelData.shiptositeid == 823 ~ 'BROAD_CLISH',
                        BICLabelData.shiptositeid == 830 ~ 'STAN',
                        BICLabelData.shiptositeid == 840 ~ 'EMORY',
                        BICLabelData.shiptositeid == 841 ~ 'EMORY_GEORGIATECH',
                        BICLabelData.shiptositeid == 850 ~ 'MSSM',
                        BICLabelData.shiptositeid == 860 ~ 'MAYO',
                        BICLabelData.shiptositeid == 870 ~ 'UMICH')) %>%
                mutate(Key.anirandgroup = case_when(
                        (Key.intervention == 'Control' & Key.sacrificetime == '1 week of training or control time (Phase 1B)') ~ 'Control 1wk',
                        (Key.intervention == 'Training' & Key.sacrificetime == '1 week of training or control time (Phase 1B)') ~ 'Training 1wk',
                        (Key.intervention == 'Training' & Key.sacrificetime == '2 weeks of training (Phase 1B)') ~ 'Training 2wk',
                        (Key.intervention == 'Training' & Key.sacrificetime == '4 weeks of training (Phase 1B)') ~ 'Training 4wk',
                        (Key.intervention == 'Control' & Key.sacrificetime == '8 weeks of training or control time (Phase 1B)') ~ 'Control 8wk',
                        (Key.intervention == 'Training' & Key.sacrificetime == '8 weeks of training or control time (Phase 1B)') ~ 'Training 8wk'))
        
        #Key.protocol: Phase 1B
        #Key.agegroup: 6 months
        #Key.intervention: Control
        #Key.sacrificetime: 8 weeks of training or control time (Phase 1B)
        
        # Adjust column objects
        ########################
        # To factors
        factor_cols <- c('pid','bid','labelid','viallabel','Key.agegroup','Key.protocol',
                         'Key.intervention','Key.sacrificetime','Key.anirandgroup','Key.sitename',
                         'Registration.formname','Registration.siteid','Registration.ratid',
                         'Registration.batchnumber','Registration.cagenumber',
                         'Registration.comments','Familiarization.staffid',
                         'Familiarization.siteid','Registration.sex','Familiarization.formname',
                         'Familiarization.versionnbr','Familiarization.compliant',
                         'Familiarization.comments','Training.staffid','Training.siteid',
                         'Training.formname','Training.versionnbr','Training.day1_comments',
                         'Training.day2_comments','Training.day3_comments','Training.day4_comments',
                         'Training.day5_comments','Training.day6_comments','Training.day7_comments',
                         'Training.day8_comments','Training.day9_comments',
                         'Training.day10_comments','Training.day11_comments',
                         'Training.day12_comments','Training.day13_comments',
                         'Training.day14_comments','Training.day15_comments',
                         'Training.day16_comments','Training.day17_comments',
                         'Training.day18_comments','Training.day19_comments',
                         'Training.day20_comments','Training.day21_comments',
                         'Training.day22_comments','Training.day23_comments',
                         'Training.day24_comments','Training.day25_comments',
                         'Training.day26_comments','Training.day27_comments',
                         'Training.day28_comments','Training.day29_comments',
                         'Training.day30_comments','Training.day31_comments',
                         'Training.day32_comments','Training.day33_comments',
                         'Training.day34_comments','Training.day35_comments',
                         'Training.day36_comments','Training.day37_comments',
                         'Training.day38_comments','Training.day39_comments',
                         'Training.day40_comments','VO2.Max.Test.staffid_1',
                         'VO2.Max.Test.staffid_2','VO2.Max.Test.siteid_1',
                         'VO2.Max.Test.siteid_2','VO2.Max.Test.visitcode_1',
                         'VO2.Max.Test.visitcode_2','VO2.Max.Test.formname_1',
                         'VO2.Max.Test.formname_2','VO2.Max.Test.versionnbr_1',
                         'VO2.Max.Test.versionnbr_2','VO2.Max.Test.vo2_comments_1',
                         'VO2.Max.Test.vo2_comments_2','NMR.Testing.staffid_1',
                         'NMR.Testing.staffid_2','NMR.Testing.visitcode_1',
                         'NMR.Testing.visitcode_2','NMR.Testing.formname_1',
                         'NMR.Testing.formname_2','NMR.Testing.versionnbr_1',
                         'NMR.Testing.versionnbr_2','NMR.Testing.nmr_comments_1',
                         'NMR.Testing.nmr_comments_2','Specimen.Collection.staffid',
                         'Specimen.Collection.visitcode','Specimen.Collection.formname',
                         'Specimen.Collection.versionnbr','Specimen.Collection.anesthesiaid',
                         'Specimen.Collection.anesthesiacomments','Specimen.Collection.bloodtype',
                         'Specimen.Collection.bloodtube','Specimen.Collection.bloodtechid',
                         'Specimen.Collection.bloodcomplete','Specimen.Collection.bloodcomments',
                         'Specimen.Collection.uterustype','Specimen.Collection.uterustechid',
                         'Specimen.Collection.uteruscomplete','Specimen.Collection.uteruscomments',
                         'Specimen.Processing.versionnbr','Specimen.Processing.formname',
                         'Specimen.Processing.visitcode',
                         'Specimen.Processing.sampletypedescription',
                         'Specimen.Processing.aliquotdescription','Specimen.Processing.volume',
                         'Specimen.Processing.comments','Specimen.Processing.techid',
                         'BICLabelData.shiptositeid','sample_key','Omic_batch','Key.anirandgroup')
        for(fc in factor_cols){
                col_data[[fc]] <- as.factor(col_data[[fc]])
        }
        
        # Custom factor levels
        ec_levels <- c('Training 1wk','Training 2wk','Training 4wk','Training 8wk','Control 8wk')
        col_data$Key.anirandgroup <- factor(col_data$Key.anirandgroup, levels = ec_levels)
        
        # To Numeric
        num_cols <- c('Familiarization.fat','Familiarization.lean',
                      'Specimen.Collection.uterusweight','Specimen.Processing.partialamt')
        for(nc in num_cols){
                col_data[[nc]] <- as.factor(col_data[[nc]])
        }
        
        # From Characters: 03JUL2018 to Dates
        date_cols <- c('Key.d_arrive','Key.d_sacrifice','Registration.d_visit',
                       'Registration.d_arrive','Registration.d_reverselight',
                       'Familiarization.d_visit','Familiarization.d_treadmillbegin',
                       'Familiarization.d_treadmillcomplete','Training.d_visit',
                       'Training.day1date','Training.day2date','Training.day3date',
                       'Training.day4date','Training.day5date','Training.day6date',
                       'Training.day7date','Training.day8date','Training.day9date',
                       'Training.day10date','Training.day11date','Training.day12date',
                       'Training.day13date','Training.day14date','Training.day15date',
                       'Training.day16date','Training.day17date','Training.day18date',
                       'Training.day19date','Training.day20date','Training.day21date',
                       'Training.day22date','Training.day23date','Training.day24date',
                       'Training.day25date','Training.day26date','Training.day27date',
                       'Training.day28date','Training.day29date','Training.day30date',
                       'Training.day31date','Training.day32date','Training.day33date',
                       'Training.day34date','Training.day35date','Training.day36date',
                       'Training.day37date','Training.day38date','Training.day39date',
                       'Training.day40date','VO2.Max.Test.d_visit_1','VO2.Max.Test.d_visit_2',
                       'VO2.Max.Test.d_vo2_1','VO2.Max.Test.d_vo2_2','NMR.Testing.d_visit_1',
                       'NMR.Testing.d_visit_2','NMR.Testing.d_nmr_1','NMR.Testing.d_nmr_2',
                       'Specimen.Collection.d_visit')
        for(dc in date_cols){
                col_data[[dc]] <- dmy(col_data[[dc]])
        }
        
        # From Dates: Jan 2018
        date_cols <- c('Registration.d_birth')
        for(dc in date_cols){
                col_data[[dc]] <- as.yearmon(col_data[[dc]], "%b%Y")
        }
        
        # From Strings: "10:14" to Times: 10:14:00, then to seconds (better for PCA)
        time_cols <- c('Training.day1time','Training.day2time','Training.day3time',
                       'Training.day4time','Training.day5time','Training.day6time',
                       'Training.day7time','Training.day8time','Training.day9time',
                       'Training.day10time','Training.day11time','Training.day12time',
                       'Training.day13time','Training.day14time','Training.day15time',
                       'Training.day16time','Training.day17time','Training.day18time',
                       'Training.day19time','Training.day20time','Training.day21time',
                       'Training.day22time','Training.day23time','Training.day24time',
                       'Training.day25time','Training.day26time','Training.day27time',
                       'Training.day28time','Training.day29time','Training.day30time',
                       'Training.day31time','Training.day32time','Training.day33time',
                       'Training.day34time','Training.day35time','Training.day36time',
                       'Training.day37time','Training.day38time','Training.day39time',
                       'Training.day40time','VO2.Max.Test.t_vo2_1','VO2.Max.Test.t_vo2_2',
                       'VO2.Max.Test.t_complete_1','VO2.Max.Test.t_complete_2',
                       'Specimen.Collection.t_anesthesia','Specimen.Collection.t_bloodstart',
                       'Specimen.Collection.t_bloodstop','Specimen.Collection.t_edtafill',
                       'Specimen.Collection.t_uterusstart','Specimen.Collection.t_uterusstop',
                       'Specimen.Collection.t_death','Specimen.Processing.t_collection',
                       'Specimen.Processing.t_edtaspin','Specimen.Processing.t_freeze')
        #tc <- 'VO2.Max.Test.t_complete_1'
        for(tc in time_cols){
                col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% 
                        parse_time() %>% as.numeric()
        }
        
        # Generate Additional Variables
        ##################################
        
        # Generate a non-zero variable
        # The fraction of non-zero counts (function)
        nonzero <- function(x) sum(x != 0)
        # Non zero counts
        count_data2 <- count_data %>% select(-gene_id)
        nzc <- apply(count_data2,2,nonzero)
        ginic <- apply(count_data2,2,reldist::gini)
        # Non zero fraction 
        nzf <- nzc/nrow(count_data2)
        # Mean gene expresison
        m <- apply(count_data2,2,mean)
        # SD gene expression
        sd <- apply(count_data2,2,sd)
        # Coefficient of variation
        cv <- sd/m
        # Join the dataframe into meta data
        nz_tis_df <- data.frame('viallabel' = colnames(count_data2),
                                'NZF' = nzf,
                                'GINI' = ginic,
                                'CV' = cv)
        col_data <- left_join(col_data, nz_tis_df, by = c('viallabel'))
        # Save data as an R objects
        # ################################################################################
        # To determine object size
        sl <- object.size(count_data)
        print(sl, units = 'auto')
        # Meta Data
        meta_file <- paste0(WD,'/data/20200718_pass1b-rnaseq-meta-proc_steep.rds')
        saveRDS(col_data, file = meta_file)
        
        # Create objects for Key.anirandgroup colors and levels
        ec_levels <- c('Training 1wk','Training 2wk','Training 4wk','Training 8wk','Control 8wk')
        ec_colors <- c('gold',
                       'darkgoldenrod1',
                       'darkorange2',
                       'darkorange4',
                       'steelblue4')
        
        meta_file <- paste0(WD,'/data/20200718_pass1b-rnaseq-meta-proc_steep.rds')
        #' #### Polished metadata saved as:  
        #' `r meta_file`  
        #'  
        count_file <- paste0(WD,'/transcriptomics/',tissue,'/transcript-rna-seq/results/motrpac_pass1b-06_',tissue,'_transcript-rna-seq_rsem-genes-count.txt')
        #' #### Polished read counts saved as:  
        #' `r count_file`  
        
        # Load Metadata and count data as R objects
        ################################################################################
        # Restore the metadata object
        col_data <- readRDS(file = meta_file)
        # Tissue read counts (raw)
        count_data <- read.table(file = count_file, 
                                 header = TRUE, sep = '\t', check.names = FALSE)
        
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
        # bic_animal_tissue_code
        # SUbset the tissue information
        tod_cols <- col_data %>%
                filter(Specimen.Processing.sampletypedescription == TISSUE) %>%
                as.data.frame()
        rownames(tod_cols) <- tod_cols$viallabel
        
        # Collect samples without NA values in TOD
        nona_sams <- tod_cols %>%
                filter(!is.na(Registration.sex)) %>%
                select(viallabel) %>% unlist() %>% as.character()
        nona_sams <- nona_sams[nona_sams %in% colnames(count_data)]
        tod_cols <- tod_cols[nona_sams,]
        # Collect tissue specific counts
        row.names(count_data) <- count_data$gene_id
        tod_counts <- count_data[,nona_sams]
        tod_counts <- as.data.frame(sapply(tod_counts, as.integer))
        
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
        #' Variance Stabilizing Transform (vst)
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
        pca <- prcomp(t(assay(rld)), scale. = F)
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        PC <- pca$x %>% as.data.frame()
        PC$sample_key <- row.names(pca$x)
        df_plot <- data.frame(PC = seq_along(pca$sdev^2/sum(pca$sdev^2)),
                              VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
                mutate(ELBOW = factor(ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE'),
                                      levels = c('UPPER','FORE')))
        p2 <- df_plot %>%  
                ggplot(aes(PC,VAR, color = ELBOW)) +
                geom_point() +
                geom_vline(xintercept = elbow_finder(df_plot$PC,df_plot$VAR)[1], linetype='dashed') +
                ggtitle(paste0(TISSUE,': Variance by PC')) +
                ylab('Variance')
        PC_meta <- t(pca$x) %>%
                as_tibble() %>%
                mutate(PC = seq_along(pca$sdev^2/sum(pca$sdev^2))) %>%
                mutate(VAR = pca$sdev^2/sum(pca$sdev^2)) %>%
                mutate(ELBOW = factor(ifelse(PC < elbow_finder(PC,VAR)[1], 'UPPER', 'FORE'), 
                                      levels = c('UPPER','FORE'))) %>%
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
                mutate(MAHALANOBIS_OUTLIER = ifelse(MAHALANOBIS_FDR <= 0.05,'Outlier', 'Normal')) %>%
                mutate(LOF = LOF(pca$x, log = F)) %>%
                left_join(y = col_data, by = 'sample_key') %>%
                left_join(y = join_meta, by = 'SD_DIST_PC') %>%
                arrange(dplyr::desc(SD_DIST))
        # Plot the PCA plot
        p1 <- PC %>%
                ggplot(aes(x = PC1, y = PC2, color = Key.anirandgroup)) +
                geom_point(size = 3) +
                coord_equal() +
                xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
                ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) + 
                ggtitle(paste0(TISSUE)) +
                scale_color_manual(values=ec_colors, name = "Training/Control Group")
        plot(p1)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-PCA-naive_steep.pdf'),
            width = 6, height = 4)
        plot(p1)
        dev.off()
        PC <- PC %>%
                mutate(OUT_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                                                  log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, ''))
        plot(p2)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-varPC_steep.pdf'),
            width = 6, height = 4)
        plot(p2)
        dev.off()
        
        p <- PC %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = Key.anirandgroup)) +
                geom_label_repel(aes(label=OUT_LABEL), colour = "black") +
                geom_point(size = 2) +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle(paste0(TISSUE,': Multidimensional Distance')) +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)') +
                scale_color_manual(values=ec_colors, name = "Training/Control Group")
        # Samples plotted against 2 multidimension distance metrics--local
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-outlier-distance-cohorts_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        # Examine the variance from of the PC that is associated with the outlier
        p <- PC %>%
                mutate(PC_LABEL = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                                                 log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), SD_DIST_PC, '')) %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = PC_VAR)) +
                geom_point(size = 2) +
                geom_label_repel(aes(label=PC_LABEL), color = 'black') +
                #geom_label_repel(aes(label=sample_key)) +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle(paste0(TISSUE,': Multidimensional Distance')) +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)') +
                labs(color='Variance per PC') 
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-outlier-distance-PC_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        
        mod <- lm(GINI ~ log(MAHALANOBIS_DIST), data = PC)
        print(summary(mod))
        mod <- lm(GINI ~ log(LOF), data = PC)
        print(summary(mod))
        mod <- lm(log(MAHALANOBIS_DIST) ~ log(LOF), data = PC)
        print(summary(mod))
        
        
        # Visualize the gini index
        p <- PC %>%
                ggplot(aes(x = log(MAHALANOBIS_DIST), y = log(LOF), color = GINI)) +
                geom_point(size = 2) +
                geom_label_repel(aes(label=OUT_LABEL), color = "black") +
                geom_vline(xintercept = tukey_mc_up(log(PC$MAHALANOBIS_DIST)), color = 'red') +
                geom_hline(yintercept = tukey_mc_up(log(PC$LOF)),  color = 'red') +
                ggtitle(paste0(TISSUE,': Multidimensional Distance')) +
                ylab('Local Outlier Factor (log)') +
                xlab('Mahalanobis Distance (log)') +
                labs(color='Gini Index') +
                scale_color_continuous(high = "#132B43", low = "#56B1F7",
                                       name = "Gini Index")
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-outlier-distance-gini_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        # Plot the Gini Index for these samples
        p <- ggplot() +
                geom_violin(data = PC, aes(x = TISSUE, y= GINI)) +
                geom_jitter(data = PC, aes(x = TISSUE, y= GINI, 
                                           color = log(MAHALANOBIS_DIST)), 
                            height = NULL, width = 0.15, 
                            alpha = 0.9, size = 2) +
                geom_label_repel(data = PC, aes(x = TISSUE, y= GINI,label=OUT_LABEL), color = 'black') +
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
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-gini-outliers_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        
        # PCA Mahalanobis
        p <- PC %>%
                mutate(OUTLIER = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                                                log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, '')) %>%
                ggplot(aes(x = PC1, y = PC2, color = log(MAHALANOBIS_DIST))) +
                geom_point(size = 3) +
                coord_equal() +
                geom_label_repel(aes(label=OUTLIER), color = "black") +
                xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
                ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) + 
                ggtitle(paste0(TISSUE)) +
                scale_color_continuous(high = "#132B43",low = "#56B1F7",
                                       name = "Mahalanobis Distance (log)")
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-PCA-naive-mahalanobis_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        # PCA Gini
        p <- PC %>%
                mutate(OUTLIER = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                                                log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, '')) %>%
                ggplot(aes(x = PC1, y = PC2, color = GINI)) +
                geom_point(size = 3) +
                coord_equal() +
                geom_label_repel(aes(label=OUTLIER), color = "black") +
                xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
                ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) + 
                ggtitle(paste0(TISSUE)) +
                scale_color_continuous(high = "#132B43", low = "#56B1F7",
                                       name = "Gini Index")
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-PCA-naive-gini_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        # Training/Control Group
        p <- PC %>%
                mutate(OUTLIER = ifelse(log(LOF) > tukey_mc_up(log(PC$LOF)) | 
                                                log(MAHALANOBIS_DIST) > tukey_mc_up(log(PC$MAHALANOBIS_DIST)), sample_key, '')) %>%
                ggplot(aes(x = PC1, y = PC2, color = Key.anirandgroup)) +
                geom_point(size = 3) +
                coord_equal() +
                geom_label_repel(aes(label=OUTLIER), color = "black") +
                xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
                ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) + 
                ggtitle(paste0(TISSUE)) +
                scale_color_manual(values=ec_colors, name = "Training/Control Group")
        plot(p)
        pdf(paste0(WD,'/plots/20200718_pass1b-rnaseq-',tissue,'-PCA-naive-cohorts-labeled_steep.pdf'),
            width = 6, height = 4)
        plot(p)
        dev.off()
        
}

#' ## Qualitative Assessment of Variance in Tissue

# Qualitative Assessment of Variance in Tissue
################################################################################
#####     Qualitative Assessment of Variance in Tissue      ####################
################################################################################

# Examine the Annotations most strongly associated with PCs
################################################################################

# Annotations to exclude
exclude <- c('Training.day40_posttrainlact', 'Registration.cagenumber', 'Training.day29_score',
             'Training.day36time','NMR.Testing.nmr_weight_1','Terminal.Weight.sol',
             'Familiarization.weight','Registration.weight',
             'Training.day2_weight','Training.day3_weight',
             'Training.day4_weight','Training.day5_weight','Training.day6_weight',
             'Training.day7_weight','Training.day8_weight','Training.day9_weight',
             'Training.day10_weight','Training.day11_weight','Training.day12_weight',
             'Training.day13_weight','Training.day14_weight','Training.day15_weight',
             'Training.day16_weight','Training.day17_weight','Training.day18_weight',
             'Training.day19_weight','Training.day20_weight','Training.day21_weight',
             'Training.day22_weight','Training.day23_weight','Training.day24_weight',
             'Training.day25_weight','Training.day26_weight','Training.day27_weight',
             'Training.day28_weight','Training.day29_weight','Training.day30_weight',
             'Training.day31_weight','Training.day32_weight','Training.day33_weight',
             'Training.day34_weight','Training.day35_weight','Training.day36_weight',
             'Training.day37_weight','Training.day38_weight','Training.day39_weight',
             'Training.day2date','Training.day3date',
             'Training.day4date','Training.day5date','Training.day6date',
             'Training.day7date','Training.day8date','Training.day9date',
             'Training.day10date','Training.day11date','Training.day12date',
             'Training.day13date','Training.day14date','Training.day15date',
             'Training.day16date','Training.day17date','Training.day18date',
             'Training.day19date','Training.day20date','Training.day21date',
             'Training.day22date','Training.day23date','Training.day24date',
             'Training.day25date','Training.day26date','Training.day27date',
             'Training.day28date','Training.day29date','Training.day30date',
             'Training.day31date','Training.day32date','Training.day33date',
             'Training.day34date','Training.day35date','Training.day36date',
             'Training.day37date','Training.day38date','Training.day39date',
             'Training.day2_score','Training.day3_score',
             'Training.day4_score','Training.day5_score','Training.day6_score',
             'Training.day7_score','Training.day8_score','Training.day9_score',
             'Training.day10_score','Training.day11_score','Training.day12_score',
             'Training.day13_score','Training.day14_score','Training.day15_score',
             'Training.day16_score','Training.day17_score','Training.day18_score',
             'Training.day19_score','Training.day20_score','Training.day21_score',
             'Training.day22_score','Training.day23_score','Training.day24_score',
             'Training.day25_score','Training.day26_score','Training.day27_score',
             'Training.day28_score','Training.day29_score','Training.day30_score',
             'Training.day31_score','Training.day32_score','Training.day33_score',
             'Training.day34_score','Training.day35_score','Training.day36_score',
             'Training.day37_score','Training.day38_score','Training.day39_score',
             'Training.day1_comments','Training.day2_comments','Training.day3_comments',
             'Training.day4_comments','Training.day5_comments','Training.day6_comments',
             'Training.day7_comments','Training.day8_comments','Training.day9_comments',
             'Training.day10_comments','Training.day11_comments','Training.day12_comments',
             'Training.day13_comments','Training.day14_comments','Training.day15_comments',
             'Training.day16_comments','Training.day17_comments','Training.day18_comments',
             'Training.day19_comments','Training.day20_comments','Training.day21_comments',
             'Training.day22_comments','Training.day23_comments','Training.day24_comments',
             'Training.day25_comments','Training.day26_comments','Training.day27_comments',
             'Training.day28_comments','Training.day29_comments','Training.day30_comments',
             'Training.day31_comments','Training.day32_comments','Training.day33_comments',
             'Training.day34_comments','Training.day35_comments','Training.day36_comments',
             'Training.day37_comments','Training.day38_comments','Training.day39_comments',
             'Training.day40_comments',
             'Training.day2_posttrainlact','Training.day3_posttrainlact',
             'Training.day4_posttrainlact','Training.day5_posttrainlact',
             'Training.day6_posttrainlact','Training.day7_posttrainlact',
             'Training.day8_posttrainlact','Training.day9_posttrainlact',
             'Training.day10_posttrainlact','Training.day11_posttrainlact',
             'Training.day12_posttrainlact','Training.day13_posttrainlact',
             'Training.day14_posttrainlact','Training.day15_posttrainlact',
             'Training.day16_posttrainlact','Training.day17_posttrainlact',
             'Training.day18_posttrainlact','Training.day19_posttrainlact',
             'Training.day20_posttrainlact','Training.day21_posttrainlact',
             'Training.day22_posttrainlact','Training.day23_posttrainlact',
             'Training.day24_posttrainlact','Training.day25_posttrainlact',
             'Training.day26_posttrainlact','Training.day27_posttrainlact',
             'Training.day28_posttrainlact','Training.day29_posttrainlact',
             'Training.day30_posttrainlact','Training.day31_posttrainlact',
             'Training.day32_posttrainlact','Training.day33_posttrainlact',
             'Training.day34_posttrainlact','Training.day35_posttrainlact',
             'Training.day36_posttrainlact','Training.day37_posttrainlact',
             'Training.day38_posttrainlact','Training.day39_posttrainlact')

# Find variables associated with PCs
pc_cor_df <- cor_PC_1_6(rld = rld, ntop = 20000, intgroups = names(tod_cols)) %>%
        filter(Adjusted_R_Sq > 0.3) %>%
        filter(Condition %!in% exclude) %>%
        arrange(PC)
# Examine Variables of interest for PCs 1 - 4
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


# v <- pc_cor_df %>%
#         filter(PC %in% c(1:2)) %>%
#         select(Condition) %>% dplyr::slice(8) %>% 
#         unlist() %>% as.character()
# v <- 'VO2.Max.Test.d_vo2_2'
# title_parts <- pc_cor_df %>%
#         filter(Condition == v)
# title <- paste0(title_parts$Condition,'\nPC = ',
#                 title_parts$PC,'; R2 = ',
#                 round(title_parts$Adjusted_R_Sq,3))
# p <- DESeq2::plotPCA(rld, intgroup =v, ntop = 20000) +
#         ggtitle(title) +
#         coord_equal()
# plot(p)

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
                        mis_id_var <- 'Registration.sex'
                }
                #mis_id_var <- MIS_ID_VAR[1]
                print(mis_id_var)
                # Examine metadata variables most correlated with suspected misidentified samples to determine if there is an effect are genuine
                # cor_out_df <- cor_outlier2(rld = rld, ntop = 20000, intgroups = names(colData(rld)))
                # cor_out_df <- cor_out_df %>% 
                #         filter(Condition != 'outlier')
                # cor_out_df
                
                # Show a PCA plot labeled with suspected misidentified samples
                mypar()
                pcaData <- DESeq2::plotPCA(rld, 
                                           intgroup=c(mis_id_var, 'sample_key'), 
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
                        #scale_color_manual(values=ec_colors) +
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
        pc_cor_df <- cor_PC_1_4(rld = rld, ntop = 20000, intgroups = names(colData(rld))) %>%
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
                p <- DESeq2::plotPCA(rld, intgroup =pri_var) +
                        guides(color=guide_legend(title=pri_var))
                plot(p)
                # After
                p <- DESeq2::plotPCA(rld_final, intgroup = pri_var, ntop = 500) +
                        guides(color=guide_legend(title=pri_var))
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



