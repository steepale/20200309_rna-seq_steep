#'---
#' title: "PASS1A Rat Tissue: -- "
#' author: "Alec Steep"
#' date: "20200823"
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
#BiocManager::install('edgeR')
#install.packages('tidyverse')

# Load dependencies
pacs...man <- c('tidyverse','GenomicRanges', 'DESeq2','devtools','rafalib','GO.db','vsn','hexbin','ggplot2', 'GenomicFeatures','Biostrings','BSgenome','AnnotationHub','plyr','dplyr', 'org.Rn.eg.db','pheatmap','sva','formula.tools','pathview','biomaRt', 'PROPER','SeqGSEA','purrr','BioInstaller','RColorBrewer','lubridate', 'hms','ggpubr', 'ggrepel','genefilter','qvalue','ggfortify','som', 'vsn','org.Mm.eg.db','VennDiagram','EBImage','reshape2','xtable','kohonen','som','caret','enrichR','gplots','tiff','splines','gam','EnhancedVolcano','mvoutlier','multtest','bigutilsr','bigstatsr','magrittr','corrplot','edgeR','gridExtra','ggdendro')
lapply(pacs...man, FUN = function(X) {
        do.call('library', list(X)) })

#' ## Load Custom Functions

#+ Functions, message=FALSE, results='hide', warning = FALSE
################################################################################
######################### Functions ############################################
################################################################################

# Resolve conflicts
counts <- DESeq2::counts
map <- purrr::map
summarize <- dplyr::summarize
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
select <- dplyr::select
spread <- tidyr::spread

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

# Load Data files
corr_min_tissue_file <- paste0(WD,'/docs/20200823_rnaseq-pass1a-outlier-corr-mins_steep.txt')
cm_df <- read.table(corr_min_tissue_file, sep = '\t', header = T)

# Files last saved in: 20200309_exploration-rna-seq-phase1_steep.R
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

#' ## Declare Variables

#+ Declare Variables
################################################################################
#####     Declare Variables     ################################################
################################################################################
#' #### Section to generate a table of major decisions

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

# Test Conditions
TISSUE <- 'Testes'
ZeroGenes <- "NFilter"
TForm <- 'CPM'
Corr <- 'Pearson'
BatchOrder <- 'YBatch'
TForder <- 'F_T'
PrintPlot <- 'All'

# For loops
# for(TForm in c('VST','Log2','CPM')){
# for(ZeroGenes in c('YFilter','NFilter')){
# for(TForder in c('T_F','F_T')){
# for(Corr in c('Pearson','Spearman')){
#for(BatchOrder in c('YBatch','NBatch')){

# Condition to make development easier
DEV <- 'Yes'

# All tissues
all_tissue_df <- data.frame()

#microbenchmark(
for(TISSUE in c('Adrenal', 'Aorta', 'Brown Adipose', 'Cortex', 'Gastrocnemius',
          'Heart', 'Hippocampus', 'Hypothalamus', 'Kidney','Liver','Lung','Ovaries','Spleen', 'Testes', 'White Adipose')){
        #for(TISSUE in c('Adrenal')){
        # i is tissue specific
        #for(PrintPlot in c('TFormvsBatch','TFormvsFilter','TFormvsTF','TFormvsCorr','All')){
        for(PrintPlot in c('All')){
                i = 1
                # Lists of plots
                elbow_list <- list()
                elbow_final_list <- list()
                heatmap_list <- list()
                NNpca_list <- list()
                NPpca_list <- list()
                NPpca_batch_list <- list()
                NPpca_noo_list <- list()
                
                # Lists for dataframe construction
                SxW_hm_list_r <- list()
                SxW_hm_list_s <- list()
                SxW_hm_list_w <- list()
                
                # Operations
                TForm_vec <- c('VST','Log2','CPM')
                BatchOrder_vec <- c('YBatch','NBatch')
                ZeroGenes_vec <- c('NFilter','YFilter')
                TForder_vec <- c('F_T','T_F')
                Corr_vec <- c('Pearson','Spearman') 
                if(PrintPlot == 'TFormvsBatch'){
                        ZeroGenes_vec <- 'NFilter'
                        TForder_vec <- 'F_T'
                        Corr_vec <- 'Pearson' 
                }else if(PrintPlot == 'TFormvsFilter'){
                        BatchOrder_vec <- 'YBatch'
                        ZeroGenes_vec <- c('NFilter','YFilter')
                        TForder_vec <- 'F_T'
                        Corr_vec <- 'Pearson'
                }else if(PrintPlot == 'TFormvsTF'){
                        BatchOrder_vec <- 'YBatch'
                        ZeroGenes_vec <- 'NFilter'
                        TForder_vec <- c('F_T','T_F')
                        Corr_vec <- 'Pearson'
                }else if(PrintPlot == 'TFormvsCorr'){
                        BatchOrder_vec <- 'YBatch'
                        ZeroGenes_vec <- 'NFilter'
                        TForder_vec <- 'F_T'
                        Corr_vec <- c('Pearson','Spearman') 
                }else if(PrintPlot == 'All'){
                        TForm_vec <- c('VST','Log2','CPM')
                        BatchOrder_vec <- c('YBatch','NBatch')
                        ZeroGenes_vec <- c('NFilter','YFilter')
                        TForder_vec <- c('F_T','T_F')
                        Corr_vec <- c('Pearson','Spearman') 
                }
                for(ZeroGenes in ZeroGenes_vec){
                        for(TForder in TForder_vec){
                                for(Corr in Corr_vec){
                                        for(BatchOrder in BatchOrder_vec){
                                                for(TForm in TForm_vec){
                                                        print(i)
                                                        
                                                        # Declare Outliers
                                                        if(TISSUE == 'Kidney'){
                                                                TIS <- 'KID'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Liver'){
                                                                TIS <- 'LIV'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Hypothalamus'){
                                                                TIS <- 'SCN'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Aorta'){
                                                                TIS <- 'AOR'
                                                                PRI_VAR <- 'None'
                                                        }else if(TISSUE == 'Gastrocnemius'){
                                                                TIS <- 'SKM'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Heart'){
                                                                TIS <- 'HAT'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Adrenal'){
                                                                TIS <- 'ADG'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Brown Adipose'){
                                                                TIS <- 'BAT'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'White Adipose'){
                                                                TIS <- 'WAT'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Cortex'){
                                                                TIS <- 'COR'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }else if(TISSUE == 'Hippocampus'){
                                                                TIS <- 'HIP'
                                                                PRI_VAR <- 'None'
                                                        }else if(TISSUE == 'Lung'){
                                                                TIS <- 'LUNG'
                                                                PRI_VAR <- c('animal.key.batch','animal.registration.sex')
                                                        }else if(TISSUE == 'Ovaries'){
                                                                TIS <- 'OVR'
                                                                PRI_VAR <- 'None'
                                                        }else if(TISSUE == 'Testes'){
                                                                TIS <- 'TES'
                                                                PRI_VAR <- 'None'
                                                        }else if(TISSUE == 'Spleen'){
                                                                TIS <- 'SPL'
                                                                PRI_VAR <- 'animal.registration.sex'
                                                        }
                                                        
order_of_operations <- paste0(TIS,'-',PrintPlot,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)
print(paste0(PrintPlot,': ',order_of_operations))

#' ## Collect Samples of Interest and Normalize

#+ Collect Samples of Interest and Normalize
################################################################################
#####     Collect Samples of Interest and Normalize      #######################
################################################################################
if(TISSUE == c('Gastrocnemius')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Gastrocnemius') %>%
                filter(Seq_batch == 'MSSM_1') %>%
                filter(!is.na(animal.registration.sex))
}else if(TISSUE == c('Lung')){
        tod_cols <- col_data %>%
                filter(Tissue == 'Lung') %>%
                filter(Seq_batch == 'MSSM_3') %>%
                filter(!is.na(animal.registration.sex))
}else{
        # Filter Samples (meta)
        tod_cols <- col_data %>%
                filter(Tissue == TISSUE) %>%
                filter(!is.na(animal.registration.sex))
}
rownames(tod_cols) <- tod_cols$sample_key

# Collect samples without NA values in TOD
# nona_sams <- tod_cols %>%
#         filter(!is.na(specimen.collection.t_death_hour)) %>%
#         filter(!is.na(animal.registration.sex)) %>%
#         select(sample_key) %>% unlist() %>% as.character()
# SPEED CHECK: No need for this command, reduced for speed

# Collect tissue specific counts
tod_counts <- count_data[,rownames(tod_cols)]

#' ##### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(tod_cols) == colnames(tod_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
dds <- DESeqDataSetFromMatrix(countData = tod_counts,
                               colData = tod_cols,
                               design = ~ 1)

################################################################################
####################### Filter Low Count Genes  ################################
################################################################################
if(TForder == 'T_F'){
        # Transform
        ########################################################################
        if(TForm == 'VST'){
                #' #### Transform:
                rld <- DESeq2::vst(dds, blind = FALSE)
                cmat <- assay(rld)
        }else if(TForm == 'Log2'){
                cmat <- log2(assay(dds, normalized = T) + 1)
        }else if(TForm == 'CPM'){
                # filter
                cmat <- cpm(assay(dds), log = T)
        }
        # Filter
        ########################################################################
        if(ZeroGenes == 'NFilter'){
                #print('Keeping all genes')  
                genes_removed_n <- 0
                genes_removed_p <- '0.00%'
        }else if(ZeroGenes == 'YFilter'){
                reads_n <- 1
                #print(paste0('Remove genes with an average sequencing depth (across samples) of less than ',reads_n))
                keep <- rowSums(cmat)/ncol(cmat) >= reads_n
                genes_removed_n <- as.numeric(table(keep)[1])
                genes_removed_p <- paste0(as.character(round(genes_removed_n/nrow(tod_counts),4)*100),'%')
                genes_kept_n <- as.numeric(table(keep)[2])
                genes_kept_p <- paste0(as.character(round(genes_kept_n/nrow(tod_counts),4)*100),'%')
                #print(paste0('Removed ',genes_removed_n,' of ',nrow(tod_counts),' genes (',genes_removed_p,')'))
        }
}
if(TForder == 'F_T'){
        #print('Transformaing counts after filtering genes (if any)')
if(ZeroGenes == 'NFilter'){
   #print('Keeping all genes')  
   genes_removed_n <- 0
   genes_removed_p <- '0.00%'
}else if(ZeroGenes == 'YFilter'){
        reads_n <- 1
        #print(paste0('Remove genes with an average sequencing depth (across samples) of less than ',reads_n))
        keep <- rowSums(counts(dds))/ncol(dds) >= reads_n
        genes_removed_n <- as.numeric(table(keep)[1])
        genes_removed_p <- paste0(as.character(round(genes_removed_n/nrow(tod_counts),4)*100),'%')
        genes_kept_n <- as.numeric(table(keep)[2])
        genes_kept_p <- paste0(as.character(round(genes_kept_n/nrow(tod_counts),4)*100),'%')
        #print(paste0('Removed ',genes_removed_n,' of ',nrow(tod_counts),' genes (',genes_removed_p,')'))
        dds <- dds[keep,]
}

#' #### Size Factor Estimates (summary across samples):
#' Size factors are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
# estimateSizeFactors() gives us a robust estimate in sequencing depth.
dds <- estimateSizeFactors(dds)
#summary(sizeFactors(dds))

################################################################################
####################### Transformation Technique ###############################
################################################################################
if(TForm == 'VST'){
        # Transform:
        #print('VST Transformation')
        rld <- DESeq2::vst(dds, blind = FALSE)
        cmat <- assay(rld)
}else if(TForm == 'Log2'){
        #print('Log2 Transformation')
        cmat <- log2(assay(dds, normalized = T) + 1)
}else if(TForm == 'CPM'){
        #print('CPM transformation (Log2)')
        cmat <- cpm(assay(dds), log = T)
}
}
################################################################################
####################### Batch Correction #######################################
################################################################################
if(BatchOrder == 'YBatch'){
        if(PRI_VAR == 'None'){
                #print('No adjustment for batch effect')
        }else if(PRI_VAR != 'None'){
                #print(paste0('Adjusting batch effect for ',PRI_VAR))
                if(TISSUE != 'Lung'){
                        cmat <- limma::removeBatchEffect(cmat, dds[[PRI_VAR]])
                }else if(TISSUE == 'Lung'){
                        cmat <- limma::removeBatchEffect(cmat, dds[[PRI_VAR[1]]])
                        cmat <- limma::removeBatchEffect(cmat, dds[[PRI_VAR[2]]])
                }  
        }
}

################################################################################
####################### Sample by Sample Correlations ##########################
################################################################################
# Create an NxN Matrix of correlation
if(Corr == 'Pearson'){
        res <- stats::cor(cmat, method = 'pearson') %>% round(digits = 3)
}else if(Corr == 'Spearman'){
        res <- stats::cor(cmat, method = 'spearman') %>% round(digits = 3)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
        cormat[lower.tri(cormat)]<- NA
        return(cormat)
}
# Reorder the correlation matrix
reorder_cormat <- function(cormat){
        # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(res)

################################################################################
####### Determine the Suspected Outlier Samples for each tissue #################
################################################################################
# Tissue specific
if(TISSUE == 'Adrenal'){
        sus_samples = c()
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Aorta"){
        sus_samples = c()
        out_samples = c('90159016502_SN2')
        either_samples = c('90041016502_SN2')
}else if(TISSUE == "Brown Adipose"){
        sus_samples = c()
        out_samples = c('90128016905_SF1','90047016905_SF1')
        either_samples = c()
}else if(TISSUE == "Cortex"){
        sus_samples = c()
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Gastrocnemius"){
        sus_samples = c('90139015502_SN1','90045015502_SN1')
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Heart"){
        sus_samples = c()
        out_samples = c()
        either_samples = c('90133015802_SN1', '90109015802_SN1')
}else if(TISSUE == "Hippocampus"){
        sus_samples = c()
        out_samples = c('90040015202_SF2', '90115015202_SF2', '90041015202_SF2',
                        '90009015202_SF2')
        either_samples = c()
}else if(TISSUE == "Hypothalamus"){
        sus_samples = c('90150015402_SF2')
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Kidney"){
        sus_samples = c()
        out_samples = c()
        either_samples = c('90146015902_SN1')
}else if(TISSUE == "Liver"){
        sus_samples = c('90142016803_SF1', '90031016803_SF1')
        out_samples = c()
        either_samples = c('90042016803_SF1')
}else if(TISSUE == "Lung"){
        sus_samples = c()
        out_samples = c()
        either_samples = c('90041016604_SN3','90143016604_SN3')
}else if(TISSUE == "Ovaries"){
        sus_samples = c()
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Spleen"){
        sus_samples = c('90031016202_SF2')
        out_samples = c()
        either_samples = c()
}else if(TISSUE == "Testes"){
        sus_samples = c('90129016302_SN2')
        out_samples = c('90017016302_SN2','90031016302_SN2')
        either_samples = c()
}else if(TISSUE == "White Adipose"){
        sus_samples = c()
        out_samples = c()
        either_samples = c()
}

################################################################################
####################### Create Elbow Plot of Correlation Matrixes ##############
################################################################################
# Create the dataframe for the elbow plot
df_plot <- data.frame(sample_key = names(sort(rowMeans(cormat))),
                      Index = seq_along(sort(rowMeans(cormat))),
                      r = sort(rowMeans(cormat))) %>%
        mutate(nr =  r * -1) %>%
        mutate(ELBOW = ifelse(Index < elbow_finder(Index,nr)[1], 'Below', 'Above')) %>%
        mutate(OUTLIER = ifelse(ELBOW == 'Below', sample_key, '')) %>%
        mutate(OUTLIER_SR = case_when(sample_key %in% sus_samples ~ 'Suspected',
                                      sample_key %in% either_samples ~ 'Bubble',
                                      sample_key %in% out_samples ~ 'Outlier',
                                      sample_key %!in% c(sus_samples, either_samples, out_samples) ~ 'Normal')) %>%
        mutate(OUTLIER_SR = factor(OUTLIER_SR, levels = c('Normal', 'Suspected', 'Bubble', 'Outlier')))

# Final elbow plot: less samples labeled
elbow_final_list[[i]] <- df_plot %>%
        ggplot(aes(Index,r, color = OUTLIER_SR)) +
        geom_label_repel(aes(label=ifelse(OUTLIER_SR != 'Normal', sample_key, '')),
                         size = 1.8,
                         nudge_x      = 10,
                         direction    = "y",
                         hjust        = 0,
                         segment.size = 0.2) +
        geom_point(size = 2) +
        geom_vline(xintercept = elbow_finder(df_plot$Index,df_plot$nr)[1], linetype='dashed') +
        ylab('Correlation Value') +
        xlab('Sample Index') +
        ylim(0.65,1) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) +
        theme_bw() +
        theme(plot.title = element_text(size=6)) +
        theme(legend.title = element_text(size = 4), 
              legend.text = element_text(size = 4)) +
        labs(color = "Outlier_Status") +
        theme(legend.key.size = unit(0.5, 'lines')) +
        theme(legend.key.width = unit(0.5, 'lines')) +
        theme(axis.title.x = element_text(size = rel(0.75))) +
        theme(axis.title.y = element_text(size = rel(0.75))) +
        theme(aspect.ratio = 1)

elbow_list[[i]] <- df_plot %>%
        ggplot(aes(Index,r, color = ELBOW)) +
        geom_label_repel(aes(label=OUTLIER),
                         size = 1.8,
                         nudge_x      = 10,
                         direction    = "y",
                         hjust        = 0,
                         segment.size = 0.2) +
        geom_point(size = 2) +
        geom_vline(xintercept = elbow_finder(df_plot$Index,df_plot$nr)[1], linetype='dashed') +
        ylab('Correlation Value') +
        xlab('Sample Index') +
        ylim(0.65,1) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) +
        theme_bw() +
        theme(plot.title = element_text(size=6)) +
        theme(legend.position="none") +
        theme(axis.title.x = element_text(size = rel(0.75))) +
        theme(axis.title.y = element_text(size = rel(0.75))) +
        theme(aspect.ratio = 1)

################################################################################
####################### Create Heatmaps of Correlation Matrices ################
################################################################################

# Melt the correlation matrix
melted_cormat <- melt(cormat, na.rm = TRUE)
melted_cormat <- melted_cormat %>%
        mutate(Vial1 = substr(Var1, start = 1, stop = 6)) %>%
        mutate(Vial2 = substr(Var2, start = 1, stop = 6))
melted_cormat$Vial1 <- as.factor(melted_cormat$Vial1)
melted_cormat$Vial2 <- as.factor(melted_cormat$Vial2)
melted_cormat$Tissue <- TISSUE
melted_cormat$Transform <- TForm
melted_cormat$Correlation <- Corr
melted_cormat$Gene_Filter <- ZeroGenes
melted_cormat$Filter_Order <- TForder
melted_cormat$value <- round(melted_cormat$value, digits = 3)

# Tissue and iteration specific minimum correlation values (manual annotation)
min_corr <- cm_df %>%
        filter(Tissue == TISSUE) %>%
        filter(Print_Plot == PrintPlot) %>%
        select(Corr_Min) %>%
        unlist() %>% min()

# Tissue and iteration specific mean correlation values (manual annotation)
min_corr <- cm_df %>%
        filter(Tissue == TISSUE) %>%
        filter(Print_Plot == PrintPlot) %>%
        select(Corr_Min) %>%
        unlist() %>% min()

# Set the breaks
brks <- round(seq(min_corr, 1, length.out=(1-min_corr)*30), digits = 3)
# Set the colors
col <- colorRampPalette(c("blue", "white", "red"))((1-min_corr)*1000)
# Plot the heatmaps
heatmap_list[[i]] <- melted_cormat %>%
        ggplot(aes(x = Var2,
                   y = Var1,
                   fill = value)) +
        geom_tile(color = "gray") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                         size = 1, hjust = 1)) +
        theme(axis.text.y = element_text(size = 3)) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        #guides(fill = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 4),
              legend.text = element_text(size = 4)) +
        theme(legend.key.size = unit(0.5, 'lines')) +
        theme(legend.key.width = unit(0.5, 'lines')) +
        scale_fill_gradientn(colors = col,
                              name="Correlation",
                             limits=c(min_corr,1),
                             breaks = brks) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) +
        theme(plot.title = element_text(size=6)) +
        coord_fixed() +
        theme(aspect.ratio=1)

################################################################################
################# Tissue-Specific Heatmap (Sample x Workflow) ##################
################################################################################

# Store the heatmap values, samples, and workflows in lists to be stroed in final dataframe later (faster)
SxW_hm_list_r[[i]] <- as.numeric(rowMeans(cormat))
SxW_hm_list_s[[i]] <- names(rowMeans(cormat))
SxW_hm_list_w[[i]] <- rep(order_of_operations, nrow(cormat))


################################################################################
########################## PCA Analysis (via NxN) ##############################
################################################################################

# Plot a PCA with outliers labeled
pca <- prcomp(cormat, scale. = F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC <- pca$x %>% as.data.frame()
PC$sample_key <- as.character(row.names(pca$x))
PC <- left_join(PC, tod_cols, by = 'sample_key')
NNpca_list[[i]] <- PC %>%
        mutate(OUTLIER = ifelse(sample_key %in% c(out_samples,either_samples,sus_samples), sample_key, '')) %>%
        mutate(OUTLIER_SR = case_when(sample_key %in% sus_samples ~ 'Suspected',
                                      sample_key %in% either_samples ~ 'Bubble',
                                      sample_key %in% out_samples ~ 'Outlier',
                                      sample_key %!in% c(sus_samples, either_samples, out_samples) ~ 'Normal')) %>%
        mutate(OUTLIER_SR = factor(OUTLIER_SR, levels = c('Normal', 'Suspected', 'Bubble', 'Outlier'))) %>%
        ggplot(aes(x = PC1, y = PC2, color=animal.key.anirandgroup, shape = OUTLIER_SR)) +
        geom_label_repel(aes(label=ifelse(OUTLIER_SR != 'Normal', sample_key, '')), 
                         size = 1.5, min.segment.length = unit(0, 'lines'),
                         nudge_y = .02) +
        geom_point(size = 2) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) + 
        theme(plot.title = element_text(size=6)) +
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) +
        scale_color_manual(values=ec_colors) +
        guides(color = guide_legend(override.aes = list(size = 1))) +
        guides(shape = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6)) +
        theme(legend.title=element_blank()) +
        theme(axis.title.x = element_text(size = rel(0.5))) +
        theme(axis.title.y = element_text(size = rel(0.5))) +
        theme(axis.text.x = element_text(size = 4)) +
        theme(axis.text.y = element_text(size = 4)) +
        theme(aspect.ratio=1) +
        theme(legend.position = c(1.4, 0.5)) +
        theme(legend.key.size = unit(0.75, 'lines')) +
        theme(plot.margin = unit(c(1,3,1,1), "cm"))


################################################################################
########################## PCA Analysis (via NxP) ##############################
################################################################################

# Plot a PCA with outliers labeled
pca <- prcomp(t(cmat), scale. = F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC <- pca$x %>% as.data.frame()
PC$sample_key <- as.character(row.names(pca$x))
PC <- left_join(PC, tod_cols, by = 'sample_key')
NPpca_list[[i]] <- PC %>%
        mutate(OUTLIER = ifelse(sample_key %in% c(out_samples,either_samples,sus_samples), sample_key, '')) %>%
        mutate(OUTLIER_SR = case_when(sample_key %in% sus_samples ~ 'Suspected',
                                      sample_key %in% either_samples ~ 'Bubble',
                                      sample_key %in% out_samples ~ 'Outlier',
                                      sample_key %!in% c(sus_samples, either_samples, out_samples) ~ 'Normal')) %>%
        mutate(OUTLIER_SR = factor(OUTLIER_SR, levels = c('Normal', 'Suspected', 'Bubble', 'Outlier'))) %>%
        ggplot(aes(x = PC1, y = PC2, color=animal.key.anirandgroup, shape = OUTLIER_SR)) +
        # ggforce::geom_mark_ellipse(aes(x=PC1, y=PC2,fill=animal.registration.sex, group=animal.registration.sex), color = 'black', alpha = 0.2) +
        geom_label_repel(aes(label=ifelse(OUTLIER_SR != 'Normal', sample_key, '')), 
                         size = 1.5, min.segment.length = unit(0, 'lines'),
                         nudge_y = .04) +
        geom_point(size = 2) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) + 
        theme(plot.title = element_text(size=6)) +
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) +
        #ggtitle(paste0('PCA of ',TISSUE)) +
        scale_color_manual(values=ec_colors) +
        guides(shape = guide_legend(override.aes = list(size = 1))) +
        guides(color = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6)) +
        theme(legend.title=element_blank()) +
        theme(axis.title.x = element_text(size = rel(0.5))) +
        theme(axis.title.y = element_text(size = rel(0.5))) +
        theme(axis.text.x = element_text(size = 4)) +
        theme(axis.text.y = element_text(size = 4)) +
        theme(aspect.ratio=1) +
        theme(legend.position = c(1.4, 0.5)) +
        theme(legend.key.size = unit(0.75, 'lines')) +
        theme(plot.margin = unit(c(1,3,1,1), "cm"))

# Plot the same PCA but with sex-annotation
if(TISSUE == 'Lung'){
        colcol = 'animal.key.batch'
}else{
        colcol = 'animal.registration.sex'
}

NPpca_batch_list[[i]] <- PC %>%
        mutate(OUTLIER = ifelse(sample_key %in% c(out_samples,either_samples,sus_samples), sample_key, '')) %>%
        mutate(OUTLIER_SR = case_when(sample_key %in% sus_samples ~ 'Suspected',
                                      sample_key %in% either_samples ~ 'Bubble',
                                      sample_key %in% out_samples ~ 'Outlier',
                                      sample_key %!in% c(sus_samples, either_samples, out_samples) ~ 'Normal')) %>%
        mutate(OUTLIER_SR = factor(OUTLIER_SR, levels = c('Normal', 'Suspected', 'Bubble', 'Outlier'))) %>%
        ggplot(aes_string(x = "PC1", y = "PC2", color=colcol, shape = "OUTLIER_SR")) +
        # ggforce::geom_mark_ellipse(aes(x=PC1, y=PC2,fill=animal.registration.sex, group=animal.registration.sex), color = 'black', alpha = 0.2) +
        geom_label_repel(aes(label=ifelse(OUTLIER_SR != 'Normal', sample_key, '')), 
                         size = 1.5, min.segment.length = unit(0, 'lines'),
                         nudge_y = .04) +
        geom_point(size = 2) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) + 
        #theme(plot.title = element_text(face = "bold")) +
        theme(plot.title = element_text(size=6)) +
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) +
        #ggtitle(paste0('PCA of ',TISSUE)) +
        guides(shape = guide_legend(override.aes = list(size = 1))) +
        guides(color = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6)) +
        theme(legend.title=element_blank()) +
        theme(axis.title.x = element_text(size = rel(0.5))) +
        theme(axis.title.y = element_text(size = rel(0.5))) +
        theme(axis.text.x = element_text(size = 4)) +
        theme(axis.text.y = element_text(size = 4)) +
        theme(aspect.ratio=1) +
        theme(legend.position = c(1.4, 0.5)) +
        theme(legend.key.size = unit(0.75, 'lines')) +
        theme(plot.margin = unit(c(1,3,1,1), "cm"))


################################################################################
########################## PCA Analysis (via NxP) W/o Outliers #################
################################################################################

# Only apply if there are actually outliers
# Plot the samples without the outliers
in_samples <- colnames(cmat)[colnames(cmat) %!in% out_samples]
pca <- prcomp(t(cmat[, in_samples]), scale. = F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
PC <- pca$x %>% as.data.frame()
PC$sample_key <- as.character(row.names(pca$x))
PC <- left_join(PC, tod_cols, by = 'sample_key')
NPpca_noo_list[[i]] <- PC %>%
        filter(sample_key %!in% out_samples) %>%
        ggplot(aes(x = PC1, y = PC2, color=animal.key.anirandgroup)) +
        geom_point(size = 2) +
        xlab(paste0('PC1: ',round(percentVar[1]*100,2),'% variance')) +
        ylab(paste0('PC2: ',round(percentVar[2]*100,2),'% variance')) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) + 
        theme(plot.title = element_text(size=6)) +
        scale_color_manual(values=ec_colors) +
        guides(shape = guide_legend(override.aes = list(size = 1))) +
        guides(color = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 8), 
              legend.text = element_text(size = 6)) +
        theme(legend.title=element_blank()) +
        theme(axis.title.x = element_text(size = rel(0.5))) +
        theme(axis.title.y = element_text(size = rel(0.5))) +
        theme(axis.text.x = element_text(size = 4)) +
        theme(axis.text.y = element_text(size = 4)) +
        theme(aspect.ratio=1) +
        theme(legend.position = c(1.4, 0.5)) +
        theme(legend.key.size = unit(0.75, 'lines')) +
        theme(plot.margin = unit(c(1,3,1,1), "cm"))

################################################################################
########################## Create a Final Report ###############################
################################################################################

solo_tissue_df <- data.frame(Tissue = TISSUE,
                             Transformation = TForm,
                             Filter_Genes = ZeroGenes,
                             Total_Genes = length(keep),
                             Genes_Filtered_N = genes_removed_n,
                             Genes_Filtered_P = genes_removed_p,
                             Correlation_Metric = Corr,
                             Batch_Adjustment = paste(as.character(PRI_VAR), collapse=";"),
                             Batch_Correction = BatchOrder,
                             Suspected_N = length(sus_samples),
                             Suspected_Samples = paste(as.character(sus_samples), collapse=";"),
                             Outlier_N = length(out_samples),
                             Outlier_Samples = paste(as.character(out_samples), collapse=";"),
                             Order = order_of_operations,
                             Print_Plot = PrintPlot,
                             Corr_Min = min(cormat),
                             Corr_Max = max(cormat))
all_tissue_df <- rbind(all_tissue_df, solo_tissue_df) %>% unique()

i = i + 1
}}}}}

################################################################################
############### Combine iterative lists to dataframes ##########################
################################################################################

# Create the heatmap df (in long/tidy format)
SxW_long <- data.frame('r' = unlist(SxW_hm_list_r, use.names = F),
                        'sample' = unlist(SxW_hm_list_s, use.names = F),
                        'workflow' = unlist(SxW_hm_list_w, use.names = F))

# Convert the data to wide format
SxW_wide_y <- spread(SxW_long, sample, r)
SxW_wide_x <- spread(SxW_long, workflow, r)

# Create a dendrogram from the wide format
SxW_wide_y_mat <- as.matrix(SxW_wide_y %>% select(-workflow))
SxW_wide_x_mat <- as.matrix(SxW_wide_x %>% select(-sample))
rownames(SxW_wide_y_mat) <- SxW_wide_y$workflow
rownames(SxW_wide_x_mat) <- SxW_wide_x$sample
SxW_wide_y_dendro <- as.dendrogram(hclust(d = dist(x = SxW_wide_y_mat)))
SxW_wide_x_dendro <- as.dendrogram(hclust(d = dist(x = SxW_wide_x_mat)))
# Create dendrogram
SxW_wide_y_dendroplot <- ggdendrogram(data = SxW_wide_y_dendro, rotate = TRUE)
SxW_wide_x_dendroplot <- ggdendrogram(data = SxW_wide_x_dendro)

# Extract the order of the tips from the dendrogram
SxW_order_y <- order.dendrogram(SxW_wide_y_dendro)
SxW_order_x <- order.dendrogram(SxW_wide_x_dendro)
# Order the heatmap axes by dendrogram tips
SxW_long$workflow <- factor(x = SxW_long$workflow,
                               levels = SxW_wide_y$workflow[SxW_order_y], 
                               ordered = TRUE)
SxW_long$sample <- factor(x = SxW_long$sample,
                            levels = SxW_wide_x$sample[SxW_order_x], 
                            ordered = TRUE)

################################################################################
########################## Tissue Specific Plots ###############################
################################################################################
# Elbow plots        
elbow_out <- marrangeGrob(elbow_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-elbow-plots.pdf'), elbow_out,  width=11, height=8.5)

# Elbow plots final     
elbow_final_out <- marrangeGrob(elbow_final_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-elbow-final-plots.pdf'), elbow_final_out,  width=11, height=8.5)

# Heatmap plots
heatmap_out <- marrangeGrob(heatmap_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-heatmap-plots.pdf'), heatmap_out,  width=11, height=8.5)

# PCA Plots (NxN)
NNpca_out <- marrangeGrob(NNpca_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-NxN-pca-plots.pdf'), NNpca_out,  width=11, height=8.5)

# PCA Plots (NxP) with outlier labeled
NPpca_out <- marrangeGrob(NPpca_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-NxP-pca-outliers-labeled-plots.pdf'), NPpca_out,  width=11, height=8.5)
        
# PCA Plots NxP with outlier labeled (batch-annotated)
NPpca_batch_out <- marrangeGrob(NPpca_batch_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-NxP-pca-batch-outliers-labeled-plots.pdf'), NPpca_batch_out,  width=11, height=8.5)

# PCA Plots (NxP) with outlier removed
NPpca_noo_out <- marrangeGrob(NPpca_noo_list, nrow=2, ncol=3, layout_matrix = matrix(1:6, 2, 3, TRUE))
ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-NxP-pca-outliers-removed-plots.pdf'), NPpca_noo_out,  width=11, height=8.5)

# Heatmap of sample by workflow
# Set the breaks & colors
x <- median(SxW_long$r) - (3*sd(SxW_long$r))
brks <- round(seq(x, 1, length.out=5), digits = 3)
col <- colorRampPalette(c("blue", "white", "red"))((1-x)*1000)
SxW_hm_p <- ggplot(data = SxW_long, aes(x = sample, y = workflow, fill = r)) +
        geom_tile(color = "gray") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                         size = 6, hjust = 1)) +
        theme(axis.text.y = element_text(size = 6)) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        #guides(fill = guide_legend(override.aes = list(size = 1))) +
        ggtitle(TISSUE) +
        theme(legend.title = element_text(size = 8),
              legend.text = element_text(size = 8)) +
        theme(legend.key.size = unit(0.5, 'lines')) +
        theme(legend.key.width = unit(0.5, 'lines')) +
        scale_fill_gradientn(colors = col,name="Correlation",limits=c(x,1),breaks = brks) +
        coord_fixed() +
        theme(aspect.ratio=1)

ggsave(paste0(WD,'/plots/20200823_rnaseq-pass1a-',TIS,'-',PrintPlot,'-heatmap-workflow-sample.pdf'), SxW_hm_p,  width=11, height=8.5)

}
}
################################################################################
############# Final Report With All Tissues and Conditions #####################
################################################################################
# , times = 3)
# Create the corr_min df
# corr_min_tissue_df <- all_tissue_df %>%
#         select(Tissue,Print_Plot,Corr_Min) %>% unique() %>%
#         arrange(Tissue,Print_Plot,Corr_Min)

# Save the corr_min df
# corr_min_tissue_file <- paste0(WD,'/docs/20200823_rnaseq-pass1a-outlier-corr-mins_steep.txt')
# write.table(corr_min_tissue_df, file = corr_min_tissue_file, row.names = F, sep = '\t', quote = F)

# Create the final stats file
all_tissue_df <- all_tissue_df %>%
        arrange(Tissue, Transformation,Filter_Genes,Correlation_Metric,Order)

# Save the final table
all_tissue_file <- paste0(WD,'/docs/20200823_rnaseq-pass1a-outlier-decisions-stats_steep.txt')
write.table(all_tissue_df, file = all_tissue_file, row.names = F, sep = '\t', quote = F)



