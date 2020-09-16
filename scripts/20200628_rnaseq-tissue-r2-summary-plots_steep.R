'---
#' title: "PASS1A Rat Tissue: -- Summary R2 Plots"
#' author: "Alec Steep & Jiayu Zhang" 
#' date: "20200628"
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
#' * 2D Heatmap for all R2 plots
#' * Apply Manova to all genes all tissues
#' * Apply SOM clustering to manova tissues
#' * Bin genes with manual binning
#' * Histogram of R2 differential
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
#BiocManager::install("ggplot2")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR","gplots","tiff","splines","gam",'ggplot2')
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select
slice <- dplyr::slice
counts <- DESeq2::counts
map <- purrr::map

source(paste0(WD,'/functions/not_in.R'))
source(paste0(WD,'/functions/mouse2rat_ortho.R'))
source(paste0(WD,'/functions/sin.R'))
source(paste0(WD,'/functions/cos.R'))
source(paste0(WD,'/functions/circleFun.R'))
source(paste0(WD,'/functions/rat_mouse_ortho.R'))

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

TISSUE <- "Kidney"
# for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
        print(TISSUE)
        TISSUE1 <- TISSUE
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
        # TODO: Needs updating
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
        
        # # To load the annotation
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
        
        #' ## Load Annotations for PVE Plots
        
        #+ Load Annotations for PVE Plots
        ################################################################################
        #####     Load Annotations for PVE Plots       ##################################
        ################################################################################
        # Load the PVE values from modeling
        models_file <- paste0(WD,'/data/test/20200603_rnaseq-tissue-models-pve-table_steep.txt')
        models_df <- read.table(file = models_file ,sep = '\t', header = T, check.names = F) %>%
                as_tibble()
        
        # Load the DE values from different DE comparisons
        de_file <- paste0(WD,'/data/20200624_rnaseq-tissue-all-de-table_steep.txt')
        de_df <- read.table(file = de_file ,sep = '\t', header = T, check.names = F) %>%
                as_tibble()
        
        # Load the KEGG circadian genes TODO: Come back and calucalte orthologs, for now, manually done
        kegg_file <- paste0(WD,'/data/KEGG-04710-mouse-genes_steep.csv')
        kegg_df <- read.table(file = kegg_file ,sep = ',', header = F) %>%
          as_tibble()
        names(kegg_df) <- c('ENSEMBL_MOUSE', 'SYMBOL_MOUSE')
        # Skp1a is mouse ortholog to Skp1 in rat 1:1 high confidence
        gs <- models_df$SYMBOL_RAT %>% unique() %>% unlist() %>% as.character()
        kegg_df$SYMBOL_MOUSE[kegg_df$SYMBOL_MOUSE %!in% gs]
        kegg_df <- kegg_df %>%
          mutate(SYMBOL_RAT = ifelse(SYMBOL_MOUSE == 'Skp1a', 'Skp1', SYMBOL_MOUSE))
        keggs <- kegg_df$SYMBOL_RAT %>% unlist() %>% as.character()
        
        # Load the Manova file
        # manova_file <- paste0(WD,'/data/20200413_rnaseq-tissue-manova-table_steep.txt')
        # manova_df <- read.table(file = manova_file ,sep = '\t', header = T, check.names = F) %>%
        #         as_tibble()
        
        # Load the SOM and K-means cluster file
        clusters_file <- paste0(WD,'/data/20200413_rnaseq-tissue-cluster-table_steep.txt')
        clusters_df <- read.table(file = clusters_file ,sep = '\t', header = T, check.names = F) %>%
                as_tibble()
        
        # Join the tables
        # bin_df <- full_join(models_df, manova_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
        #         full_join(clusters_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
        #         mutate(BIN_TYPE = case_when(SIN1_R2 > GAM1_R2 + 0.15 ~ 'Circadian',
        #                                     GAM1_R2 > SIN1_R2 + 0.15 ~ 'Exercise',
        #                                     (((SIN1_R2 <= GAM1_R2 + 0.15) & (GAM1_R2 <= SIN1_R2 + 0.15)) & SIN1_R2 >= 0.3 & GAM1_R2 >= 0.3) ~ 'Ambiguous'))
        #bin_df <- left_join(models_df, manova_df, by = c("TISSUE","ENSEMBL_RAT"))
        bin_df <- full_join(models_df, de_df, by = c("TISSUE","ENSEMBL_RAT")) 
        
                # full_join(manova_df, by = c("TISSUE","ENSEMBL_RAT")) %>%
                # full_join(clusters_df, by = c("TISSUE","ENSEMBL_RAT"))
        
        # Adjust the bintype
        ######################
        # Calculate the angle from x y coordinates
        bin_df <- bin_df %>%
                mutate(c_ce = sqrt( (pve_ce2_e^2) + (pve_ce2_c^2) )) %>%
                mutate(rad = atan(pve_ce2_c/pve_ce2_e)) %>%
                mutate(GROUP_ce = case_when(
                       c_ce <= 0.2 ~ '0',
                       (c_ce > 0.2 & rad <= (1/2*pi) & rad > atan(1/0.5)) ~ '1',
                       (c_ce > 0.2 & rad <= atan(1/0.5) & rad > (1/4*pi)) ~ '2',
                       (c_ce > 0.2 & rad <= (1/4*pi) & rad > atan(0.5/1)) ~ '3',
                       (c_ce > 0.2 & rad <= atan(0.5/1) & rad > 0) ~ '4')) %>%
                mutate(GROUP_ce_n = as.numeric(GROUP_ce))
      
        # bin_df %>%
        #   filter(SYMBOL_RAT == 'Egr1') %>%
        #   select(TISSUE, SYMBOL_RAT,GROUP_ce,pve_e2, pve_c, pve_ce2_c, pve_ce2_e)
        
        
        # Collect the top 500 DE genes for each comparison 
        for(G1G2 in c('C0C7','C7E7','C0E0','E7E0')){
          lfc <- paste0(G1G2, '_log2FoldChange')
          padj <- paste0(G1G2, '_padj')
          sig <- paste0(G1G2, '_SIG')
          bin_df <- bin_df %>%
            mutate(!!sym(sig) := ifelse(((!!sym(lfc) >= 0.25 | !!sym(lfc) <= -0.25) 
                                         & !!sym(padj) <= 0.05),'SIG','NON-SIG'))
        }
# Collect top 500 E7E0 genes in Lung
bin_df$C0C7_TOP500 <- '0'
bin_df$C7E7_TOP500 <- '0'
bin_df$C0E0_TOP500 <- '0'
bin_df$E7E0_TOP500 <- '0'
names(bin_df)

for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
    (TISSUE1 <- TISSUE)
  print(TISSUE)
    for(G1G2 in c('C0C7','C7E7','C0E0')){
      print(G1G2)
      lfc <- paste0(G1G2, '_log2FoldChange')
      padj <- paste0(G1G2, '_padj')
      sig <- paste0(G1G2, '_SIG')
      t500 <- paste0(G1G2, '_TOP500')
      top500 <- c()
      top500 <- bin_df %>%
        filter(TISSUE == TISSUE1) %>%
        filter(!!sym(padj) <= 0.05) %>%
        filter(!!sym(lfc) >= 0.25) %>%
        dplyr::arrange(desc(abs(!!sym(lfc)))) %>%
        head(n =500) %>%
          select(ENSEMBL_RAT) %>% unlist() %>% as.character()
      bin_df <- bin_df %>%
        mutate(!!sym(t500) := ifelse((TISSUE == TISSUE1 & ENSEMBL_RAT %in% top500),'1',!!sym(t500)))
}  
}

# Add the KEGG annotation
bin_df <- bin_df %>% 
  mutate(KEGG_04710 = ifelse(SYMBOL_RAT %in% keggs, '1', '0'))

#' ## Tallied Genes by Bin Type

#+ Tallied Genes by Bin Type
################################################################################
#####     Tallied Genes by Bin Type       ######################################
################################################################################

# Gene Groups
########################################################
# Tally Genes by Bin Type
bin_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce) %>% 
  mutate(GROUP_N = n()) %>%
  mutate(GROUP_GENE_N = n_distinct(ENSEMBL_RAT))

tissue_counts_df <- bin_df %>%
  group_by(TISSUE) %>% 
  mutate(TISSUE_N = n()) %>%
  mutate(TISSUE_GENE_N = n_distinct(ENSEMBL_RAT)) %>% 
  select(TISSUE, TISSUE_N, TISSUE_GENE_N) %>%
  unique()
bin_counts_df <- left_join(bin_counts_df, tissue_counts_df, by = "TISSUE")

# Tally the genes and proportions
bin_tally_df <- bin_counts_df %>%
  select(TISSUE, GROUP_ce, GROUP_GENE_N, TISSUE_GENE_N) %>%
  group_by(TISSUE) %>% 
  mutate(TISSUE_N = n()) %>%
  mutate(GROUP_TYPE_FREQ = round(GROUP_GENE_N/TISSUE_GENE_N,3)) %>%
  arrange(TISSUE, GROUP_ce) %>%
  select(-TISSUE_N) %>%
  unique()

bin_tally_df %>% filter(TISSUE == TISSUE1)

View(bin_tally_df %>% arrange(GROUP_ce))

# Kegg Genes
########################################################
# Tallt genes by Tissue, gene group, and kegg status
kegg_counts_df <- bin_df %>%
    group_by(TISSUE, GROUP_ce, KEGG_04710) %>% 
    mutate(TIS_GG_KEGG_N = n()) %>%
    mutate(TIS_GG_GENE_N = n_distinct(ENSEMBL_RAT)) %>%
    select(TIS_GG_KEGG_N, TIS_GG_GENE_N) %>%
    unique() %>% 
    ungroup()

gg_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce) %>% 
  mutate(GROUP_N = n()) %>%
  mutate(GROUP_GENE_N = n_distinct(ENSEMBL_RAT)) %>%
  select(GROUP_N) %>%
  unique() %>%
  ungroup()

kegg_counts_df <- left_join(gg_counts_df, kegg_counts_df, by = c('TISSUE', 'GROUP_ce')) %>%
  filter(KEGG_04710 == '1') %>%
  select(-TIS_GG_GENE_N) %>%
  mutate(TIS_GG_KEGG_FREQ = TIS_GG_KEGG_N/GROUP_N)

# Save file and load into excel to visualize
kegg_counts_file <- paste0(WD,'/data/20200628_kegg-counts-gene-groups_steep.txt')
write.table(kegg_counts_df,kegg_counts_file,
            sep = '\t', quote = FALSE, row.names = FALSE)

# DE Genes (C0C7) 
########################################################
# Tally genes by Tissue, gene group, and kegg status
C0C7_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce, C0C7_TOP500) %>% 
  mutate(TIS_GG_C0C7_N = n()) %>%
  select(TIS_GG_C0C7_N) %>%
  unique() %>% 
  ungroup()

C0C7_counts_df <- left_join(gg_counts_df, C0C7_counts_df, by = c('TISSUE', 'GROUP_ce')) %>%
  filter(C0C7_TOP500 == '1') %>%
  mutate(TIS_GG_C0C7_FREQ = TIS_GG_C0C7_N/GROUP_N) %>%
  select(TISSUE, GROUP_ce, TIS_GG_C0C7_FREQ, GROUP_N, C0C7_TOP500, TIS_GG_C0C7_N)

# Save file and load into excel to visualize
C0C7_counts_file <- paste0(WD,'/data/20200628_C0C7-counts-gene-groups_steep.txt')
write.table(C0C7_counts_df,C0C7_counts_file,
            sep = '\t', quote = FALSE, row.names = FALSE)

# DE Genes (C7E7)
########################################################
# Tally genes by Tissue, gene group, and kegg status
C7E7_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce, C7E7_TOP500) %>% 
  mutate(TIS_GG_C7E7_N = n()) %>%
  select(TIS_GG_C7E7_N) %>%
  unique() %>% 
  ungroup()

C7E7_counts_df <- left_join(gg_counts_df, C7E7_counts_df, by = c('TISSUE', 'GROUP_ce')) %>%
  filter(C7E7_TOP500 == '1') %>%
  mutate(TIS_GG_C7E7_FREQ = TIS_GG_C7E7_N/GROUP_N)

# Save file and load into excel to visualize
C7E7_counts_file <- paste0(WD,'/data/20200628_C7E7-counts-gene-groups_steep.txt')
write.table(C7E7_counts_df,C7E7_counts_file,
            sep = '\t', quote = FALSE, row.names = FALSE)

# DE Genes (C0E0) 
#'E7E0')
########################################################
# Tally genes by Tissue, gene group, and kegg status
C0E0_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce, C0E0_TOP500) %>% 
  mutate(TIS_GG_C0E0_N = n()) %>%
  select(TIS_GG_C0E0_N) %>%
  unique() %>% 
  ungroup()

C0E0_counts_df <- left_join(gg_counts_df, C0E0_counts_df, by = c('TISSUE', 'GROUP_ce')) %>%
  filter(C0E0_TOP500 == '1') %>%
  mutate(TIS_GG_C0E0_FREQ = TIS_GG_C0E0_N/GROUP_N)

# Save file and load into excel to visualize
C0E0_counts_file <- paste0(WD,'/data/20200628_C0E0-counts-gene-groups_steep.txt')
write.table(C0E0_counts_df,C0E0_counts_file,
            sep = '\t', quote = FALSE, row.names = FALSE)

# DE Genes (E7E0) 
########################################################
# Tally genes by Tissue, gene group, and kegg status
E7E0_counts_df <- bin_df %>%
  group_by(TISSUE, GROUP_ce, E7E0_TOP500) %>% 
  mutate(TIS_GG_E7E0_N = n()) %>%
  select(TIS_GG_E7E0_N) %>%
  unique() %>% 
  ungroup()

E7E0_counts_df <- left_join(gg_counts_df, E7E0_counts_df, by = c('TISSUE', 'GROUP_ce')) %>%
  filter(E7E0_TOP500 == '1') %>%
  mutate(TIS_GG_E7E0_FREQ = TIS_GG_E7E0_N/GROUP_N)

# Save file and load into excel to visualize
E7E0_counts_file <- paste0(WD,'/data/20200628_E7E0-counts-gene-groups_steep.txt')
write.table(E7E0_counts_df,E7E0_counts_file,
            sep = '\t', quote = FALSE, row.names = FALSE)



        #' ## Save the PVE Summary File
        
        #+ Save the PVE Summary File
        ################################################################################
        #####     Save the PVE Summary File       #######################################
        ################################################################################
        # Label the gene symbols
        bin_out <- bin_df %>%
                select(TISSUE, ENSEMBL_RAT, SIN1_R2, GAM1_R2, 
                       MANOVA_PVAL, CLUSTER_ANN, BIN_TYPE) %>%
                unique()
        bin_out$ENSEMBL_RAT <- as.character(bin_out$ENSEMBL_RAT)
        bin_out <- bin_out %>%
                mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))
        bin_out <- bin_out %>%
                select(TISSUE, ENSEMBL_RAT,SYMBOL_RAT, SIN1_R2, GAM1_R2, 
                       MANOVA_PVAL, CLUSTER_ANN, BIN_TYPE)
        
        r2_file <- 
                paste0(WD,'/data/20200628_rnaseq-tissue-stats_steep.txt')
        write.table(bin_out, file=r2_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        
        #' ## Visualize PVE Plots: Bin Type
        
        #+ Visualize PVE Plots: Bin Type
        ################################################################################
        #####     Visualize PVE Plots: Bin Type       ###################################
        ################################################################################
        # Visualize Bin Type
        TISSUE1 <- 'Kidney'
        # Visualize the clusters (Tissue-specific)
        d <- bin_df %>%
                unique()
        
        p <- ggplot() +
                geom_point(data = d, aes(x = pve_ce2_e, y=pve_ce2_c, 
                                         color = GROUP_ce_n) ,alpha = 0.3) +
                xlim(0,1) + ylim(0,1) +
                # geom_abline(intercept = 0, slope = 1) +
                # geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                # geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = 0, slope = 2, linetype = 'dashed') +
                geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = 0, slope = 0.5, linetype = 'dashed') +
                geom_path(data = circleFun(c(0,0),0.40,npoints = 1000), aes(x,y)) +
                xlab("PVE Exercise") +
                ylab("PVE Circadian") +
                ggtitle("PVE Comparisons within Joint Model (Exercise + Circadian)") +
                # theme(strip.text = element_text(size=18),
                #       axis.text.x = element_text(size = 16),
                #       axis.text.y = element_text(size = 14),
                #       axis.title.x = element_text(size = 20),
                #       axis.title.y = element_text(size = 20),
                #       legend.text = element_text(size = 16)) +
                facet_wrap(vars(TISSUE)) +
                # scale_color_manual(values = c('0' = 'grey',
                #                               '1' = 'orange',
                #                               '2' = 'gold',
                #                               '3' = 'steelblue1',
                #                               '4' = 'steelblue4')) +
                scale_color_continuous(high = "#56B1F7", low = "#132B43",
                                       name = "Group") +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-',TISSUE1,'-theoretical-bins-anova2_steep.pdf'),
            width=26,height=14)
        plot(p)
        dev.off()
        ec_colors
        # Visualize the clusters across tissues
        d <- bin_df %>%
                filter(TISSUE == TISSUE1)
        #filter(MANOVA_PVAL <= 0.05)
        
        # PVE (between models)
        p <- ggplot() +
                geom_point(data = bin_df, 
                           aes(x = pve_ce2_e, y= pve_ce2_c, color = BIN_TYPE_ce), 
                           alpha = 0.8) +
                xlim(0,1) + 
                ylim(0,1) +
                scale_x_continuous() +
                geom_abline(intercept = 0, slope = 1) +
                geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("PVE Polynomial Model (Exercise)") +
                ylab("PVE SIN/COS Model (Circadian)") +
                ggtitle("Circadian vs. Exercise:\nPVE Comparisons Between Models") +
                theme_bw() +
                theme(strip.text = element_text(size=18),
                      axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.text = element_text(size = 16)) +
                facet_wrap(vars(TISSUE)) +
                scale_color_manual(values = c('Ambiguous_High' = "red",
                                              'Ambiguous_Low' = "grey",
                                              'Circadian' = "green",
                                              'Exercise' = 'blue')) +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-theoretical-bins_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        # PVE (within models -- circadian first)
        p <- ggplot() +
                geom_point(data = d, 
                           aes(x = pve_ce_e, y= pve_ce_c, color = BIN_TYPE), 
                           alpha = 0.8) +
                # geom_density_2d_filled(data = d, 
                #                        aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.5) +
                # geom_density_2d(data = d, 
                #                 aes(x = GAM1_R2, y= SIN1_R2), 
                #                 size = 0.25, colour = "black") +
                xlim(0,1) + ylim(0,1) +
                geom_abline(intercept = 0, slope = 1) +
                geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("PVE Exercise") +
                ylab("PVE Circadian") +
                ggtitle("Circadian vs. Exercise:\nPVE Comparisons within Model") +
                theme_bw() +
                theme(strip.text = element_text(size=18),
                      axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.text = element_text(size = 16)) +
                facet_wrap(vars(TISSUE)) +
                scale_color_manual(values = c('Ambiguous_High' = "red",
                                              'Ambiguous_Low' = "grey",
                                              'Circadian' = "green",
                                              'Exercise' = 'blue')) +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-theoretical-bins-ce_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        # PVE (within models -- exercise first)
        p <- ggplot() +
                geom_point(data = d, 
                           aes(x = pve_ec_e, y= pve_ec_c, color = BIN_TYPE), 
                           alpha = 0.8) +
                # geom_density_2d_filled(data = d, 
                #                        aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.5) +
                # geom_density_2d(data = d, 
                #                 aes(x = GAM1_R2, y= SIN1_R2), 
                #                 size = 0.25, colour = "black") +
                xlim(0,1) + ylim(0,1) +
                geom_abline(intercept = 0, slope = 1) +
                geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("PVE Exercise") +
                ylab("PVE Circadian") +
                ggtitle("Circadian vs. Exercise:\nPVE Comparisons within Model") +
                theme_bw() +
                theme(strip.text = element_text(size=18),
                      axis.text.x = element_text(size = 16),
                      axis.text.y = element_text(size = 14),
                      axis.title.x = element_text(size = 20),
                      axis.title.y = element_text(size = 20),
                      legend.text = element_text(size = 16)) +
                facet_wrap(vars(TISSUE)) +
                scale_color_manual(values = c('Ambiguous_High' = "red",
                                              'Ambiguous_Low' = "grey",
                                              'Circadian' = "green",
                                              'Exercise' = 'blue')) +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-theoretical-bins-ec_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        #' ## Visualize PVE Plots: Plane Jane Heroine
        
        #+ Visualize PVE Plots: Plane Jane Heroine
        ################################################################################
        #####     Visualize PVE Plots: Plane Jane Heroine       ##############################
        ################################################################################
        dim(bin_df)
        # Single Tissue
        p <- bin_df %>%
                filter(TISSUE == TISSUE1) %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                #sample_n(2000) %>%
                ggplot(aes(x = pve_e2, y= pve_c, alpha = 0.005)) +
                geom_point() +
                #scale_fill_gradientn(colours=r, trans = "log") +
                xlim(0,1) + ylim(0,1) +
                geom_abline(intercept = 0, slope = 1) +
                # geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                # geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("PVE Polynomial Model (Exercise)") +
                ylab("PVE SIN/COS Model (Circadian)") +
                ggtitle(paste0("PVE Comparisons Between Models")) +
                #labs(fill = 'Log(Gene Counts)') +
                theme(legend.position = "none") +
                coord_equal() +
                facet_wrap( ~ TISSUE)
        pdf(paste0(WD,'/plots/20200628_rnaseq-',TISSUE1,'-plane-jane-heroine_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        #' ## Visualize R2 Plots: 2D Histograms
        
        #+ Visualize R2 Plots: 2D Histograms
        ################################################################################
        #####     Visualize R2 Plots: 2D Histograms       ##############################
        ################################################################################
        
        # Color Housekeeping
        rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
        r <- rf(32)
        
        # Visualize a 2D Heatmap for every tissue
        # Single Tissue
        p <- bin_df %>%
                filter(TISSUE == TISSUE1) %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                #sample_n(200) %>%
                ggplot(aes(x = pve_ce2_e, y= pve_ce2_c)) +
                stat_bin2d(bins = 50) +
                scale_fill_gradientn(colours=r, trans = "log") +
                xlim(0,1) + ylim(0,1) +
          geom_abline(intercept = 0, slope = 2, linetype = 'dashed') +
          geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
          geom_abline(intercept = 0, slope = 0.5, linetype = 'dashed') +
          geom_path(data = circleFun(c(0,0),0.40,npoints = 1000), aes(x,y)) +
                xlab("PVE Exercise") +
                ylab("PVE Circadian") +
                ggtitle(paste0("PVE Comparisons within Combined Model")) +
                labs(fill = 'Log(Gene Counts)') +
                coord_equal() +
                facet_wrap( ~ TISSUE)
        pdf(paste0(WD,'/plots/20200628_rnaseq-',TISSUE1,'-faceted-2Dhistogram_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        
        
        # All Tissues
        p <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                #sample_n(2000) %>%
                ggplot(aes(x = pve_ce2_e, y= pve_ce2_c)) +
                stat_bin2d(bins = 100) +
                scale_fill_gradientn(colours=r, trans = "log") +
                xlim(0,1) + ylim(0,1) +
          geom_abline(intercept = 0, slope = 2, linetype = 'dashed') +
          geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
          geom_abline(intercept = 0, slope = 0.5, linetype = 'dashed') +
          geom_path(data = circleFun(c(0,0),0.40,npoints = 1000), aes(x,y)) +
          xlab("PVE Exercise") +
          ylab("PVE Circadian") +
          ggtitle("PVE Comparisons within Joint Model (Exercise + Circadian)") +
                theme(plot.title = element_text(size = 20),
                      strip.text = element_text(size=12),
                       axis.text.x = element_text(size = 10),
                       axis.text.y = element_text(size = 10),
                       axis.title.x = element_text(size = 18),
                       axis.title.y = element_text(size = 18),
                       legend.text = element_text(size = 10),
                      legend.title = element_text(size = 18)) +
                labs(fill = 'Log(Gene Counts)') +
                coord_equal() +
                facet_wrap( ~ TISSUE)
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-2Dhistogram_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        #' ## Visualize R2 Plots: K-means Clusters
        
        #+ Visualize R2 Plots: K-means Clusters
        ################################################################################
        #####     Visualize R2 Plots: K-means Clusters       ###########################
        ################################################################################
        
        for(cluster in c('Early_Up', 'Middle_Up','Late_Up','Delayed_Up',
                           'Early_Down', 'Middle_Down','Late_Down','Delayed_Down')){
                print(cluster)
                # Visualize the clusters across tissues
                d <- bin_df %>%
                        #filter(TISSUE == TISSUE1) %>%
                        filter(CLUSTER_ANN == cluster) %>%
                        filter(MANOVA_PVAL <= 0.05)
                p <- ggplot() +
                        #geom_point(data = d, aes(x = GAM1_R2, y= SIN1_R2), alpha= 0.01, size = 0.5) +
                        #geom_density_2d_filled(contour_var = "count") +
                        geom_density_2d_filled(data = d, aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.5) +
                        geom_density_2d(data = d, aes(x = GAM1_R2, y= SIN1_R2), size = 0.25, colour = "black") +
                        #geom_density_2d() +
                        # geom_point(data = filter(bin_df, CLUSTER_ANN == 'Early_Up'),
                        #            aes(x = GAM1_R2, y= SIN1_R2), color = 'red') +
                        # stat_bin2d(data = filter(bin_df, CLUSTER_ANN == 'Early_Up',
                        #                          aes(x = GAM1_R2, y= SIN1_R2),
                        #                          bins = 100)
                        xlim(0,1) + ylim(0,1) +
                        geom_abline(intercept = 0, slope = 1) +
                        geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                        geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                        xlab("R^2 Natural Spline Model (Exercise)") +
                        ylab("R^2 SIN/COS Model (Circadian)") +
                        ggtitle(paste0(TISSUE1,":\nR2 Comparisons Between Models")) +
                        facet_wrap(vars(TISSUE)) +
                        coord_equal()
                pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-',cluster,'-annotated_steep.pdf'),width=26,height=14)
                plot(p)
                dev.off()
        }
        
        #' ## Visualize R2 Plots: DE Genes (Different Comparisons)
        
        #+ Visualize R2 Plots: DE Genes (Different Comparisons)
        ################################################################################
        #####     Visualize R2 Plots: DE Genes (Different Comparisons)       ###########
        ################################################################################
        
        # For each Tissue, generate a bar plot of DE gene frequencies
        bin_long <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                select(C0C7_SIG,C0E0_SIG,C0E0.5_SIG,C0E1_SIG,C0E4_SIG,C0E7_SIG,
                       C0E24_SIG,C0E48_SIG,C7E0_SIG,C7E0.5_SIG,C7E1_SIG,C7E4_SIG,C7E7_SIG,
                       C7E24_SIG,C7E48_SIG,BIN_TYPE,TISSUE,ENSEMBL_RAT) %>%
                unique() %>%
                unite(TISSUE_ENSEMBL_RAT_BIN_TYPE, c("TISSUE","ENSEMBL_RAT","BIN_TYPE")) %>%
                pivot_longer(-TISSUE_ENSEMBL_RAT_BIN_TYPE, 
                             names_to = 'COMP', values_to = 'SIG') %>%
                filter(SIG == 'SIG') %>% unique()
        bin_long <- bin_long %>% 
                separate(TISSUE_ENSEMBL_RAT_BIN_TYPE, 
                                 c('TISSUE','ENSEMBL_RAT','BIN_TYPE'), '_')
        # Adjust the order of factors
        bin_long$COMP <- factor(bin_long$COMP, 
                                levels = c('C0E0_SIG','C7E0_SIG',
                                           'C0E0.5_SIG','C7E0.5_SIG',
                                           'C0E1_SIG','C7E1_SIG',
                                           'C0E4_SIG','C7E4_SIG',
                                           'C0E7_SIG','C7E7_SIG',
                                           'C0E24_SIG','C7E24_SIG',
                                           'C0E48_SIG','C7E48_SIG','C0C7_SIG'))
        # CvsE
        for(TISSUE in c('Lung','Hypothalamus','Aorta','Liver', 'Kidney', 'Adrenal', 'Brown Adipose', 'Cortex','Gastrocnemius', 'Heart', 'Hippocampus','Ovaries','Spleen','Testes', 'White Adipose')){
                print(TISSUE)
                TISSUE1 <- TISSUE
                p <- bin_long %>%
                filter(TISSUE == TISSUE1) %>%
                ggplot(aes(x = COMP, fill = BIN_TYPE)) +
                geom_bar() +
                xlab("") +
                ylab("Gene Count") +
                ggtitle(paste0('Frequency of Differentially Expressed Genes (',TISSUE1,')'))
                pdf(paste0(WD,'/plots/20200628_rnaseq-',str_replace(TISSUE1,' ','_'),'-gene-groups_steep.pdf'),width=14,height=8)
                plot(p)
                dev.off()
        }
        
        # Iterate through different comparisons and visualize
        comp <- 'C0C7'
        for(comp in c('C0C7','C0E7','C7E7')){
                print(comp)
                # Visualize the clusters across tissues
                d1 <- bin_df %>%
                        #filter(TISSUE == TISSUE1) %>%
                        filter(MANOVA_PVAL <= 0.05)
                        #sample_n(20000)
                d2 <- bin_df %>%
                        #filter(TISSUE == TISSUE1) %>%
                        filter(MANOVA_PVAL <= 0.05) %>%
                        filter(!!as.name(paste0(comp,'_SIG')) == 'SIG')
                        #sample_n(20000)
                p <- ggplot() +
                        geom_point(data = d1, 
                                   aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.001) +
                        geom_point(data = d2, 
                                   aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.005, color = 'red') +
                        #geom_density_2d_filled(data = d, 
                        #                       aes(x = GAM1_R2, y= SIN1_R2), alpha = 0.5) +
                        # geom_density_2d(data = d, 
                        #                 aes(x = GAM1_R2, y= SIN1_R2), 
                        #                 size = 0.25, colour = "black") +
                        xlim(0,1) + ylim(0,1) +
                        geom_abline(intercept = 0, slope = 1) +
                        geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                        geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                        xlab("R^2 Natural Spline Model (Exercise)") +
                        ylab("R^2 SIN/COS Model (Circadian)") +
                        ggtitle(paste0("R2 Comparisons Between Models:
Differentially Expressed Genes (",comp,")")) +
                        facet_wrap(vars(TISSUE)) +
                        coord_equal()
                pdf(paste0(WD,'/plots/20200628_rnaseq-',TISSUE1,'-faceted-',comp,'-annotated_steep.pdf'),width=26,height=14)
                plot(p)
                dev.off()
        }

        
        
        #' ## Tallied Genes by Bin Type
        
        #+ Tallied Genes by Bin Type
        ################################################################################
        #####     Tallied Genes by Bin Type       ######################################
        ################################################################################
        
        bin_df %>%
                filter(TISSUE == TISSUE1) %>%
                select(GROUP_ce) %>%
                table()
        
        
        # Tally Genes by Bin Type
        bin_counts_df <- bin_df %>%
                group_by(TISSUE, GROUP_ce) %>% 
                mutate(GROUP_N = n()) %>%
                mutate(GROUP_GENE_N = n_distinct(ENSEMBL_RAT))
        # bin_counts_df %>% select(TISSUE, GROUP_ce, GROUP_N) %>% 
        #         filter(TISSUE == 'Kidney') %>%
        #         unique()
        
        tissue_counts_df <- bin_df %>%
                group_by(TISSUE) %>% 
                mutate(TISSUE_N = n()) %>%
                mutate(TISSUE_GENE_N = n_distinct(ENSEMBL_RAT)) %>% 
                select(TISSUE, TISSUE_N, TISSUE_GENE_N) %>%
                unique()
        bin_counts_df <- left_join(bin_counts_df, tissue_counts_df, by = "TISSUE")
        
        # Tally the genes and proportions
        bin_tally_df <- bin_counts_df %>%
                select(TISSUE, GROUP_ce, GROUP_GENE_N, TISSUE_GENE_N) %>%
                group_by(TISSUE) %>% 
                mutate(TISSUE_N = n()) %>%
                mutate(GROUP_TYPE_FREQ = round(GROUP_GENE_N/TISSUE_GENE_N,3)) %>%
                arrange(TISSUE, GROUP_ce) %>%
                select(-TISSUE_N) %>%
                unique()

        bin_tally_df %>% filter(TISSUE == TISSUE1)
        
        View(bin_tally_df %>% arrange(GROUP_ce))
        
        # Save the exercise counts
        e_tally_out <- bin_tally_df %>%
                filter(GROUP_ce == '4')
        e_tally_file <- 
                paste0(WD,'/data/20200628_rnaseq-tissue-exercise-gene-tallies_steep.txt')
        write.table(e_tally_out, file=e_tally_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        # Save the circadian counts
        c_tally_out <- bin_tally_df %>%
                filter(GROUP_ce == '1')
        c_tally_file <- 
                paste0(WD,'/data/20200628_rnaseq-tissue-circadian-gene-tallies_steep.txt')
        write.table(c_tally_out, file=c_tally_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        
        #' ## Visualize Tallied Genes
        
        #+ Visualize Tallied Genes
        ################################################################################
        #####     Visualize Tallied Genes       ########################################
        ################################################################################
        
        # Label the gene symbols
        bin_df$ENSEMBL_RAT <- as.character(bin_df$ENSEMBL_RAT)
        bin_df <- bin_df %>%
                mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))
        
        # Tally the circadian genes across tissues
        g0_tally_df <- bin_df %>%
          #filter(MANOVA_PVAL <= 0.05) %>%
          filter(GROUP_ce == '0') %>%
          filter(!is.na(SYMBOL_RAT)) %>%
          group_by(GROUP_ce,SYMBOL_RAT) %>% 
          select(TISSUE, GROUP_ce, SYMBOL_RAT) %>%
          #mutate(GENE_BIN_N = n()) %>%
          mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
          select(-TISSUE) %>%
          unique() %>%
          #filter(GENE_TISSUE_N >= 8) %>%
          arrange(desc(GENE_TISSUE_N))
        g1_tally_df <- bin_df %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                filter(GROUP_ce == '1') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(GROUP_ce,SYMBOL_RAT) %>% 
                select(TISSUE, GROUP_ce, SYMBOL_RAT, ENSEMBL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                #filter(GENE_TISSUE_N >= 8) %>%
                arrange(desc(GENE_TISSUE_N))
        g1_file <- paste0(WD,'/data/20200628_g1-ranked-list-tissue_steep.txt')
        write.table(g1_tally_df, file=g1_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        g2_tally_df <- bin_df %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                filter(GROUP_ce == '2') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(GROUP_ce,SYMBOL_RAT) %>% 
                select(TISSUE, GROUP_ce, SYMBOL_RAT, ENSEMBL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                #filter(GENE_TISSUE_N >= 8) %>%
                arrange(desc(GENE_TISSUE_N))
        g2_file <- paste0(WD,'/data/20200628_g2-ranked-list-tissue_steep.txt')
        write.table(g2_tally_df, file=g2_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        g3_tally_df <- bin_df %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                filter(GROUP_ce == '3') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(GROUP_ce,SYMBOL_RAT) %>% 
                select(TISSUE, GROUP_ce, SYMBOL_RAT, ENSEMBL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                #filter(GENE_TISSUE_N >= 8) %>%
                arrange(desc(GENE_TISSUE_N))
        g3_file <- paste0(WD,'/data/20200628_g3-ranked-list-tissue_steep.txt')
        write.table(g3_tally_df, file=g3_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        g4_tally_df <- bin_df %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                filter(GROUP_ce == '4') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(GROUP_ce,SYMBOL_RAT) %>% 
                select(TISSUE, GROUP_ce, SYMBOL_RAT, ENSEMBL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                #filter(GENE_TISSUE_N >= 8) %>%
                arrange(desc(GENE_TISSUE_N))
        g4_file <- paste0(WD,'/data/20200628_g4-ranked-list-tissue_steep.txt')
        write.table(g4_tally_df, file=g4_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        
        tally_df <- rbind(g1_tally_df, g2_tally_df, g3_tally_df, g4_tally_df) %>%
          ungroup()
        
        # Find mouse orthologs for enrichr
        orth_df <- tally_df %>%
          ungroup() %>%
          select(ENSEMBL_RAT) %>% 
          unique() %>% as.data.frame() %>%
          rat_mouse_ortho(column = 'ENSEMBL_RAT', direction = 'rat2mouse') %>%
          as_tibble()
        tally_df <- left_join(tally_df, orth_df, by='ENSEMBL_RAT')
        
        enrichr_final_df <- data.frame()
        for(gg in c('1','2','3','4')){
          # Select the Group 1 Genes
          df <- tally_df %>%
            filter(GROUP_ce == gg) %>%
            filter(GENE_TISSUE_N >= 2)
          
          # Pathway Analysis
          ###############################
          
          # Collect gene symbols
          symbols <- df %>%
            filter(!is.na(SYMBOL_MOUSE)) %>%
            select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()
          
          # Collect the desired pathways to query for enrichment
          desired <- c('KEGG_2019_Mouse')
          enriched <- enrichr(symbols, desired)
          
          n<-1
          enrichr_df <- data.frame()
          for(db in desired){
            # Collapse the data into flat format
            if(!is.null(nrow(enriched[[desired]]))){
              enriched[[n]]$Data.Base <- db
              enrichr_df <- rbind(enrichr_df,enriched[[n]])
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
              dplyr::slice(1:10) %>%
              group_by(Term) %>%
              arrange(FDR_nlog) %>%
              ungroup() %>%
              mutate(Gene_Group = gg)
            enrichr_final_df <- rbind(enrichr_final_df, enrichr_plot)
          }
        }
          # Create a special column that will allow for proper ordering of Pathways in Plot:
          # Ordering will have 2 priorities in this order: Up/Down regulation & -Log(FDR)
          # Visualize results with barplot
          #pdf(paste0(WD,'/plots/20200426_rnaseq-',TIS,'-pathways-sexmod-C0vsC7_steep.pdf'),
          #    width = 12, height = 6)
        names(enrichr_final_df)
        ggplot(enrichr_final_df, aes(x = reorder(Term, FDR_nlog), y = FDR_nlog, size = Odds.Ratio, col = Gene_Group)) +
          geom_point() +
          coord_flip() +
          ylab('-Log(FDR)') +
          xlab('') +
          theme(
                strip.background = element_blank(),
                axis.text = element_text(size=12, colour = 'black'),
                axis.ticks = element_line(colour = 'black'),
                axis.title=element_text(size=12,face='bold'),
                strip.text = element_blank())
        # scale_color_gradient(low = "red",
        #                        high = "black") 
        theme(axis.text = element_text(size = rel(1)))
          p <- ggplot(enrichr_final_df, 
                      aes(x=reorder(Term, FDR_nlog), y=FDR_nlog, fill = Gene_Group)) +
            geom_bar(stat='identity') +
            coord_flip() +
            ylab('-Log(FDR)') +
            xlab('') +
            theme_linedraw() +
            guides(fill=guide_legend(title='Gene Group')) +
            theme(panel.background = element_blank(),
                  plot.background = element_blank(),
                  strip.background = element_blank(),
                  axis.text = element_text(size=12, colour = 'black'),
                  axis.ticks = element_line(colour = 'black'),
                  axis.title=element_text(size=12,face='bold'),
                  strip.text = element_blank()) +
            theme(panel.grid.major = element_blank()) +
            ggtitle(paste0('Pathway Enrichment of Gene Groups (',desired,')'))
          pdf(paste0(WD,'/plots/20200628_rnaseq-',TISSUE1,'-',gene_group,'-pathways_steep.pdf'),
              width = 12, height = 6)
          plot(p)
          dev.off()
        }
        


enrichr_final_df$TISSUE <- as.factor(enrichr_final_df$TISSUE)

# Circadian Plot
#########################################

enrichr_circ_df <- enrichr_final_df %>%
  filter(Total.N <= 500) %>%
  filter(Total.N >= 10) %>%
  filter(Adjusted.P.value <= 0.1) %>%
  arrange(desc(FDR_nlog)) %>%
  filter(BIN_TYPE == 'Circadian') %>%
  group_by(TISSUE) %>%
  dplyr::slice(1:10)

enrichr_circ_df %>%
  filter(BIN_TYPE == 'Circadian') %>%
  ggplot(aes(x=reorder(Term, FDR_nlog), y=FDR_nlog)) +
  geom_bar(stat='identity') +
  coord_flip() +
  ylab('-Log(FDR)') +
  xlab('') +
  theme_linedraw() +
  guides(fill=guide_legend(title='Gene Expression')) +
  # theme(panel.background = element_blank(),
  #        plot.background = element_blank()) +
  # #       strip.background = element_blank(),
  # #       axis.text = element_text(size=12, colour = 'black'),
  # #       axis.ticks = element_line(colour = 'black'),
  # #       axis.title=element_text(size=12,face='bold'),
  # #       strip.text = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ggtitle(paste0('Circadian-Modeled Genes
Pathway Enrichment (',desired,')')) +
  facet_wrap(vars(TISSUE))


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        c_tally_out <- bin_df %>%
                #filter(MANOVA_PVAL <= 0.05) %>%
                filter(BIN_TYPE == 'Circadian') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(BIN_TYPE,SYMBOL_RAT) %>% 
                select(TISSUE, BIN_TYPE, SYMBOL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                filter(GENE_TISSUE_N >= 4) %>%
                arrange(desc(GENE_TISSUE_N))
        c_tally_stat <- bin_df %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                filter(BIN_TYPE == 'Circadian') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(BIN_TYPE,SYMBOL_RAT) %>% 
                select(TISSUE, BIN_TYPE, SYMBOL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                filter(GENE_TISSUE_N >= 5) %>%
                arrange(desc(GENE_TISSUE_N))
        dim(c_tally_stat)
        summary(c_tally_stat)
        
        # Save the tallies
        c_tally_file <- 
                paste0(WD,'/data/20200628_rnaseq-tissue-shared-circadian-gene-tallies_steep.txt')
        write.table(c_tally_out, file=c_tally_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        
       
        
        c_tally_stat <- bin_df %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                filter(BIN_TYPE == 'Circadian') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(BIN_TYPE,SYMBOL_RAT) %>% 
                select(TISSUE, BIN_TYPE, SYMBOL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                filter(GENE_TISSUE_N >= 1) %>%
                arrange(desc(GENE_TISSUE_N))
        e_tally_stat <- bin_df %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                filter(BIN_TYPE == 'Exercise') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(BIN_TYPE,SYMBOL_RAT) %>% 
                select(TISSUE, BIN_TYPE, SYMBOL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                filter(GENE_TISSUE_N >= 1) %>%
                arrange(desc(GENE_TISSUE_N))
        dim(e_tally_stat)
        
        # Visualize the number of genes per tissue
        rolling_tissue_df <- data.frame()
        for(n in 15:1){
                tissue_n_e <- e_tally_stat %>%
                        filter(GENE_TISSUE_N >= n) %>% nrow()
                tissue_n_c <- c_tally_stat %>%
                        filter(GENE_TISSUE_N >= n) %>% nrow()
                df_e <- data.frame('Tissues' = n, 'Genes' = tissue_n_e,'Gene_Group' = 'Group_3')
                df_c <- data.frame('Tissues' = n, 'Genes' = tissue_n_c,'Gene_Group' = 'Group_1')
                rolling_tissue_df <- rbind(rolling_tissue_df, df_e, df_c)
                rolling_tissue_df <- rolling_tissue_df %>%
                        arrange(desc(Gene_Group), desc(Tissues))
        }
        
        # Plot the gene tallies
        rolling_tissue_df %>%
                ggplot(aes(x = Tissues, y = Genes, color = Gene_Group)) +
                geom_line(size =1) +
                #geom_smooth(se = F, lm = ns(4)) +
                scale_color_manual(values = c('Group_1' = "green",
                                             'Group_3' = "blue")) +
                scale_x_reverse() +
                ylim(c(0,41)) +
                xlim(c(8,15)) +
                xlab("Number of Tissues") +
                ylab("Number of Shared Genes Across Tissues") +
                ggtitle("Number of Shared Genes Across Tissues by Gene Group")
        
        
        e_tally_out <- bin_df %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                filter(BIN_TYPE == 'Exercise') %>%
                filter(!is.na(SYMBOL_RAT)) %>%
                group_by(BIN_TYPE,SYMBOL_RAT) %>% 
                select(TISSUE, BIN_TYPE, SYMBOL_RAT) %>%
                #mutate(GENE_BIN_N = n()) %>%
                mutate(GENE_TISSUE_N = n_distinct(TISSUE)) %>%
                select(-TISSUE) %>%
                unique() %>%
                filter(GENE_TISSUE_N >= 8) %>%
                arrange(desc(GENE_TISSUE_N))
        
        # Save the tallies
        e_tally_file <- 
                paste0(WD,'/data/20200628_rnaseq-tissue-shared-exercise-gene-tallies_steep.txt')
        write.table(e_tally_out, file=e_tally_file,
                    sep = '\t', quote = FALSE, row.names = FALSE)
        
        # Visualize tallied genes
        d <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                mutate(circ_gene_label = ifelse(SYMBOL_RAT %in% c(c_tally_df$SYMBOL_RAT),
                                           SYMBOL_RAT, ''))
        d2 <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                mutate(circ_gene_label = ifelse(SYMBOL_RAT %in% c(c_tally_df$SYMBOL_RAT),
                                                SYMBOL_RAT, '')) %>%
                filter(circ_gene_label != '' & BIN_TYPE == 'Circadian') %>%
                select(TISSUE, SYMBOL_RAT, circ_gene_label, BIN_TYPE,GAM1_R2,SIN1_R2) %>%
                unique()
        dim(d2)
        
        p <- ggplot() +
                geom_point(data = d, 
                           aes(x = GAM1_R2, y= SIN1_R2, color = BIN_TYPE), 
                           alpha = 0.5) +
                geom_label_repel(data = d2,
                                 aes(x = GAM1_R2, y= SIN1_R2),
                                 label=d2$circ_gene_label, hjust=0, vjust=0) +
                xlim(0,1) + ylim(0,1) +
                geom_abline(intercept = 0, slope = 1) +
                geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("R^2 Natural Spline Model (Exercise)") +
                ylab("R^2 SIN/COS Model (Circadian)") +
                ggtitle("Circadian vs. Exercise:\nR2 Comparisons Between Models") +
                facet_wrap(vars(TISSUE)) +
                scale_fill_manual(values = c('Ambiguous' = "green",
                                             'Circadian' = "red",
                                             'Exercise' = 'blue')) +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-circ-genes_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
        
        # Visualize tallied genes
        d <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                mutate(ex_gene_label = ifelse(SYMBOL_RAT %in% c(e_tally_stat$SYMBOL_RAT),
                                                SYMBOL_RAT, ''))
        d2 <- bin_df %>%
                #filter(TISSUE == TISSUE1) %>%
                filter(MANOVA_PVAL <= 0.05) %>%
                mutate(ex_gene_label = ifelse(SYMBOL_RAT %in% c(e_tally_stat$SYMBOL_RAT),
                                                SYMBOL_RAT, '')) %>%
                filter(ex_gene_label != '' & BIN_TYPE == 'Exercise') %>%
                select(TISSUE, SYMBOL_RAT, ex_gene_label, BIN_TYPE,GAM1_R2,SIN1_R2) %>%
                unique()
        dim(d2)
        
        p <- ggplot() +
                geom_point(data = d, 
                           aes(x = GAM1_R2, y= SIN1_R2, color = BIN_TYPE), 
                           alpha = 0.5) +
                geom_label_repel(data = d2,
                                 aes(x = GAM1_R2, y= SIN1_R2),
                                 label=d2$ex_gene_label, hjust=0, vjust=0) +
                xlim(0,1) + ylim(0,1) +
                geom_abline(intercept = 0, slope = 1) +
                geom_abline(intercept = 0.15, slope = 1, linetype = 'dashed') +
                geom_abline(intercept = -0.15, slope = 1, linetype = 'dashed') +
                xlab("R^2 Natural Spline Model (Exercise)") +
                ylab("R^2 SIN/COS Model (Circadian)") +
                ggtitle("Circadian vs. Exercise:\nR2 Comparisons Between Models") +
                facet_wrap(vars(TISSUE)) +
                scale_fill_manual(values = c('Ambiguous' = "green",
                                             'Circadian' = "red",
                                             'Exercise' = 'blue')) +
                coord_equal()
        pdf(paste0(WD,'/plots/20200628_rnaseq-all-tissues-faceted-shared-exercise-genes_steep.pdf'),width=26,height=14)
        plot(p)
        dev.off()
       # 
        
        
        
        
        
        
        
        
        
        
        
        

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





################################################################################
# bin_test <- bin_df %>%
#         filter(MANOVA_PVAL <= 0.05) %>%
#         sample_n(100)
# 
# bin_test$R2_DIST <- vdist2d(a1=bin_test$GAM1_R2, a2=bin_test$GAM1_R2,
#                             b1=rep(0,nrow(bin_test)), 
#                             b2=rep(0,nrow(bin_test)),
#                             c1=rep(1,nrow(bin_test)),
#                             c2=rep(1,nrow(bin_test)))
# select(R2_DIST,BIN_TYPE)
# ggplot(aes(x = R2_DIST, fill = BIN_TYPE)) +
#         geom_freqpoly() +
#         facet_wrap(~ TISSUE)
# dist2d <- function(a1,a2,b1=0,b2=0,c1=1,c2=1) {
#         a <- c(a1, a2)
#         b <- c(b1,b2)
#         c <- c(c1,c2)
#         v1 <- b - c
#         v2 <- a - b
#         m <- cbind(v1,v2)
#         d <- abs(det(m))/sqrt(sum(v1*v1))
#         d
# } 
# vdist2d <- Vectorize(dist2d)
# 
# vdist2d(bin_df$GAM1_R2[1], bin_df$SIN1_R2[1],0,0,1,1)

# library(sf)
# library(sp)
# data(meuse)
# 
# # create new line - please note the small changes relative to your code : 
# # x = first column and cbind to have a matrix instead of a data.frame
# newline = cbind(x = seq(from = 178000, to = 181000, length.out = 100),
#                 y = seq(from = 330000, to = 334000, length.out = 100))
# 
# newline = cbind(x = seq(from = 0, to = 1, length.out = 100),
#                 y = seq(from = 0, to = 1, length.out = 100))
# 
# # transform the points and lines into spatial objects
# bin_test <- bin_df %>%
#         select(ENSEMBL_RAT, TISSUE, GAM1_R2, SIN1_R2, BIN_TYPE) %>%
#         sample_n(10000) %>%
#         unique()
# 
# bin_test <- st_as_sf(bin_test, coords = c("GAM1_R2", "SIN1_R2"))
# newline <- st_linestring(newline)
# 
# # Compute the distance - works also for non straight lines !
# bin_test$R2_DIST <- st_distance(bin_test, newline) %>% as.numeric()
# ## [1] 291.0 285.2 409.8 548.0 647.6 756.0 510.0 403.8 509.4 684.8
# 
# bin_test %>%
#         ggplot(aes(x = R2_DIST, fill = BIN_TYPE)) +
#         geom_density(alpha=0.4) +
#         xlim(0,0.5) +
#         facet_wrap(~ TISSUE)








