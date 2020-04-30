#'---
#' title: "PASS1A Rat (Female) Kidney: -- Control 0 hr vs Control 7 hr Differential Gene Expression"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200426"
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
#' * Perform a differential gene expression analysis between Control 0 and Control 7 to determine DE genes
#' * Determine if DE gene list is enriched for circadian genes in the kidney
#' * Determine if DE gene list is enriched for circadian rhythm pathways
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
#BiocManager::install("enrichR")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret","enrichR")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

################################################################################
######################### Functions ############################################
################################################################################

# Set select
select <- dplyr::select

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and Author
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

# Functions to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
rat_mouse_ortho <- function(x, column = "ENSEMBL_RAT", direction = 'rat2mouse') {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        if ( column %!in% c('ENSEMBL_RAT','ENSEMBL_MOUSE') ){
                stop("'column' must be either 'ENSEMBL_RAT' or 'ENSEMBL_MOUSE'", class.= FALSE)
        }
        if ( direction %!in% c('rat2mouse', 'mouse2rat')){
                stop("'direction' must be either 'rat2mouse' or 'mouse2rat'")
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
        # Create ortholog table
        if (direction == 'rat2mouse'){
                ortho_df <- getLDS(attributes=c("ensembl_gene_id",
                                             "mmusculus_homolog_orthology_confidence"),
                                   filters="ensembl_gene_id", 
                                   values = x[[column]], 
                                   mart=mart_rn_ens,
                                   attributesL=c("ensembl_gene_id"), 
                                   martL=mart_mm_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the HUGO symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_RAT') %>%
                        filter(!is.na(ENSEMBL_RAT)) %>%
                        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))
                
        }else{
                ortho_df <- getLDS(attributes=c("ensembl_gene_id",
                                                "rnorvegicus_homolog_orthology_confidence"),
                                   filters="ensembl_gene_id", 
                                   values = x[[column]], 
                                   mart=mart_mm_ens,
                                   attributesL=c("ensembl_gene_id"), 
                                   martL=mart_rn_ens) # Use biomart to get orthologs
                # Filter out any low confidence orthologs and any genes 
                #that are not one-to-one orthologs in both directions
                ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
                ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
                ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
                names(ortho_df) <- c('ENSEMBL_MOUSE','CONFIDENCE','ENSEMBL_RAT') 
                ortho_df <- ortho_df %>%
                        select(-CONFIDENCE)
                # Assign the HUGO symbols to a new column
                x <- left_join(x, ortho_df, by = 'ENSEMBL_MOUSE') %>%
                        filter(!is.na(ENSEMBL_MOUSE)) %>%
                        mutate(SYMBOL_RAT = mapIds(org.Rn.eg.db, ENSEMBL_RAT, "SYMBOL", "ENSEMBL"))
        }
        # Return the output
        x
}

#########################################################################################

# Function to take lowercase strings and convert the first letter to uppercase
################################################################################
firstup <- function(x) {
        substr(x, 1, 1) <- toupper(substr(x, 1, 1))
        x
}
################################################################################


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
factor_cols <- c("labelid",
                 "vial_label",
                 "animal.registration.sex",
                 "animal.key.exlt4")
for(fc in factor_cols){
        col_data[[fc]] <- as.factor(col_data[[fc]])
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
               "specimen.collection.d_visit",
               "animal.registration.d_birth")
for(dc in date_cols){
        col_data[[dc]] <- ymd(col_data[[dc]])
}

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
               "specimen.processing.t_freeze",
               "acute.test.howlongshock")
for(tc in time_cols){
        col_data[[tc]] <- col_data[[tc]] %>% as.character() %>% parse_time()
}

# Releveling factors
col_data$animal.key.anirandgroup <- as.character(col_data$animal.key.anirandgroup)
col_data$animal.key.anirandgroup <- factor(col_data$animal.key.anirandgroup, 
                                           levels = c("Control - IPE",
                                                      "Control - 7 hr",
                                                      "Exercise - IPE",
                                                      "Exercise - 0.5 hr",
                                                      "Exercise - 1 hr",
                                                      "Exercise - 4 hr",
                                                      "Exercise - 7 hr",
                                                      "Exercise - 24 hr",
                                                      "Exercise - 48 hr"))


#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script 20200413_rnaseq-kidney-female-temporal-clusters_steep.R

# Circadian Genes (Kidney)
in_file <- paste0(WD,'/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt')
circ_kid <- read.table(file=in_file, sep = '\t', header = TRUE)

#' ## Subset Data for Samples of Interest

#+ Subset Data for Samples of Interest
################################################################################
###########     Subset Data for Samples of Interest   ##########################
################################################################################

#' ##### Female Kidney Samples (39 unique) from 1 batch were filtered prior to normalization
#' # TODO: Create a seperate script that explains why sample 90109015902_SN1 is considered in outlier
fkid <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(sample_key != '90109015902_SN1') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
fkid_counts <- count_data[,fkid] %>% as.matrix()

# SUbset kidney metadata
fkid_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(sample_key != '90109015902_SN1') %>%
        #filter(animal.registration.sex == 'Female') %>%
        filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        filter(!is.na(animal.registration.sex))

row.names(fkid_cols) <- fkid_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))

#' ### Differential Gene Expression: Control 0 hr vs Control 7 hr

#+ Differential Gene Expression: Control 0 hr vs Control 7 hr
################################################################################
###########Differential Gene Expression: Control 0 hr vs Control 7 hr  #########
################################################################################
# Generate an annotation specific to this analysis
fkid_cols <- fkid_cols %>%
        mutate(controls_binary = case_when(animal.key.anirandgroup == 'Control - IPE' ~ 'Control_IPE',
                                           animal.key.anirandgroup == 'Control - 7 hr' ~ 'Control_7hr'))
# Make sure the "control" or "untreated" level is first
fkid_cols$controls_binary <- factor(fkid_cols$controls_binary, 
                                    levels = c("Control_IPE","Control_7hr"))
row.names(fkid_cols) <- fkid_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))

# Create a design formula and load counts and supporting annotation into an S4 object (DESeq infrastructure)
design = ~animal.registration.sex + controls_binary # Primary variable needs to be last.
design = ~controls_binary # Primary variable needs to be last.
title = paste0('Design: ',as.character(design))
dds1 <- DESeqDataSetFromMatrix(countData = fkid_counts,
                               colData = fkid_cols,
                               design = design)

# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
# dds1
zero_n <- dds1[(rowSums(counts(dds1))/ncol(dds1) < 1), ] %>% nrow() %>% as.character()
reads_n <- 10
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

#rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for(n in 1){
        start_time <- Sys.time()
        rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Naive Model (~ 1)") +
        guides(color=guide_legend(title="animal.key.anirandgroup"))

#' ### Adjust for Between Sex Variance

#+ Adjust for Between Sex Variance
################################################################################
########### Adjust for Between Sex Variance  ###################################
################################################################################

# "To adjust for batch effects, we median- centered the expression levels of each transcript within each batch and confirmed, using the correlation matrices, that the batch effects were removed after the adjustment." 
#~ Li, J. Z. et al. Circadian patterns of gene expression in the human brain and disruption in major depressive disorder. Proc. Natl. Acad. Sci. U. S. A. 110, 9950–9955 (2013).

# Here we have 2 Groups: Control - IPE and Control 7 hr; we'll median center these groups to combine the sexes.

M_samples <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Male') %>%
        filter(sample_key != '90109015902_SN1') %>%
        filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
F_samples <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(!is.na(animal.registration.sex)) %>%
        filter(animal.registration.sex == 'Female') %>%
        filter(animal.key.anirandgroup %in% c('Control - IPE', 'Control - 7 hr')) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Select the counts
M_counts <- assay(rld[, M_samples])
F_counts <- assay(rld[, F_samples])

# Median Center data
# Collects median of each row, then subtracts by row medians
M_medians <- apply(M_counts,1,median)
M_centered <- M_counts - M_medians
F_medians <- apply(F_counts,1,median)
F_centered <- F_counts - F_medians
counts_centered <- cbind(M_centered, F_centered)
counts_centered <- counts_centered[, colnames(assay(rld))]
rld_medcen <- rld
assay(rld_medcen) <- counts_centered

#' #### We see just how well duplicate samples correlate regardless of sequencing batch
mypar()
pcaData <- DESeq2::plotPCA(rld_medcen, 
                           intgroup=c("animal.key.anirandgroup",
                                      "animal.registration.sex",
                                      "sample_key"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup,shape=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_label_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("PCA: Naive Model (~ 1)") +
        guides(color=guide_legend(title="animal.key.anirandgroup"))


# Generate a DESeq2 Model
dds <- DESeq(dds)

# Generate a results table and sort by fdr
res <- results(dds, alpha = 0.05,lfcThreshold=0)
res <- res[order(res$padj),]
head(res)

# Number of differentially expressed genes based on q value
table(res$padj < 0.05)

# Generate a summary of the results
summary(res_old) # res_old has design for only primary variable of interest
summary(res)

#' ### Visualize Results: Control 0 hr vs Control 7 hr

#+ Visualize Results: Control 0 hr vs Control 7 hr
################################################################################
########### Visualize Results: Control 0 hr vs Control 7 hr  ###################
################################################################################

# Generate an MA-plot: plots log2 fold change (y-axis) over the mean of normalized counts (x-axis)
DESeq2::plotMA(res, ylim =c(-4,4))

# Generate a pvalue histogram
hist(res$pvalue[res$baseMean > 1], 
     col="grey", border="white", xlab="", ylab="", main="")

# Grab the differentially expressed genes
res_df <- as.data.frame(res)
res_df$ENSEMBL_RAT <- row.names(res_df)
de_df <- res_df %>%
        filter(padj <= 0.05)
de_df$SYMBOL_RAT <- mapIds(org.Rn.eg.db, de_df$ENSEMBL_RAT, "SYMBOL", "ENSEMBL")
de_df$ENTREZ_RAT <- mapIds(org.Rn.eg.db, de_df$ENSEMBL_RAT, "ENTREZID", "ENSEMBL")

# Convert gene symbols to mouse orthologs
de_df <- rat_mouse_ortho(de_df, column = 'ENSEMBL_RAT', direction = 'rat2mouse')

output_n <- de_df$ENSEMBL_RAT %>% unique() %>% length() %>% as.character()
#' ##### High confidence ortholgos were collected from Yan et. al. (Supplementary Table 2) and converted to high confidence rat orthologs with Ensembl.
#' ###### Mouse annotation: GRCm38.p6
#' ###### Rat annotation: Rnor_6.0
#' ###### Stats:
#' * Mouse genes input: `r input_n`
#' * High confidence Rat orthologs output: `r output_n`
#' ###### Steps in ortholog selection:
#' * Demonstration of circadian oscilliations in mouse kidney (from Yan et. al.)
#' * Demonstration of circadian oscilliations in 1 (min Kidney) or more tissues (from Yan et. al.)
#' * Mouse genes were removed that did not have Ensembl gene symbol (or alias) to gene id 
#'     * First choice of 1:many were included
#' * Genes were removed if they did not have a high orthology confidence score between mouse gene id and rat gene id  (binary value 0|1)
#' * Duplicate orthologs were removed
#'     * First choice of 1:many were included

# COllect the gene symbols
de_df$SYMBOL_MOUSE <- mapIds(org.Mm.eg.db, de_df$ENSEMBL_MOUSE, "SYMBOL", "ENSEMBL")
symbol_up <- de_df %>%
        filter(!is.na(SYMBOL_MOUSE)) %>%
        filter(log2FoldChange > 0) %>%
        select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()
symbol_down <- de_df %>%
        filter(!is.na(SYMBOL_MOUSE)) %>%
        filter(log2FoldChange < 0) %>%
        select(SYMBOL_MOUSE) %>% unlist() %>% unique() %>% as.character()

#' ### Pathway Enrichment Analysis (w/ Enrichr)

#+ Visualize Results: Control 0 hr vs Control 7 hr
################################################################################
########### Visualize Results: Control 0 hr vs Control 7 hr  ###################
################################################################################

# To check available databases from enrichr
#dbs <- listEnrichrDbs()
#if (is.null(dbs)) websiteLive <- FALSE
#dbs$libraryName

# Collect the desired pathways to query for enrichment
desired <- c("KEGG_2019_Mouse",
             "GO_Biological_Process_2018",
             "GO_Molecular_Function_2018")

enriched_up <- enrichr(symbol_up, desired)
enriched_down <- enrichr(symbol_down, desired)
n<-1
enrichr_df <- data.frame()
for(db in desired){
        # Collapse the data into flat format
        enriched_up[[n]]$Data.Base <- db
        enriched_up[[n]]$UPDOWN <- 'Up'
        enriched_down[[n]]$Data.Base <- db
        enriched_down[[n]]$UPDOWN <- 'Down'
        enrichr_df <- rbind(enrichr_df,enriched_up[[n]],enriched_down[[n]])
        n <- n + 1
}

# Generate a new column for the number of genes per enriched term
enrichr_df  <- enrichr_df %>% separate(Overlap, sep = "/", into = c("Overlap.N","Total.N"))
enrichr_df$FDR_nlog <- -log(enrichr_df$Adjusted.P.value)

# Filter the results
enrichr_plot <- enrichr_df %>%
        filter(Total.N <= 500) %>%
        filter(Total.N >= 10) %>%
        filter(Adjusted.P.value <= 0.05) %>%
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
ggplot(enrichr_plot, aes(x=reorder(Term, FDR_ord), y=FDR_nlog, fill = UPDOWN)) +
        geom_bar(stat='identity') +
        coord_flip() +
        ylab("-Log(FDR)") +
        xlab("") +
        theme_linedraw() +
        guides(fill=guide_legend(title="Up/Down")) +
        scale_fill_manual(values = c('Down' = "#00BFC4",'Up' = "#F8766D")) +
        theme(panel.background = element_blank(),
              plot.background = element_blank(),
              strip.background = element_blank(),
              axis.text = element_text(size=12, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              axis.title=element_text(size=12,face="bold"),
              strip.text = element_blank()) +
        theme(panel.grid.major = element_blank())

#' ## Additional Test for Enrichment: Circadian Gene Set

#+ Additional Test for Enrichment: Circadian Gene Set
################################################################################
########### Additional Test for Enrichment: Circadian Gene Set  ################
################################################################################

# Steps:
# Logistic regression for enrichment of circadian gene set
# Generate updated list of circadian genes
# Geomtile plot of shared genes (bi-heatmap)

# Create a dataframe for logistic regression
ENSEMBL_RAT <- row.names(assay(rld))
de_genes <- de_df %>%
        filter(!is.na(ENSEMBL_RAT)) %>%
        select(ENSEMBL_RAT) %>% unlist() %>% unique() %>% as.character()
logreg_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
logreg_df <- mutate(logreg_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 1, 0))
circ_genes <- ENSEMBL_RAT[ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT]
# Generate a column for differential gene epression status
logreg_df <- mutate(logreg_df, DE = ifelse(ENSEMBL_RAT %in% de_genes, 1, 0))
# Remove gene ids and set them to rownames
row.names(logreg_df) <- logreg_df$ENSEMBL_RAT
logreg_df <-  logreg_df %>% select(-ENSEMBL_RAT)

# Adjust column objects
factor_cols <- colnames(logreg_df)
for(fc in factor_cols){
        logreg_df[[fc]] <- as.factor(logreg_df[[fc]])
}

# Fit a logistic regression model
mod <- glm(formula = CIRC ~ DE, 
           family = binomial(link = logit), 
           data = logreg_df)
summary(mod)

# Perform a Chi-square analysis
#######################
# Create a proper data structure for chi squared

de_circ <- logreg_df %>%
        filter(CIRC == 1 & DE == 1) %>%
        nrow()
de_no_circ <- logreg_df %>%
        filter(CIRC == 0 & DE == 1) %>%
        nrow()
no_de_circ <- logreg_df %>%
        filter(CIRC == 1 & DE == 0) %>%
        nrow()
no_de_no_circ <- logreg_df %>%
        filter(CIRC == 0 & DE == 0) %>%
        nrow()
x2_df <- data.frame("CIRC" = c(de_circ, no_de_circ), 
                    "NON-CIRC" = c(de_no_circ, no_de_no_circ))
row.names(x2_df) <- c('DE','NON-DE')

# Create the proper data structure for visualization of chi square with mosaic plot
mosaic_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
mosaic_df <- mutate(mosaic_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 'CIRC', "NON-CIRC"))
# Generate a column for DE gene status
mosaic_df <- mutate(mosaic_df, DE = ifelse(ENSEMBL_RAT %in% de_genes, 'DE', "NON-DE"))

# Generate a table of results
mosaic_tbl <- table(mosaic_df$CIRC, mosaic_df$DE, dnn=c("CIRC","DE"))

# Perform Chi-squared
#sink('')
(( c=chisq.test(x2_df, simulate.p.value = FALSE) ))
(( fisher.test(x2_df, simulate.p.value = FALSE) ))

#sink()
c$observed
round(c$expected)

# Visualize the chisquare analysis with a mosaic plot
mosaic(~ CIRC + DE, data = mosaic_tbl,
       shade = TRUE, legend = TRUE)

# Update the shared list of genes
##############################################
logreg_df$ENSEMBL_RAT <- row.names(logreg_df)
de_kid_circ <- logreg_df %>%
        filter(CIRC == 1 & DE == 1) %>%
        select(ENSEMBL_RAT) %>% unlist() %>% as.character()

kegg_circ <- enrichr_df %>%
        filter(Term == 'Circadian rhythm') %>%
        select(Genes) %>%
        unlist() %>% as.character() %>% 
        str_split(';') %>% unlist() %>% unique() %>%
        tolower() %>% firstup() %>% as.data.frame()
names(kegg_circ) <- 'SYMBOL_MOUSE'
# Convert the mouse symbols to rat ensembl gene ids
kegg_circ$ENSEMBL_MOUSE = mapIds(org.Mm.eg.db, 
                                 as.character(kegg_circ$SYMBOL_MOUSE), 
                                 "ENSEMBL", "SYMBOL")
kegg_circ <- rat_mouse_ortho(kegg_circ, column = 'ENSEMBL_MOUSE', direction = 'mouse2rat')
CIRC_GENES <- c(kegg_circ$ENSEMBL_RAT[kegg_circ$ENSEMBL_RAT %!in% de_kid_circ], de_kid_circ) %>% 
        unique()

#' ## Visualize Circadian Genes Differentially Expressed

#+ Visualize Circadian Genes Differentially Expressed
################################################################################
########### Visualize Circadian Genes Differentially Expressed  ################
################################################################################

# Steps:
#-Barplot of DE expression between C0 and C7
#-Geom tile heatmap of circadian genes in other pathways
#-Heatmap and k-means/SOM clustered circadian genes in control contrast and exercsie groups
#- Same analysis but for downstream targets of CIRC genes
#-Next to barplot show expression of identical gene in exercise groups
#-Perform a ttest between peaks of both genes
#-Venn diagram and supporting qqplot of DE genes between analyses
#Examination of genes lost
#Can it be corrected for?







