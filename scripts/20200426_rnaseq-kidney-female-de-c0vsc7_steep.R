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

# Function to relabel RNASeq read names to orthologs
###############################################################################################
# mmusculus_gene_ensembl: Mouse genes (GRCm38.p6)
# rnorvegicus_gene_ensembl: Rat genes (Rnor_6.0)
rat2mouse_ortho <- function(x, column = "ENSEMBL_RAT") {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        if ( column != 'ENSEMBL_RAT' ){
                stop("'column' must be 'ENSEMBL_RAT'", class.= FALSE)
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
        # Create ortholog table
        ortho_df <- getLDS(attributes=c("ensembl_gene_id","mmusculus_homolog_orthology_confidence"),
                           filters="ensembl_gene_id", 
                           values = x[[column]], 
                           mart=mart_rn_ens,
                           attributesL=c("ensembl_gene_id"), 
                           martL=mart_mm_ens) # Use biomart to get orthologs
        # Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
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
        x
}
#########################################################################################

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
design = ~ controls_binary # Primary variable needs to be last less.
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
de_df <- rat2mouse_ortho(de_df, column = 'ENSEMBL_RAT')

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
symbol_out <- de_df %>%
        filter(!is.na(SYMBOL_MOUSE)) %>%
        select(SYMBOL_MOUSE) %>% unlist() %>% as.character()

out_file <- paste0(WD,'/data/20200426_rnaseq-kidney-C0vsC7-DE-genes-q05-lfc0-mouse-symbol_steep.txt')
write.table(symbol_out, file = out_file, sep = '\n', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Filter and format input for LRPath
#LRPath_out <- de_df %>%
#        filter(!is.na(RAT_ENTREZ)) %>%
#        select(geneid = RAT_ENTREZ,
#               PValue = pvalue,
#               logFC = log2FoldChange,
#               norm_avg_readcount = baseMean)

#out_file <- paste0(WD,'/data/20200426_rnaseq-kidney-C0vsC7-DE-genes-q05-lfc0_steep.txt')
#write.table(LRPath_out, file = out_file, sep = '\t', 
#            quote = FALSE, row.names = FALSE, col.names = TRUE)


#' ### Pathway Enrichment Analysis (w/ Enrichr)

#+ Visualize Results: Control 0 hr vs Control 7 hr
################################################################################
########### Visualize Results: Control 0 hr vs Control 7 hr  ###################
################################################################################

# Manipulate the file in linux
# sed "s/\#//g" 20200426_rnaseq-kidney-C0vsC7-DE-genes-q05-lfc0-lrpath_steep.txt | tr "'" "___" | sed 's/___//g' > 20200426_rnaseq-kidney-C0vsC7-DE-genes-q05-lfc0-lrpath_steep.txt

# Load the Enrichr results
in_file <- paste0(WD,'/data/20200426_rnaseq-kidney-C0vsC7-DE-genes-q05-lfc0-lrpath_steep.txt')
lrpath_res <- read.table(file = in_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
names(lrpath_res)
head(lrpath_res)

# Filter the results
dat <- lrpath_res %>%
        select(Name,ConceptType, Genes,Direction,OddsRatio,P.Value,FDR, SigGenes) %>%
        #filter(Genes < 2000 & Genes > 30) %>%
        #filter(FDR <= 0.05) %>%
        group_by(ConceptType, Direction) %>%
        arrange(FDR) %>%
        dplyr::slice(1:30) %>%
        ungroup()

#dot plot
ggplot(dat[dat$ConceptType=="KEGG" & dat$Direction == "up",], 
       aes(x = Name, y = Genes, size = OddsRatio, col = P.Value)) +
        geom_point() +
        coord_flip() +
        scale_color_gradient(low = "red", high = "black") +
        theme(axis.text = element_text(size = rel(1)))

#dot plot
ggplot(dat[dat$ConceptType=="GOBP" & dat$Direction == "down",], 
       aes(x = Name, y = X.Genes, size = 1/OddsRatio, col = P.Value)) +
        geom_point() +
        coord_flip() +
        scale_color_gradient(low = "blue", high = "black") +
        theme(axis.text = element_text(size = rel(1)))

# Visualize results with barplot
ggplot(df, aes(x=Term, y=log, fill = up_dn)) +
        geom_bar(stat='identity') +
        coord_flip() +
        ylab("-Log(FDR q-value)") +
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

