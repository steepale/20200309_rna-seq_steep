#'---
#' title: "PASS1A Rat Kidney: -- Temporal Clusters"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200413"
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
#' TODO:
#' 
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
#BiocManager::install("EBImage")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen","som","caret")
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
mouse2rat_ortho <- function(x) {
        # Ensure 'x' is a data.frame
        if ( class(x) != "data.frame" ) {
                stop("'x' must be a data frame", class.= FALSE)
        }
        
        # Load requirements
        library(biomaRt)
        library(purrr)
        library(dplyr)
        # Load in annotations
        mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
        mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
        # Create ortholog table
        ortho_df <- getLDS(attributes=c("ensembl_gene_id","rnorvegicus_homolog_orthology_confidence"),
                           filters="ensembl_gene_id", 
                           values = x$ENSEMBL_MOUSE, 
                           mart=mart_mm_ens,
                           attributesL=c("ensembl_gene_id"), 
                           martL=mart_rn_ens) # Use biomart to get orthologs
        # Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
        ortho_df <- ortho_df[ortho_df$Rat.orthology.confidence..0.low..1.high. == '1',]
        ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
        ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
        names(ortho_df) <- c('ENSEMBL_MOUSE','CONFIDENCE','ENSEMBL_RAT') 
        ortho_df <- ortho_df %>%
                select(-CONFIDENCE)
        
        # Assumes that 'x' has ensembl chicken gene ID's as rownames
        # Ensure that only chicken genes appear in the matrix
        #x <- x[startsWith(rownames(x), "ENSGALG"),]
        # Assign the HUGO symbols to a new column
        x <- left_join(x, ortho_df, by = "ENSEMBL_MOUSE") %>%
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

### Determine which control samples are male or female
# Get the list of genes on the W chromosome

# Construct your own personal galgal5 reference genome annotation
# Construct from gtf file from Ensembl (same file used in mapping)
ens_gtf <- paste0(WD,'/data/Rattus_norvegicus.Rnor_6.0.96.gtf')
Rn_TxDb <- makeTxDbFromGFF(ens_gtf,
                           format=c("gtf"),
                           dataSource="Ensembl_Rattus6_gtf",
                           organism="Rattus norvegicus",
                           taxonomyId=NA,
                           circ_seqs=DEFAULT_CIRC_SEQS,
                           chrominfo=NULL,
                           miRBaseBuild=NA,
                           metadata=NULL)

# Define Female specific sex genes (X chromosome)
# To examine chromosome names
seqlevels(Rn_TxDb)[1:23]
# Extract genes as GRanges object, then names
X_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("X")))
# Collect ensembl gene ids for female specific genes
X_ens_id <- names(X_genes_gr)
# Examine the gene symbols
X_sym <- mapIds(org.Rn.eg.db, names(X_genes_gr), "SYMBOL", "ENSEMBL")
# Extract genes as GRanges object, then names
Y_genes_gr <- genes(Rn_TxDb, columns = "TXCHROM", filter = list(tx_chrom=c("Y")))
# Collect ensembl gene ids for female specific genes
Y_ens_id <- names(Y_genes_gr)
sex_ens_id <- c(X_ens_id,Y_ens_id)
# Examine the gene symbols
Y_sym <- mapIds(org.Rn.eg.db, names(Y_genes_gr), "SYMBOL", "ENSEMBL")

#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep

#+ Circadian Gene Sets
################################################################################
#####     Circadian Gene Sets       ############################################
################################################################################

# Circadian Genes
in_file <- paste0(WD,'/data/20080516_mouse-tissue-circadian-genes_yan.txt')
circ_df <- read.table(in_file,sep = '\t', header = TRUE,check.names = FALSE)
# Adjust gene symbol names
names(circ_df)[1] <- "SYMBOL_MOUSE"
circ_kid <- circ_df %>%
        select(SYMBOL_MOUSE, KID, NUM.TISSUE, RANGE.P, PEAK.MEAN) %>%
        filter(!is.na(KID)) %>%
        filter(NUM.TISSUE >= 1)
input_n <- circ_kid$SYMBOL_MOUSE %>% unique() %>% length() %>% as.character()

# Some genes were annotated with synonyms by Yan et al. We manually convert these synonyms to gene symbols
circ_kid$SYMBOL_MOUSE <- as.character(circ_kid$SYMBOL_MOUSE)
circ_kid$ENSEMBL_MOUSE <- mapIds(org.Mm.eg.db, circ_kid$SYMBOL_MOUSE, "ENSEMBL", "SYMBOL")
# Useful trick, search symbol in aliases and grab symbol
circ_kid$ALIAS_MOUSE <- mapIds(org.Mm.eg.db, 
                               circ_kid$SYMBOL_MOUSE, "SYMBOL", "ALIAS",
                               multiVals = "first")

# With a slow for loop, iterate through rows of dataframe and adjust nomenclature scheme as needed
for( r in 1:nrow(circ_kid)){
        ENSEMBL <- circ_kid[r,'ENSEMBL_MOUSE']
        ALIAS <- circ_kid[r,'ALIAS_MOUSE']
        if(is.na(ENSEMBL)){
                if(!is.na(ALIAS)){
                        # Remember that ALIAS is actually from Ensembl's SYMBOL column
                        ENSEMBL_RPL <- mapIds(org.Mm.eg.db, ALIAS, "ENSEMBL", "SYMBOL")
                        circ_kid[r,'ENSEMBL_MOUSE'] <- ENSEMBL_RPL
                        circ_kid[r,'SYMBOL_MOUSE'] <- ALIAS
                }
        }
}

circ_kid <- circ_kid %>%
        filter(!is.na(ENSEMBL_MOUSE)) %>%
        select(-ALIAS_MOUSE)
circ_kid2 <- mouse2rat_ortho(circ_kid)

output_n <- circ_kid2$ENSEMBL_RAT %>% unique() %>% length() %>% as.character()
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

# Cleanup
circ_kid <- circ_kid2
rm(circ_kid2)

# Make one custom adjustment (Biomart did not catch this annotation)
circ_kid <- circ_kid %>%
        mutate(SYMBOL_RAT = ifelse(
                ENSEMBL_RAT %in% c('ENSRNOG00000060956','ENSRNOG00000047309',
                                   'ENSRNOG00000000419','ENSRNOG00000018536'), 
                SYMBOL_MOUSE, SYMBOL_RAT))

# Collect kidney circadian kidney genes
circ_kid_ens <- circ_kid %>%
        select(ENSEMBL_RAT) %>% unlist()

#' #### Data Save: Mouse and Rat Orthologus Circadian Genes (Kidney)
#' Circadian Kidney File: motrpac/20200309_rna-seq_steep/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt

#+ Data Save: Mouse and Rat Orthologus Circadian Genes (Kidney)
################################################################################
######## Data Save: Mouse and Rat Orthologus Circadian Genes (Kidney) ###########
################################################################################

# Circadian Genes (Kidney)
out_file <- paste0(WD,'/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt')
# write.table(circ_kid, file=out_file, row.names = FALSE, quote = FALSE, sep = '\t')

#' ## Normalization for Sequencing Depth (Female Kidneys)

#+ Normalization for Sequencing Depth (Female Kidneys)
################################################################################
#####     Normalization for Sequencing Depth  (Female Kidneys)   ###############
################################################################################

#' ##### Female Kidney Samples (39 unique) from 1 batch were filtered prior to normalization
#' One outlier removed: TODO: Document removal of outlier in seperate script
fkid <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(animal.registration.sex == 'Male') %>%
        filter(!is.na(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
fkid_counts <- count_data[,fkid] %>% as.matrix()

mean_males <- rowMeans(fkid_counts_m)
mean_females <- rowMeans(fkid_counts_f)

mean_avg <- (mean_males + mean_females)/2

# All samples unique?
#col_data %>%
#  filter(Tissue == 'Kidney') %>%
#  select(vial_label) %>% unlist() %>% unique() %>% length()

# SUbset kidney metadata
fkid_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        #filter(animal.registration.sex == 'Female')
        filter(!is.na(animal.registration.sex))
#filter(sample_key != '90042016803_SF1')
row.names(fkid_cols) <- fkid_cols$sample_key
names(fkid_cols)
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))
colnames(fkid_cols)
#design = ~ sex + animal.key.exlt4 # Primary variable needs to be last
#design = ~ 1
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object

fkid_cols <- fkid_cols %>%
        mutate(sex = ifelse(animal.registration.sex == 'Female', 1, 0))

dds1 <- DESeqDataSetFromMatrix(countData = fkid_counts,
                               colData = fkid_cols,
                               design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
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



#' ### Transformation (Female Kidneys)

#+ Transformation (Female Kidneys)
################################################################################
#####     Transformation (Female Kidneys)      #################################
################################################################################

#' ##### Log2 normalization (with size factor and pseudo count)"
log_norm_counts <- log2(counts(dds, normalized=TRUE) + 1)
#rld <- DESeq2::vst(dds)
#' Regularized Log (rlog) Transform
for( n in 1){
        start_time <- Sys.time()
        rld <- DESeq2::rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
#' #### Here we visualize the counts comparing log2 transform and rlog transform
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] +1), cex=0.2, main = "log2 Transform", ylim=c(0,22)) # Log2 Transform
boxplot(assay(rld)[rs > 0,], cex=0.2, main = "rlog Transform", ylim=c(0,22)) # Regularized Log Transform

#' ##### Returning to our prior plot of comparing expression of 2 samples, we can see the the variance is now more stable in the rlog transform.
mypar(1,1)
plot(assay(rld)[,c(1,2)],
     xlab = "Sample 1 rlog Gene Counts", 
     ylab = "Sample 2 rlog Gene Counts", 
     cex = 0.3, main = "rlog Transform")







#' TODO: QC Check, Jun recommends demonstrating the number of features at different (perhaps dec) -tiles from different transforms. We need to double check that the transform did not remove features or shrink low expression features so extremely to demonstrate zero variance.

#' #### Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group)

#+ Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group)
################################################################################
#Unsupervised Clustering w/ SOMs & K-means (Ordered by Exercsie/Control Group) #
################################################################################

#' #### Steps in Analysis:
#' * Orient Count Matrix by sample (group) for temporal visualization
#' * Perform a variant test to remove genes that do not vary across groups
#'     * Ensure that variant filter parameters are met if test is parametric
#' * Scale and center data to have a mean of zero and a variance of 1
#' * Perform SOM clustering

# Create a vector for exercise/group levels
ec_levels <- c("Exercise - IPE",
               "Exercise - 0.5 hr",
               "Exercise - 1 hr",
               "Exercise - 4 hr",
               "Exercise - 7 hr",
               "Exercise - 24 hr",
               "Exercise - 48 hr",
               "Control - IPE",
               "Control - 7 hr")
# Create a vector for exercise/group time of death levels
ec_time_levels <- c("13:15 - 13:35",
                    "11:55 - 12:15",
                    "13:55 - 14:15",
                    "14:35 - 14:55",
                    "17:35 - 17:55",
                    "10:35 - 10:55",
                    "11:15 - 11:35",
                    "9:55 - 10:15",
                    "16:55 - 17:15")

# Order samples by exercise group
group_order <- fkid_cols %>%
        arrange(factor(animal.key.anirandgroup, 
                       levels = ec_levels), 
                desc(specimen.collection.t_death), 
                desc(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Order of exercise groups
group_order2 <- fkid_cols %>%
        arrange(factor(animal.key.anirandgroup, 
                       levels = ec_levels), 
                desc(specimen.collection.t_death), 
                desc(animal.registration.sex)) %>%
        select(animal.key.anirandgroup) %>% unlist() %>% as.character()

# Arrange matrix in order of exercise and control groups
# All expression values per gene are normalized to generate a level playing field. Means are subtracted by values and the resulting difference is divided by the standard deviation.
som_xgroup <- assay(rld)[,group_order]
som_xgroup <- som_xgroup %>% t() %>% data.frame()

# Test for multiv
#x <- as.numeric(MVN::mvn(som_xgroup, mvnTest = "hz")$univariateNormality$`p value`)
#(p.adjust(x, method = "BH") > 0.05) %>% table()
som_xgroup$animal.key.anirandgroup <- group_order2

# MANOVA test
# Generate a formula
dependents <- colnames(som_xgroup)[colnames(som_xgroup) %!in% c("animal.key.anirandgroup")]
form <- as.formula(paste0("cbind(",paste(dependents, collapse=","),")", "~animal.key.anirandgroup"))
# Perform MANOVA
res_man <- manova(form, data = som_xgroup)
res_sum <- summary.aov(res_man)

# Organize the index of genes by pvalue
df_x <- data.frame(1:(ncol(som_xgroup)-1))
names(df_x) <- 'idx'
pvec <- vector()
for(nr in 1:(ncol(som_xgroup)-1)){
        pval <- res_sum[[nr]]$`Pr(>F)`[1]
        pvec <- c(pvec, pval)
}
df_x$pval <- pvec

# Choose row index with significant p value
c_idx <- df_x %>%
        filter(pval <= 0.05) %>%
        select(idx) %>% unlist() %>% as.numeric()
# Select only genes that show significant variance across timepoints
som_xgroup <- som_xgroup[,c_idx]
dim(som_xgroup)
# Data should be normalized to have a mean zero and a variance of 1 (between -1 and 1)
preproc1 <- preProcess(som_xgroup, method=c("center", "scale"))
norm1 <- predict(preproc1, som_xgroup)
som_mat <- t(norm1) %>% as.matrix()

# Create a SOM (with som)
set.seed(666)
som_group <- som::som(som_mat,6,5)
plot(som_group, ylim=c(-2,2))

#som::som(som_xgroup,10,10) %>% plot()
table(som_group$visual[,1:2])

# Perform k-means clustering
k6<-kmeans(som_mat,6)$cluster

# Create an annotation for the heatmap
ann_df <- col_data %>%
        filter(sample_key %in% colnames(som_mat)) %>%
        select(animal.key.anirandgroup)
ann_df$animal.key.anirandgroup <- factor(ann_df$animal.key.anirandgroup, levels = ec_levels)
ann_df <- ann_df %>% arrange(animal.key.anirandgroup)
row.names(ann_df) <- colnames(som_mat)
ann_colors = list(
        animal.key.anirandgroup = 
                c("Exercise - IPE" = "gold",
                "Exercise - 0.5 hr" = "darkgoldenrod1",
                "Exercise - 1 hr" = "orange",
                "Exercise - 4 hr" = "darkorange",
                "Exercise - 7 hr" = "darkorange2",
                "Exercise - 24 hr" = "darkorange3",
                "Exercise - 48 hr" = "darkorange4",
                "Control - IPE" = "steelblue1",
                "Control - 7 hr" = "steelblue4"))

# Create heatmaps from kmeans
heat_plots <- vector('list', 6)
heat_plots[[1]] <- pheatmap(som_mat[k6==1,], 
                                  annotation_col = ann_df,
                                  annotation_colors = ann_colors,
                                  fontsize = 8,
                                  show_rownames=FALSE,
                                  show_colnames=FALSE,
                                  color = rev(brewer.pal(n = 9, name ="RdBu")),
                                  cluster_cols = FALSE,
                                  cluster_rows = FALSE,
                                  legend=F,
                                  annotation_legend = FALSE)
heat_plots[[2]] <- pheatmap(som_mat[k6==2,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[3]] <- pheatmap(som_mat[k6==3,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[4]] <- pheatmap(som_mat[k6==4,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[5]] <- pheatmap(som_mat[k6==5,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)
heat_plots[[6]] <- pheatmap(som_mat[k6==6,], 
                            annotation_col = ann_df,
                            annotation_colors = ann_colors,
                            fontsize = 8,
                            show_rownames=FALSE,
                            show_colnames=FALSE,
                            color = rev(brewer.pal(n = 9, name ="RdBu")),
                            cluster_cols = FALSE,
                            cluster_rows = FALSE,
                            legend=F,
                            annotation_legend = FALSE)

# Clusters are numbers in numerical order from 1 to 30 from bottom-left to top-right
# Identify the circadian rhythm genes that belong to each cluster
i <- 1
cluster <- list()
circ_cl <- list()
circ_cldf <- list()
ENSEMBL_RAT <- vector()
for(yn in 0:4){
        for(xn in 0:5){
                # Collect gene ids for each cluster
                c_rn <- som_group$visual[(som_group$visual$x == xn & som_group$visual$y == yn),] %>%
                        row.names() %>% as.numeric()
                cluster[[i]] <- som_group$data[c_rn,] %>% row.names()
                # Collect all gene ids
                ENSEMBL_RAT <- c(ENSEMBL_RAT, cluster[[i]])
                # collect circadian gene ids per cluster
                i <- i + 1
        }
}

# Create a geom_tile to mimic the som and visualize the kmeans clusters
grad_plots <- vector('list', 6)
kmeans_cluster <- list()
for(kn in 1:6){
        # Kmeans cluster 1
        x <- 0:5
        y <- 0:4
        tile_df <- expand.grid(X=x,Y=y)
        # Kmeans cluster genes are in which som clusters
        kgenes <- names(k6[k6==kn])
        kmeans_cluster[[kn]] <- kgenes
        circ_cl[[kn]] <- circ_kid %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(ENSEMBL_RAT) %>% unlist(use.names = FALSE)
        # collect additional information for circadian genes in cluster
        circ_cldf[[kn]] <- circ_kid %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(NUM.TISSUE, ENSEMBL_RAT, SYMBOL_RAT)
        Z <- vector()
        for(i in 1:30){
                Z <- c(Z, (kgenes[kgenes %in% cluster[[i]]] %>% length()))
        }
        tile_df$Z <- Z
        # Generate the kmeans plot
        grad_plots[[kn]] <- local({
                kn <- kn
                kmsom <- ggplot(tile_df, aes(x = X, y = Y, fill =Z)) +
                        geom_raster(interpolate=TRUE) +
                        scale_fill_gradient2(low="navy", mid="white", high="red", 
                                             midpoint=30, limits=range(tile_df$Z),
                                             guide = FALSE) +
                        theme(axis.text        = element_blank(),
                              axis.ticks       = element_blank(),
                              axis.title       = element_blank(),
                              panel.background = element_blank())
                print(kmsom)
        })
}

# Create a dataframe for logistic regression
logreg_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
logreg_df <- mutate(logreg_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 1, 0))
# Generate 6 columns for each kmeans cluster status
for(i in 1:6){
        new_col_name <- paste0('C',i)
        logreg_df <- logreg_df %>% 
                mutate(!!sym(new_col_name) := ifelse(ENSEMBL_RAT %in% kmeans_cluster[[i]], 1, 0))
}

# Remove gene ids and set them to rownames
row.names(logreg_df) <- logreg_df$ENSEMBL_RAT
logreg_df <-  logreg_df %>% select(-ENSEMBL_RAT)

# Adjust column objects
factor_cols <- colnames(logreg_df)
for(fc in factor_cols){
        logreg_df[[fc]] <- as.factor(logreg_df[[fc]])
}
str(logreg_df)

# Fit a logistic regression model
mod <- glm(formula = CIRC ~ C1+C2+C3+C4+C5+C6, 
           family = binomial(link = logit), 
           data = logreg_df)
summary(mod)

# Perform a Chi-square analysis
#######################
# Create a proper data structure for chi squared
x2_list <- list()
for(i in 1:6){
        n_clust <- length(kmeans_cluster[[i]])
        n_circ <- length(circ_cl[[i]])
        n_no_circ <- n_clust - n_circ
        x2_list[[i]] <- c(n_circ, n_no_circ)
}

x2_df <- as.data.frame(x2_list)
row.names(x2_df) <- c('CIRC','NON-CIRC')
for(i in 1:6){
        new_col_name <- paste0('C',i)
        colnames(x2_df)[i] <- new_col_name
}

# Create the proper data structure for visualization of chi square with mosaic plot
mosaic_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
mosaic_df <- mutate(mosaic_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 'CIRC', "NON-CIRC"))
# Generate a Cluster columcluster status
mosaic_df <- mosaic_df %>%
        mutate(CLUSTER = case_when(ENSEMBL_RAT %in% kmeans_cluster[[1]] ~ 1,
                                   ENSEMBL_RAT %in% kmeans_cluster[[2]] ~ 2,
                                   ENSEMBL_RAT %in% kmeans_cluster[[3]] ~ 3,
                                   ENSEMBL_RAT %in% kmeans_cluster[[4]] ~ 4,
                                   ENSEMBL_RAT %in% kmeans_cluster[[5]] ~ 5,
                                   ENSEMBL_RAT %in% kmeans_cluster[[6]] ~ 6))
# Generate a table of results
mosaic_tbl <- table(mosaic_df$CIRC, mosaic_df$CLUSTER, dnn=c("CIRC","CLUSTER"))

# Perform Chi-squared
#sink('')
(( c=chisq.test(x2_df, simulate.p.value = TRUE) ))
(( fisher.test(x2_df, simulate.p.value = TRUE) ))

#sink()
c$observed
round(c$expected)

# Visualize the chisquare analysis with a mosaic plot
mosaic(~ CIRC + CLUSTER, data = mosaic_tbl,
       shade = TRUE, legend = TRUE)

grad_plots[[4]]

# Clusters enriched WITH circadian genes
# 4
grad_plots[[4]]
heat_plots[[4]]
mypar()
plot(som_group, ylim=c(-2,2))
plot(som_group)

# Cluster enriched WITHOUT circadian genes
# 5
grad_plots[[5]]
heat_plots[[5]]
plot(som_group, ylim=c(-2,2))

# Add a column for time of death binned specifically to exercise group
# TODO: Add this to original metadata generation script, then update excel documentation
#col_data <- col_data %>%
#        mutate(specimen.collection.t_death_bins_ec = case_when(
#                animal.key.anirandgroup == "Exercise - IPE" ~ "13:15 - 13:35",
#                animal.key.anirandgroup == "Exercise - 0.5 hr" ~ "11:55 - 12:15",
#                animal.key.anirandgroup == "Exercise - 1 hr" ~ "13:55 - 14:15",
#                animal.key.anirandgroup == "Exercise - 4 hr" ~ "14:35 - 14:55",
#                animal.key.anirandgroup == "Exercise - 7 hr" ~ "17:35 - 17:55",
#                animal.key.anirandgroup == "Exercise - 24 hr" ~ "10:35 - 10:55",
#                animal.key.anirandgroup == "Exercise - 48 hr" ~ "11:15 - 11:35",
#                animal.key.anirandgroup == "Control - IPE" ~ "9:55 - 10:15",
#                animal.key.anirandgroup == "Control - 7 hr" ~ "16:55 - 17:15"))



# Plot the circadian genes in cluster 4
# Create a dataframe for the plot
df1 <- data.frame(t(som_mat))
df1$sample_key<- row.names(df1)
df2 <- col_data %>%
        select(specimen.collection.t_death, 
               animal.key.anirandgroup)
df2$sample_key <- row.names(df2)
# Adjust the column names to be symbols for ease of plot interpretation
df_plot <- left_join(df1, df2, by = "sample_key")
genes <- colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))]
symbols <- mapIds(org.Rn.eg.db, genes, "SYMBOL", "ENSEMBL")
colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))] <- symbols

# Melt the plot
melt_plot <- reshape2::melt(df_plot, id.vars = c("specimen.collection.t_death", 
                                       "animal.key.anirandgroup"))
# Adjust columns and factor levels
melt_plot$value <- as.numeric(melt_plot$value)
melt_plot$animal.key.anirandgroup <- factor(melt_plot$animal.key.anirandgroup, 
                                            levels = ec_levels)

# To re-examine SOM
plot(som_group, ylim=c(-2,2))

# Plot the circadian genes for cluster of choice
# Show all plots they are enriched with and without in powerpoint
plot(grad_plots[[1]])
heat_plots[[1]]
melt_plot %>%
        filter(variable %in% circ_cldf[[2]]$SYMBOL_RAT) %>%
        ggplot(aes(x = as.integer(animal.key.anirandgroup), 
           y = value, 
           color = variable),
           group = "1") +
        geom_point(alpha = 0.05) +
        #geom_line(aes(group = "1")) +
        stat_smooth(alpha = 0.02, se = F) +
        ylab("rlog Transformed Expression") +
        xlab("Exercise/Control Groups") +
        scale_x_continuous(breaks = seq(1,9,by=1),
                           labels=ec_levels) +
        theme(legend.position = "none")


# TODO: Turn this into a function rat2mouse
################################################################################
x <- row.names(rld)
# Load in annotations
mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
# Create ortholog table
ortho_df <- getLDS(attributes=c("ensembl_gene_id","mmusculus_homolog_orthology_confidence"),
                   filters="ensembl_gene_id", 
                   values = x, 
                   mart=mart_rn_ens,
                   attributesL=c("ensembl_gene_id"), 
                   martL=mart_mm_ens) # Use biomart to get orthologs
str(ortho_df)

# Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
ortho_df <- ortho_df %>%
        select(-CONFIDENCE)
################################################################################

# Pathway enrichment for all genes, CIRC genes, and NON-CIRC genes in cluster 5
all_out <- kmeans_cluster[[1]]
circ_out <- circ_cl[[1]]

# Convert genes to mouse orthologs for pathway enrichment
all_out <- ortho_df %>%
        filter(ENSEMBL_RAT %in% all_out) %>%
        select(ENSEMBL_MOUSE) %>% unlist() %>% as.character()
all_out <- mapIds(org.Mm.eg.db, all_out, "SYMBOL", "ENSEMBL") %>% unlist(use.names = FALSE)

# Save these list to file
out_file <- paste0(WD,'/data/20200420_rnaseq-kidney-kmeans-ord6-clust1-all_steep.txt')
write.table(all_out, file = out_file, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep ='\n')

# Convert genes to mouse orthologs for pathway enrichment
circ_out <- ortho_df %>%
        filter(ENSEMBL_RAT %in% circ_out) %>%
        select(ENSEMBL_MOUSE) %>% unlist() %>% as.character()
circ_out <- mapIds(org.Mm.eg.db, circ_out, "SYMBOL", "ENSEMBL") %>% unlist(use.names = FALSE)

# Save these list to file
out_file <- paste0(WD,'/data/20200420_rnaseq-kidney-kmeans-ord2-clust4-circ_steep.txt')
write.table(circ_out, file = out_file, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep ='\n')

# Pathway enrichment for all genes, CIRC genes, and NON-CIRC genes in any cluster that resembles any sort of circadian effect


#' #### Unsupervised Clustering w/ SOMs & K-means (Ordered by TOD Bins Group)

#+ Unsupervised Clustering w/ SOMs & K-means (Ordered by TOD Bins Group)
################################################################################
#Unsupervised Clustering w/ SOMs & K-means (Ordered by TOD Bins Group) #
################################################################################

# Order samples by time of death (time of day)
tod_order <- fkid_cols %>%
        arrange(specimen.collection.t_death) %>%
        select(sample_key) %>% unlist() %>% as.character()
tod_order2 <- fkid_cols %>%
        arrange(specimen.collection.t_death_bins.type) %>%
        select(specimen.collection.t_death_bins.type) %>% unlist() %>% as.character()

# Create a vector for TOD Bin levels
tod_levels <- c("10.00-11.25",
                "11.75-12.50",
                "13.25-14.00",
                "14.25-15.00",
                "17.00-18.00")

# Arrange matrix in order of exercise and control groups
# All expression values per gene are normalized to generate a level playing field. Means are subtracted by values and the resulting difference is divided by the standard deviation.
som_tod <- assay(rld)[,tod_order]
som_tod <- som_tod %>% t() %>% data.frame()
# Test for multiv
#x <- as.numeric(MVN::mvn(som_tod, mvnTest = "hz")$univariateNormality$`p value`)
#(p.adjust(x, method = "BH") > 0.05) %>% table()
som_tod$specimen.collection.t_death_bins.type <- factor(tod_order2, 
                                                        levels = tod_levels)

# MANOVA test
# Generate a formula
dependents <- colnames(som_tod)[colnames(som_tod) %!in% c("specimen.collection.t_death_bins.type")]
form <- as.formula(paste0("cbind(",paste(dependents, collapse=","),")", "~specimen.collection.t_death_bins.type"))
# Perform the MANOVA Test
res_man <- manova(form, data = som_tod)
res_sum <- summary.aov(res_man)

# Organize the index of genes by pvalue
df_tod <- data.frame(1:(ncol(som_tod)-1))
names(df_tod) <- 'idx'
pvec <- vector()
for(nr in 1:(ncol(som_tod)-1)){
        pval <- res_sum[[nr]]$`Pr(>F)`[1]
        pvec <- c(pvec, pval)
}
df_tod$pval <- pvec

# Choose row index with significant p value
c_tod <- df_tod %>%
        filter(pval <= 0.05) %>%
        select(idx) %>% unlist() %>% as.numeric()
# Select only genes that show significant variance across timepoints
som_tod <- som_tod[,c_tod]
dim(som_tod)

# Compare Genes Between SOMs (TOD vs Exercsie/Control)
################################
colnames(som_tod) %in% colnames(som_xgroup) %>% table()
venn_file <- paste0(WD,'/plots/20200413-rnaseq-female-kidney-venn-TOD-vs-EC_steep.png')
venn.diagram(
        x = list(colnames(som_tod), colnames(som_xgroup)),
        category.names = c("Time of Death (MANOVA Genes)", 
                           "Exercise Control Group (MANOVA Genes"),
        filename = venn_file,
        imagetype = "png",
        output=TRUE
)
img = readImage(venn_file)
EBImage::display(img, method = "raster")

# Center and scale the data for SOM
proc_tod <- preProcess(som_tod, method=c("center", "scale"))
norm_tod <- predict(proc_tod, som_tod)
summary(norm_tod)
tod_mat <- t(norm_tod) %>% as.matrix()
# SOM
set.seed(666)
som_tod_plot <- som::som(tod_mat,6,5)
plot(som_tod_plot, ylim=c(-2,2))

table(som_tod_plot$visual[,1:2])

# Perform k-means cluster_toding
k_tod <-kmeans(som_mat,6)$cluster

# Create an annotation for the heatmap
ann_df <- col_data %>%
        filter(sample_key %in% colnames(som_mat)) %>%
        select(specimen.collection.t_death_bins.type)
ann_df$specimen.collection.t_death_bins.type <- factor(ann_df$specimen.collection.t_death_bins.type, levels = tod_levels)
ann_df <- ann_df %>% arrange(specimen.collection.t_death_bins.type)
row.names(ann_df) <- colnames(som_mat)

ann_colors = list(
        specimen.collection.t_death_bins.type = 
                c("10.00-11.25" = "gold",
                  "11.75-12.50" = "darkgoldenrod1",
                  "13.25-14.00" = "orange",
                  "14.25-15.00" = "darkorange2",
                  "17.00-18.00" = "darkorange4"))

# Create heatmaps from kmeans
heat_tods <- vector('list', 6)
heat_tods[[1]] <- pheatmap(som_mat[k_tod==1,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)
heat_tods[[2]] <- pheatmap(som_mat[k_tod==2,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)
heat_tods[[3]] <- pheatmap(som_mat[k_tod==3,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)
heat_tods[[4]] <- pheatmap(som_mat[k_tod==4,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)
heat_tods[[5]] <- pheatmap(som_mat[k_tod==5,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)
heat_tods[[6]] <- pheatmap(som_mat[k_tod==6,], 
                           annotation_col = ann_df,
                           annotation_colors = ann_colors,
                           fontsize = 8,
                           show_rownames=FALSE,
                           show_colnames=FALSE,
                           color = rev(brewer.pal(n = 9, name ="RdBu")),
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           legend=F,
                           annotation_legend = FALSE)

# Clusters are numbers in numerical order from 1 to 6 from bottom-left to top-right
# Identify the circadian rhythm genes that belong to each cluster
i <- 1
cluster_tod <- list()
circ_tod <- list()
circ_tod_df <- list()
ENSEMBL_RAT <- vector()
for(yn in 0:4){
        for(xn in 0:5){
                # Collect gene ids for each cluster_tod
                c_rn <- som_tod_plot$visual[(som_tod_plot$visual$x == xn & som_tod_plot$visual$y == yn),] %>%
                        row.names() %>% as.numeric()
                cluster_tod[[i]] <- som_tod_plot$data[c_rn,] %>% row.names()
                # Collect all gene ids
                ENSEMBL_RAT <- c(ENSEMBL_RAT, cluster_tod[[i]])
                # collect circadian gene ids per cluster_tod
                i <- i + 1
        }
}

# Create a geom_tile to mimic the som and visualize the kmeans cluster_tods
grad_tods <- vector('list', 6)
kmeans_cluster_tod <- list()
for(kn in 1:6){
        # Kmeans cluster_tod 1
        x <- 0:5
        y <- 0:4
        tile_df <- expand.grid(X=x,Y=y)
        # Kmeans cluster_tod genes are in which som cluster_tods
        kgenes <- names(k_tod[k_tod==kn])
        kmeans_cluster_tod[[kn]] <- kgenes
        circ_tod[[kn]] <- circ_kid %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(ENSEMBL_RAT) %>% unlist(use.names = FALSE)
        # collect additional information for circadian genes in cluster_tod
        circ_tod_df[[kn]] <- circ_kid %>%
                filter(ENSEMBL_RAT %in% kgenes) %>%
                select(NUM.TISSUE, ENSEMBL_RAT, SYMBOL_RAT)
        Z <- vector()
        for(i in 1:30){
                Z <- c(Z, (kgenes[kgenes %in% cluster_tod[[i]]] %>% length()))
        }
        tile_df$Z <- Z
        # Generate the kmeans plot
        grad_tods[[kn]] <- local({
                kn <- kn
                kmsom <- ggplot(tile_df, aes(x = X, y = Y, fill =Z)) +
                        geom_raster(interpolate=TRUE) +
                        scale_fill_gradient2(low="navy", mid="white", high="red", 
                                             midpoint=30, limits=range(tile_df$Z),
                                             guide = FALSE) +
                        theme(axis.text        = element_blank(),
                              axis.ticks       = element_blank(),
                              axis.title       = element_blank(),
                              panel.background = element_blank())
                print(kmsom)
        })
}

# Create a dataframe for logistic regression
logreg_tod <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
logreg_tod <- mutate(logreg_tod, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 1, 0))
# Generate 6 columns for each kmeans cluster_tod status
for(i in 1:6){
        new_col_name <- paste0('C',i)
        logreg_tod <- logreg_tod %>% 
                mutate(!!sym(new_col_name) := ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[i]], 1, 0))
}

# Remove gene ids and set them to rownames
row.names(logreg_tod) <- logreg_tod$ENSEMBL_RAT
logreg_tod <-  logreg_tod %>% select(-ENSEMBL_RAT)

# Adjust column objects
factor_cols <- colnames(logreg_tod)
for(fc in factor_cols){
        logreg_tod[[fc]] <- as.factor(logreg_tod[[fc]])
}
str(logreg_tod)

# Fit a logistic regression model
mod_tod <- glm(formula = CIRC ~ C1+C2+C3+C4+C5+C6, 
               family = binomial(link = logit), 
               data = logreg_tod)
summary(mod_tod)

# Perform a Chi-square analysis
#######################
# Create a proper data structure for chi squared
x2_list <- list()
for(i in 1:6){
        n_clust <- length(kmeans_cluster_tod[[i]])
        n_circ <- length(circ_tod[[i]])
        n_no_circ <- n_clust - n_circ
        x2_list[[i]] <- c(n_circ, n_no_circ)
}

x2_df <- as.data.frame(x2_list)
row.names(x2_df) <- c('CIRC','NON-CIRC')
for(i in 1:6){
        new_col_name <- paste0('C',i)
        colnames(x2_df)[i] <- new_col_name
}

# Create the proper data structure for visualization of chi square with mosaic plot
mosaic_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
mosaic_df <- mutate(mosaic_df, CIRC = ifelse(ENSEMBL_RAT %in% circ_kid$ENSEMBL_RAT, 'CIRC', "NON-CIRC"))
# Generate a Cluster columcluster_tod status
mosaic_df <- mosaic_df %>%
        mutate(CLUSTER = case_when(ENSEMBL_RAT %in% kmeans_cluster_tod[[1]] ~ 1,
                                   ENSEMBL_RAT %in% kmeans_cluster_tod[[2]] ~ 2,
                                   ENSEMBL_RAT %in% kmeans_cluster_tod[[3]] ~ 3,
                                   ENSEMBL_RAT %in% kmeans_cluster_tod[[4]] ~ 4,
                                   ENSEMBL_RAT %in% kmeans_cluster_tod[[5]] ~ 5,
                                   ENSEMBL_RAT %in% kmeans_cluster_tod[[6]] ~ 6))
# Generate a table of results
mosaic_tbl <- table(mosaic_df$CIRC, mosaic_df$CLUSTER, dnn=c("CIRC","CLUSTER"))

# Perform Chi-squared
#sink('')
(( c=chisq.test(x2_df, simulate.p.value = TRUE) ))
(( fisher.test(x2_df, simulate.p.value = TRUE) ))

#sink()
c$observed
c$expected

# Visualize the chisquare analysis with a mosaic plot
mos_tod <- mosaic(~ CIRC + CLUSTER, data = mosaic_tbl,
                  shade = TRUE, legend = TRUE)
plot(logreg_tod)

# Clusters enriched WITH circadian genes
# 4
grad_tods[[4]]
heat_tods[[4]]
mypar()
plot(som_tod_plot, ylim=c(-2,2))

# Cluster enriched WITHOUT circadian genes
# 5
grad_tods[[5]]
heat_tods[[5]]
plot(som_tod_plot, ylim=c(-2,2))

# Add a column for time of death binned specifically to exercise group
# TODO: Add this to original metadata generation script, then update excel documentation
#col_data <- col_data %>%
#        mutate(specimen.collection.t_death_bins_ec = case_when(
#                specimen.collection.t_death_bins.type == "Exercise - IPE" ~ "13:15 - 13:35",
#                specimen.collection.t_death_bins.type == "Exercise - 0.5 hr" ~ "11:55 - 12:15",
#                specimen.collection.t_death_bins.type == "Exercise - 1 hr" ~ "13:55 - 14:15",
#                specimen.collection.t_death_bins.type == "Exercise - 4 hr" ~ "14:35 - 14:55",
#                specimen.collection.t_death_bins.type == "Exercise - 7 hr" ~ "17:35 - 17:55",
#                specimen.collection.t_death_bins.type == "Exercise - 24 hr" ~ "10:35 - 10:55",
#                specimen.collection.t_death_bins.type == "Exercise - 48 hr" ~ "11:15 - 11:35",
#                specimen.collection.t_death_bins.type == "Control - IPE" ~ "9:55 - 10:15",
#                specimen.collection.t_death_bins.type == "Control - 7 hr" ~ "16:55 - 17:15"))


# Plot the circadian genes in cluster_tod 4
# Create a dataframe for the plot
df1 <- data.frame(t(som_mat), check.names = FALSE)
df1$sample_key <- row.names(df1)
df2 <- col_data %>%
        select(specimen.collection.t_death, 
               specimen.collection.t_death_bins.type)
df2$sample_key <- row.names(df2)
# Adjust the column names to be symbols for ease of plot interpretation
df_plot <- left_join(df1, df2, by = "sample_key")
genes <- colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))]
symbols <- mapIds(org.Rn.eg.db, genes, "SYMBOL", "ENSEMBL")
colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))] <- symbols

# Melt the plot
melt_plot <- reshape2::melt(df_plot, id.vars = c("specimen.collection.t_death", 
                                                 "specimen.collection.t_death_bins.type"))
# Adjust columns and factor levels
melt_plot$value <- as.numeric(melt_plot$value)
melt_plot$specimen.collection.t_death_bins.type <- factor(melt_plot$specimen.collection.t_death_bins.type, 
                                                          levels = tod_levels)

# To re-examine SOM
plot(som_tod_plot, ylim=c(-2,2))

# Plot the circadian genes for cluster_tod of choice
# Show all tods they are enriched with and without in powerpoint
plot(grad_tods[[4]])
heat_tods[[4]]
melt_plot %>%
        filter(variable %in% circ_tod_df[[4]]$SYMBOL_RAT) %>%
        ggplot(aes(x = as.integer(specimen.collection.t_death_bins.type), 
                   y = value, 
                   color = variable),
               group = "1") +
        geom_point(alpha = 0.05) +
        #geom_line(aes(group = "1")) +
        stat_smooth(alpha = 0.02, se = F) +
        ylab("rlog Transformed Expression \n(centered [0] and scaled [1])") +
        xlab("Time of Death") +
        scale_x_continuous(breaks = seq(1,5,by=1),
                           labels=tod_levels) +
        theme(legend.position = "none")


# TODO: Turn this into a function rat2mouse
################################################################################
x <- row.names(rld)
# Load in annotations
mart_mm_ens = useMart("ensembl", dataset="mmusculus_gene_ensembl")
mart_rn_ens = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
# Create ortholog table
ortho_df <- getLDS(attributes=c("ensembl_gene_id","mmusculus_homolog_orthology_confidence"),
                   filters="ensembl_gene_id", 
                   values = x, 
                   mart=mart_rn_ens,
                   attributesL=c("ensembl_gene_id"), 
                   martL=mart_mm_ens) # Use biomart to get orthologs
str(ortho_df)

# Filter out any low confidence orthologs and any genes that are not one-to-one orthologs in both directions
ortho_df <- ortho_df[ortho_df$Mouse.orthology.confidence..0.low..1.high. == '1',]
ortho_df <- ortho_df[!duplicated(ortho_df[,1]),]
ortho_df <- ortho_df[!duplicated(ortho_df[,3]),]
names(ortho_df) <- c('ENSEMBL_RAT','CONFIDENCE','ENSEMBL_MOUSE') 
ortho_df <- ortho_df %>%
        select(-CONFIDENCE)
################################################################################

# Pathway enrichment for all genes, CIRC genes, and NON-CIRC genes in cluster_tod 5
all_out <- kmeans_cluster_tod[[4]]
circ_out <- circ_tod[[4]]

# Convert genes to mouse orthologs for pathway enrichment
all_out <- ortho_df %>%
        filter(ENSEMBL_RAT %in% all_out) %>%
        select(ENSEMBL_MOUSE) %>% unlist() %>% as.character()
all_out <- mapIds(org.Mm.eg.db, all_out, "SYMBOL", "ENSEMBL") %>% unlist(use.names = FALSE)

# Save these list to file
out_file <- paste0(WD,'/data/20200420_rnaseq-kidney-kmeans-tod-clust4-all_steep.txt')
write.table(all_out, file = out_file, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep ='\n')

# Convert genes to mouse orthologs for pathway enrichment
circ_out <- ortho_df %>%
        filter(ENSEMBL_RAT %in% circ_out) %>%
        select(ENSEMBL_MOUSE) %>% unlist() %>% as.character()
circ_out <- mapIds(org.Mm.eg.db, circ_out, "SYMBOL", "ENSEMBL") %>% unlist(use.names = FALSE)

# Save these list to file
out_file <- paste0(WD,'/data/20200420_rnaseq-kidney-kmeans-tod-clust4-circ_steep.txt')
write.table(circ_out, file = out_file, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep ='\n')

# Pathway enrichment for all genes, CIRC genes, and NON-CIRC genes in any cluster_tod that resembles any sort of circadian effect



kmeans_cluster_tod[[4]]


# Create the proper data structure
somvsom_df <- data.frame(ENSEMBL_RAT)
# Generate a column for circadian gene status
somvsom_df <- mutate(somvsom_df, TOD1 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[1]], 1, 0))
somvsom_df <- mutate(somvsom_df, TOD2 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[2]], 1, 0))
somvsom_df <- mutate(somvsom_df, TOD3 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[3]], 1, 0))
somvsom_df <- mutate(somvsom_df, TOD4 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[4]], 1, 0))
somvsom_df <- mutate(somvsom_df, TOD5 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[5]], 1, 0))
somvsom_df <- mutate(somvsom_df, TOD6 = ifelse(ENSEMBL_RAT %in% kmeans_cluster_tod[[6]], 1, 0))
# Generate a column for circadian gene status
somvsom_df <- mutate(somvsom_df, EC1 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[1]], 1, 0))
somvsom_df <- mutate(somvsom_df, EC2 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[2]], 1, 0))
somvsom_df <- mutate(somvsom_df, EC3 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[3]], 1, 0))
somvsom_df <- mutate(somvsom_df, EC4 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[4]], 1, 0))
somvsom_df <- mutate(somvsom_df, EC5 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[5]], 1, 0))
somvsom_df <- mutate(somvsom_df, EC6 = ifelse(ENSEMBL_RAT %in% kmeans_cluster[[6]], 1, 0))
somvsom_cor <- somvsom_df %>% select(-ENSEMBL_RAT)




plot(cor(somvsom_cor))
str(somvsom_df)
cor(somvsom_df)

# Generate a table of results
mosaic_tbl <- table(somvsom_df$EC1, somvsom_df$TOD1, dnn=c("EC","TOD"))

#



























#' #### End of Script Notes

#+ End of Script Notes
################################################################################
################ End of Script Notes ###########################################
################################################################################

# TODO: Fix the hclust call to rld.500, either create rld.500 first or move hclust call down
# TODO: Adjust and double check all interpretations before publishing as html

# Examine character columns
col_data[, sapply(col_data, class) == 'character']
# Examine factor columns
col_data[, sapply(col_data, class) == 'factor']
# Numeric columns
col_data[, sapply(col_data, class) == 'numeric']
# Integer columns
col_data[, sapply(col_data, class) == 'integer']
# Logical columns
col_data[, sapply(col_data, class) == 'logical']

# Create a SOM (with kohonen)
set.seed(666)
# Specify grid parameters
g <- somgrid(xdim = 10, ydim = 10, topo = "rectangular")
# Specify map (alpha is the learning rate--default: 0.05 & 0.01)
komap <- kohonen::som(som_tod,
                      grid = g,
                      alpha = c(0.05, 0.01),
                      radius = 1)

# Visualize this map
# Plot the average distance to the closest unit based on iterations: tells us how many iterations we need
plot(komap, type = 'changes')
mypar(1,1)
# For a codes plot
plot(komap, type = 'codes')
# To plot a counts plot
plot(komap, type = 'count')
mypar()
# To plot elements (features) in each node plot
plot(map, type = 'mapping')
# To plot a neighbor distance plot, light colors indicate bigger distance
plot(map, type = 'dist.neighbours')



# Steps for investigating categorical variables:
# * Examine categorical variable in PCA
# * Create contingenyc table and use Chi square to determine if categorical varibale are associated with eachother

# Outlier Detection
################################################################################
#DESeq2::plotPCA(vstd, intgroup ="animal.registration.sex") +
#        guides(color=guide_legend(title="Sex"))

#par(mar=c(8,5,2,2))
#boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#assays(dds)[["cooks"]]

#https://support.bioconductor.org/p/35918/
################################################################################

#' #### Session Information
################################################################################
####################### Session Info ###########################################
################################################################################
session_info()

# Stops evaluation of code
#knitr::opts_chunk$set(eval = FALSE)

