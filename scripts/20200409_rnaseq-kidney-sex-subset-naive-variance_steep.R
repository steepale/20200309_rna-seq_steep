#'---
#' title: "PASS1A Rat Kidney: Sex Subsets -- Groupings by Naive Models"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "20200409"
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
#' TODO: Edit these ...
#' * Determine rat orthologs of circadian genes in mouse kidney
#' * Perform PCAs with kidney sex subsets, do they cluster together by:
#'     * Time of day?
#'     * Exercise cohort?
#'     * Other metadata?
#' * Subset circadian rhythm genes to determine and generate more PCAs.
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
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn","org.Mm.eg.db","VennDiagram","EBImage","reshape2","xtable","kohonen")
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
        filter(animal.registration.sex == 'Female') %>%
        select(sample_key) %>% unlist() %>% as.character()
fkid_counts <- count_data[,fkid] %>% as.matrix()

# All samples unique?
#col_data %>%
#  filter(Tissue == 'Kidney') %>%
#  select(vial_label) %>% unlist() %>% unique() %>% length()

# SUbset kidney metadata
fkid_cols <- col_data %>%
        filter(Tissue == 'Kidney') %>%
        filter(animal.registration.sex == 'Female')
#filter(sample_key != '90042016803_SF1')
row.names(fkid_cols) <- fkid_cols$sample_key

#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))

design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
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

#' ##### Note: Number of genes with average counts between zero and 1 is `r zero_n` but removing reads with less than or equal to `r reads_n` removes an additional `r filter_n` features or removes `r filter_p*100`% of the non-zero reads (total of `r total_n` reads removed).
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

#' Regularized Log (rlog) Transform
for( n in 1){
        start_time <- Sys.time()
        rld <- rlog(dds)
        end_time <- Sys.time()
        print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
#' #### Here we visualize the counts comparing log2 transform and rlog transform
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] +1), cex=0.2, main = "log2 Transform") # Log2 Transform
boxplot(assay(rld)[rs > 0,], cex=0.2, main = "rlog Transform") # Regularized Log Transform

#' ##### Returning to our prior plot of comparing expression of 2 samples, we can see the the variance is now more stable in the rlog transform.
mypar(1,1)
plot(assay(rld)[,c(1,2)],
     xlab = "Sample 1 rlog Gene Counts", 
     ylab = "Sample 2 rlog Gene Counts", 
     cex = 0.3, main = "rlog Transform")

#' TODO: QC Check, Jun recommends demonstrating the number of features at different (perhaps dec) -tiles from different transforms. We need to double check that the transform did not remove features or shrink low expression features so extremely to demonstrate zero variance.

#' ##### Another plot to examine the difference between transformations is a mean standard deviation plot, which calculates--for each gene--the mean over all samples and the standard deviation over all samples..
mypar(1,1)
#' #### rlog Transform
meanSdPlot(assay(rld)[rs > 0,], ranks=FALSE)

#' ## PCA and Hierarchical Clustering (Female Kidneys) 

#+ PCA and Hierarchical Clustering (Female Kidneys)
################################################################################
##########    PCA and Hierarchical Clustering (Female Kidneys)       ###########
################################################################################

#' TODO: With all PCA's consider creating ordinal values for groups and generating more intuitive colors

#' #### We decide to subset by sex and investigate how samples cluster annotated by exercise/control groups (first we show major 3 bin types)
#' Top 500 genes ranked by total variance
pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.key.anirandgroup.bins.1"), returnData=TRUE, 
                           ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup.bins.1)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples ") +
        guides(color=guide_legend(title="Control/Exercise Groups"))

#' #### How female kidney samples cluster annotated with exercise/control groups
#' Note that the Control - 7 hr samples cluster a bit further away (PC2) than expected (compared to liver)
pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.key.anirandgroup"), returnData=TRUE, 
                           ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples ") +
        guides(color=guide_legend(title="Control/Exercise Groups"))

# Collect the top 500 genes that demonstrate the most variance
rv <- apply(assay(rld), 1, var)
topgenes <- head(order(rv, decreasing = TRUE), 500)
# Note: prcomp expect samples to be rows
pca <- prcomp(t(assay(rld)[topgenes,]), center = TRUE, scale. = FALSE)
#Doublecheck the plot
#autoplot(pca)
rv_500 <- head(rv[order(rv, decreasing = TRUE)], 500)
# Capture the top 500 genes driving variance in the Naive PCA
top500 <- row.names(assay(rld)[topgenes,])
rld.500 <- rld[row.names(rld) %in% top500,]

#' #### An (unrooted) euclidean distance dendrogram of the top 500 genes driving variance demonstrates further supports the notion that immediate exercise groups show relatively higher variance in their gene expressions. The dendrogram adds precision to the groupings (compared to and in support of PCA).
hc <- hclust(dist(t(assay(rld.500))))
# Create a dendrogram with different annotations
myplclust(hc, labels=colData(rld.500)[["animal.key.anirandgroup"]],
          #cex=0.8, 
          lab.col=as.fumeric(as.character(colData(rld.500)[["animal.key.anirandgroup"]])), 
          main="Euclidean Distance Dendrogram: \nExercise/Control Groups")
# Including binned time of death annotation demonstrates the troubling experimental design but it may distract from the following findings; therefore it is not included in the final report but code is provided.
#myplclust(hc, labels=colData(rld.sub)[["specimen.collection.t_death_bins.type"]],
#cex=0.8, 
#          lab.col=as.fumeric(as.character(colData(rld.sub)[["specimen.collection.t_death_bins.type"]])), 
#          main="Euclidean Distance Dendrogram: \nTime of Death (binned)")

#' #### How female kidney samples cluster annotated with time of death bins
#' top 500 genes ranked by total variance
pcaData <- DESeq2::plotPCA(rld, intgroup=c("specimen.collection.t_death_bins.type"), returnData=TRUE, 
                           ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=specimen.collection.t_death_bins.type)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples ") +
        guides(color=guide_legend(title="Time of Death\n(Binned by Hour Intervals)"))

#' ##### The variances associated with PC1 is quite significant (40% variance).
summary(pca)$importance[,1:10]
plot((summary(pca)$importance[,1:10])[2,]*100, ylab = '% Variance', xlab = "PC")

#' ##### A Venn Diagram comparing the top 500 genes driving variance in female kidney samples shows more overlap between orthologous circadian genes in mouse kidney.
#' Circadian genes in top 500 genes driving variance in 
circ_kid_ens %in% top500 %>% table()
# TODO: Consider creating a Chi-square or fisher test
venn_file <- paste0(WD,'/plots/20200404-rnaseq-female-kidney-venn-500-vs-circ_steep.png')
venn.diagram(
        x = list(top500, circ_kid_ens),
        category.names = c("   --------Top 500 Variance Genes", 
                           "Rat Orthologs of \nCircadian Genes \nin Mouse Kidney"),
        filename = venn_file,
        imagetype = "png",
        output=TRUE
)
img = readImage(venn_file)
#TODO: Create a much more visually appealing venn diagram
EBImage::display(img, method = "raster")
shared_genes <- mapIds(org.Rn.eg.db, circ_kid_ens[circ_kid_ens %in% top500], "SYMBOL", "ENSEMBL")
#' ##### The shared high variance orthologous genes:
shared_genes

#' #### When we examine the variance associated with circadian genes across the entire dataset (201 circadian genes expressed), we see that circadian genes might actually be a driver of variance.
mypar()
rld.circ <- rld[row.names(rld) %in% circ_kid_ens,]
pcaData <- DESeq2::plotPCA(rld.circ, intgroup=c("animal.key.anirandgroup.bins.1"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup.bins.1)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples \n(Circadian Kidney Genes Only)") +
        guides(color=guide_legend(title="Control/Exercise Groups"))
# Time of death
pcaData <- DESeq2::plotPCA(rld.circ, 
                           intgroup=c("specimen.collection.t_death_bins.type"), 
                           returnData=TRUE, 
                           ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=specimen.collection.t_death_bins.type)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples ") +
        guides(color=guide_legend(title="Time of Death\n(Binned by Hour Intervals)"))

#' #### Instead, if we select 201 genes at random, we do not see the same pattern in variance. Note the low PC1.
#' TODO: Critic from Jun: When we select genes by random, we should downsample. Genes should have comparible total and average expression levels, then place in PCA.
set.seed(666)
mypar()
rld.sub <- rld[row.names(rld) %in% sample(row.names(rld),201),]
pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.key.anirandgroup.bins.1"), 
                           returnData=TRUE, ntop = 500)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup.bins.1)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Female Kidney Samples: \nRandomly Sampled Genes (n = 201)") +
        guides(color=guide_legend(title="Control/Exercise Groups"))

# Collect the variance for top 500 genes and for circadian rhythm genes
mypar()
set.seed(666)
rv_top_noncirc <- sample(rv_500,201)
rv_circ <- apply(assay(rld.circ), 1, var)
rv_random <- sample(rv,201)
df1 <- data.frame("VARIANCE" = rv_circ, "SET" = "CIRCADIAN", row.names = NULL)
df2 <- data.frame("VARIANCE" = rv_500, "SET" = "TOP_500", row.names = NULL)
var_df <- rbind(df1,df2)

# TODO: Edit this description:
#' ##### A histogram of the variance per gene demonstrates that Circadian genes do not capture much of the variance at all, they are drowned out by the top 500 genes. 
#' Note: This histogram is reduced on the x-axis to 1. To better put things into perspective, the sum variance from all circadian genes was `r sum(rv_circ)` while the sum variance from the top 500 genes was `r sum(rv_500)` and the maximum variance from any gene was `r max(rv_500)`.
ggplot(var_df, aes(x=VARIANCE, color=SET, fill=SET)) +
        geom_histogram(aes(y=..density..), position="identity", alpha=0.5) +
        labs(title="Variance Histogram: \nTop 500 Genes vs Circadian Rhythm Genes",x="Variance (per gene)", y = "Frequency")+
        xlim(0,1) +
        theme_classic()

# Add the randomly selected genes to the dataframe for the boxplot
df3 <- data.frame("VARIANCE" = rv_random, "SET" = "RANDOM_201", row.names = NULL)
var_df <- rbind(df1,df2,df3)
#' #### A boxplot with supporting summary stats show a similar story, except that 201 genes selected at random capture nealry identical variance as circadian rhythm genes. Note: Much of y axis has been chopped.
# Boxplot
ggplot(var_df, aes(x=SET,y=VARIANCE, color=SET)) +
        geom_boxplot(show.legend = FALSE) +
        ylim(0,0.5) +
        xlab("") +
        ylab("Variance (per gene)") +
        ggtitle("Variance Boxplot: \nPer Gene Variance of Different Gene Sets")
#' #### Top 500:
summary(rv_500)
#' #### Circadian:
summary(rv_circ)
#' #### 201 Randomly Selected Genes:
summary(rv_random)

#' ### Circadian Genes capture variance in Female Kidney Samples under a naive model (~ 1)
#' #### Model Expression of Circadian Genes
#+
################################################################################
######## Model Expression of Circadian Genes ###################################
################################################################################

# Create a dataframe for the plot
df1 <- data.frame(t(assay(rld.circ)))
df1$sample_key<- row.names(df1)
df2 <- col_data %>%
        select(specimen.collection.t_death, animal.key.anirandgroup.bins.1)
df2$sample_key <- row.names(df2)
# Adjust the column names to be symbols for ease of plot interpretation
df_plot <- left_join(df1, df2, by = "sample_key")
genes <- colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))]
symbols <- mapIds(org.Rn.eg.db, genes, "SYMBOL", "ENSEMBL")
colnames(df_plot)[grepl('ENSRNOG', colnames(df_plot))] <- symbols

# Melt the plot
melt_plot <- melt(df_plot, id.vars = c("specimen.collection.t_death", 
                                       "animal.key.anirandgroup.bins.1"))
# Change expression value to numeric
melt_plot$value <- as.numeric(melt_plot$value)

# TODO: Adjust this interpretation
#' #### When we examine the expression of a set of circadian genes selected at random, we see that they do not considerably chnage their expression over time
set.seed(666)
rando_plot <- melt_plot %>%
        filter(variable %in% sample(melt_plot$variable, 10))
ggplot(rando_plot, 
       aes(x = specimen.collection.t_death, 
           y = value, 
           color = variable)) +
        geom_point(alpha = 0.5) +
        stat_smooth(alpha = 0.5) +
        ylab("rlog Transformed Expression") +
        xlab("Specimen collection time of death")

#' #### Even when we examine the expression of circadian genes with the largest variance, we see that they do not considerably chnage their expression over time
top_plot <- melt_plot %>%
        filter(variable %in% shared_genes[1:10])
ggplot(top_plot, 
       aes(x = specimen.collection.t_death, 
           y = value, 
           color = variable)) +
        geom_point(alpha = 0.5) +
        stat_smooth(alpha = 0.5) +
        ylab("rlog Transformed Expression") +
        xlab("Specimen collection time of death")

#' ##### TODO: Calculate maximum log2 fold chnage of circadian genes and show that is doesn't equate to noise in differential expression


#' #### Differential Gene Expression Analysis

#+ Differential Gene Expression Analysis
################################################################################
######## Differential Gene Expression Analysis #############################
################################################################################

Control0
Ex0
Ex0.5
Ex1
Ex4
Ex7
Ex24
Ex48

Control7
Ex0
Ex0.5
Ex1
Ex4
Ex7
Ex24
Ex48

#' #### Unsupervised Clustering w/ SOMs & K-means

#+ Unsupervised Clustering w/ SOMs & K-means
################################################################################
######## Unsupervised Clustering w/ SOMs & K-means #############################
################################################################################

# Order samples by exercise group
group_order <- fkid_cols %>%
        arrange(animal.key.anirandgroup) %>%
        select(sample_key) %>% unlist() %>% as.character()
group_order2 <- fkid_cols %>%
        arrange(animal.key.anirandgroup) %>%
        select(animal.key.anirandgroup) %>% unlist() %>% as.character()

# Arrange matrix in order of exercise and control groups
# All expression values per gene are normalized to generate a level playing field. Means are subtracted by values and the resulting difference is divided by the standard deviation.
som_xgroup <- assay(rld)[,group_order]

som_xgroup <- som_xgroup %>% t() %>% data.frame()


# Test for multiv
x <- as.numeric(MVN::mvn(som_xgroup, mvnTest = "hz")$univariateNormality$`p value`)
(p.adjust(x, method = "BH") > 0.05) %>% table()

som_xgroup$animal.key.anirandgroup <- group_order2

# MANOVA test
colnames(som_xgroup)
res.man <- manova(cbind(ENSRNOG00000059776,
                        ENSRNOG00000020843,
                        ENSRNOG00000056940,
                        ENSRNOG00000012366) ~ animal.key.anirandgroup, data = som_xgroup)

gi <- c("ENSRNOG00000059776",
        "ENSRNOG00000020843",
        "ENSRNOG00000056940",
        "ENSRNOG00000012366")
si <- mapIds(org.Rn.eg.db, gi, "SYMBOL", "ENSEMBL")
tp <- melt_plot %>%
        filter(variable %in% si)

ggplot(tp, 
       aes(x = specimen.collection.t_death, 
           y = value, 
           color = variable)) +
        geom_point(alpha = 0.5) +
        stat_smooth(alpha = 0.5) +
        ylab("rlog Transformed Expression") +
        xlab("Specimen collection time of death")

res.man


res.man <- manova(cbind(Sepal.Length, Petal.Length) ~ Species, data = iris)
summary(res.man)

summary.aov(res.man)

# Order samples by time of death (time of day)
tod_order <- fkid_cols %>%
        arrange(specimen.collection.t_death) %>%
        select(sample_key) %>% unlist() %>% as.character()
# Arrange matrix in order of time of death
reads_n <- 10
keep <- rowSums(fkid_counts)/ncol(fkid_counts) >= reads_n
fkid_counts_10 <- fkid_counts[keep,]
fkid_counts_10 %>% scale()
som_tod <- assay(rld.circ)[,tod_order] %>% scale()

# Data should be normalized differently between

# SOM
set.seed(666)
som::som(som_xgroup,10,10) %>% plot()
som::som(som_tod,4,4) %>% plot()


# SOM
set.seed(666)
# Specify grid parameters
g <- somgrid(xdim = 10, ydim = 10, topo = "rectangular")
# Specify map (alpha is the learning rate--default: 0.05 & 0.01)
map <- som(som_tod,
           grid = g,
           alpha = c(0.05, 0.01),
           radius = 1)




# Visualize this map
# Plot the average distance to the closest unit based on iterations: tells us how many iterations we need
plot(map, type = 'changes')
mypar(1,1)
# For a codes plot
plot(map, type = 'codes')
# To plot a counts plot
plot(map, type = 'count')
mypar()
# To plot elements (features) in each node plot
plot(map, type = 'mapping')
# To plot a neighbor distance plot, light colors indicate bigger distance
plot(map, type = 'dist.neighbours')

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

