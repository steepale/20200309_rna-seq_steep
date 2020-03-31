#'---
#' title: "PASS1A Rat Liver (Stanford Batch 1) RNA-Seq: Circadian vs. Exercise: What Drives Variance?"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output:
#'     github_document: default
#'     pdf_document: default
#'     
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Goals of Analysis
#' * TODO: Continue Step-by-step goals of analysis
#' * Perform a PCA of Liver
#'     * If Liver clusters by sex, investigate clustering by autosomal and sex genes
#'         * Histogram of p values by autosomal genes and sex genes (BOTH X and Y chromosomes!)
#'         * Match with volcano plots
#'     * Either correct for sex batch, remove sex genes, or subset one sex.
#' * Perform another PCA with exercise groups, do they cluster together by:
#'     * Time of day?
#'     * Exercise cohort?
#'     * Other metadata? -- EDA (dunno yet)
#' * Build a correlation matrix -- EDA (dunno yet)
#' * Subset genes into 2 groups:
#'     * Genes associated with circadian rhythm in Rat Liver (We don't have this set, we use Mouse Liver)
#'     * Genes associated with exercise effects (TODO: Gene sets need to be collected).
#' * Investiagte batch
#'     * Examine overlap of 2 gene sets
#'     * Rank genes by p-value (dunno how to do this yet, unless by differential expression)
#'         * TODO: Determien ranking strategy
#'     * Statistical test involving rank -- dunno: wilcoxen rank sum test?
#'         * TODO: Determine proper test
#'     * Cluster: PCAs -- EDA
#' * To be continued ...
#' 
#'     
#' 
#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
WD <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("qvalue")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) })

############################################################
##### Functions ############################################
############################################################

# Set select
select <- dplyr::select

# Make the 'not in' operator
################################################################################
'%!in%' <- function(x,y) {
        !('%in%'(x,y))
}
################################################################################

# Capture the Date and AUthor
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
in_file <- paste0(WD,'/data/20200326_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv')
count_data <- read.table(in_file,sep = ',', header = TRUE,row.names = 1,check.names = FALSE)

# Meatdata table
in_file <- paste0(WD,'/data/20200326_rnaseq-meta-pass1a-stanford-sinai_steep.txt')
col_data <- read.table(in_file, header = TRUE, check.names = FALSE, sep = '\t')
row.names(col_data) <- col_data$sample_key

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

# TODO: Load Reference genome and annotations

### Determine which control samples are male and female
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
seqlevels(Rn_TxDb)[1:35]
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

################################################################################
#' ## PCA of Liver Samples 
#' Samples from Stanford Batch 1, which we suspect demonstrates a batch effect

#+ PCA of Liver Samples

################################################################################
#####     PCA of Liver Samples       ###########################################
################################################################################

design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowSums(counts(dds))/ncol(dds) >= 10
dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet
dds

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
norm_counts <- assay(rld)
counts <- as.data.frame(norm_counts)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")

#' #### Liver samples cluster by sex.
#' Grey samples represent "reference samples".
#' 
rld.sub <- rld[ , (rld$Tissue == "Liver") ]
DESeq2::plotPCA(rld.sub, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="Sex"))

# Variables of interest
male_livers <- (col_data %>% 
                        filter(Tissue == 'Liver') %>% 
                        filter(animal.registration.sex == 'Male'))$sample_key %>% 
        as.character()
female_livers <- (col_data %>% 
                          filter(Tissue == 'Liver') %>% 
                          filter(animal.registration.sex == 'Female'))$sample_key %>% 
        as.character()
ref_livers <- (col_data %>% 
                       filter(Tissue == 'Liver') %>% 
                       filter(is.na(animal.registration.sex)))$sample_key %>% 
        as.character()
livers <- c(male_livers,female_livers,ref_livers)
Y_genes <- Y_ens_id[Y_ens_id %in% row.names(norm_counts)]
X_genes <- X_ens_id[X_ens_id %in% row.names(norm_counts)]
sex <- col_data[livers,"animal.registration.sex"]
group <- col_data[livers,"animal.key.anirandgroup"]
liver_counts <- norm_counts[,livers]

#' #### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.
chryexp <- colMeans(norm_counts[Y_genes,livers])

#' ##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution. However, distinct peaks are not detected. This was a surprising result. 
mypar()
hist(chryexp, breaks = 200)
summary(chryexp)
#' We will not use this common strategy to determine sex of unknown samples, rather we will use clustering from PCA.

#' The distribution of sex by group
table(group, sex) # <- Bad idea.

#' #### Generate a heatmap of expression from 3 sets of genes:
#' * Genes from the Y chromosome
#' * The top and bottom 25 genes (50 total) associated with sex
#' * Randomly selected genes

#' ##### Males and Females demonstrate distinctly different gene expression profiles.
#' * Genes on the Y chromosome are not a good predictor of sex in Liver mRNA measures (Figure 1; suprising result)
#' * Male and female samples show distinct correlation to one another (Figure 2)

# T-test of expression associated with sex
tt <- rowttests(liver_counts,sex)

# Take genes from the Y chromosome
# Y_genes
# Take the top and bottom 25 genes associated with variable of interest (remove any genes in Y chromosome)
top <- row.names(tt[order(-tt$dm),][1:25,])
bot <- row.names(tt[order(tt$dm),][1:25,])
top_n_bot <- setdiff(c(top,bot), Y_genes)

# Randomly select 50 genes not in prior sets 
set.seed(123)
randos <- setdiff(row.names(tt[sample(seq(along=tt$dm),50),]), c(Y_genes,top_n_bot))
geneindex <- c(randos,top_n_bot,Y_genes)
# Generate the heatmap and support with a plot of a correlation matrix
mat <- liver_counts[geneindex,]
mat <- mat -rowMeans(mat)
icolors <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))(100)
mypar(1,2)
image(t(mat),xaxt="n",yaxt="n",col=icolors)
y <- liver_counts - rowMeans(liver_counts)
image(1:ncol(y),1:ncol(y),cor(y),col=icolors,zlim=c(-1,1),
      xaxt="n",xlab="",yaxt="n",ylab="")
axis(2,1:ncol(y),sex,las=2)
axis(1,1:ncol(y),sex,las=2)

#' #### A naive t-test and genes with q values less than or equal to 0.05
#' ##### The left figure represents a histogram of p values from a naive t-test (all genes). We see that a significant proportion of genes correlate with sex. To investigate if these genes are located on sex chromosomes or autosomal chromosomes, we contruct a volcano plot on the right. To our suprise, again, genes on the Y chromosome do not show signifcant correlation to sex. Rather some genes on the X chromosome demonstrate significance, however, these genes do not make up the majority of significantly assocaited genes.
#' ## The variance in gene expression is not dicated by differential expression of genes on sex chromosomes.
mypar(1,2)
# Histogram of p values associated with ttest
hist(tt$p.value,main="",ylim=c(0,1300), breaks = 100)
plot(tt$dm,-log10(tt$p.value))
points(tt[X_genes,]$dm,-log10(tt[X_genes,]$p.value),col=1,pch=16)
points(tt[Y_genes,]$dm,-log10(tt[Y_genes,]$p.value),col=2,pch=16, xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("X","Y"),col=1:2,pch=16)
p <- tt$p.value
qvals <- qvalue(tt$p.value)$qvalue
index <- which(qvals<=0.05)
abline(h=-log10(max(tt$p.value[index])))















load(file = paste0(WD,"/data/GSE5859Subset.rda"))
library(rafalib)
library(RColorBrewer)
library(genefilter)

load(file = paste0(WD,"/data/GSE5859.rda"))








pcaData <- DESeq2::plotPCA(rld.sub, intgroup=c("animal.registration.sex","Seq_batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex, shape=Seq_batch)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Gastrocnemius Samples Sequenced Across Batches")





