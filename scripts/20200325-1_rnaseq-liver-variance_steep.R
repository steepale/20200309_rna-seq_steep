#'---
#' title: "PASS1A Rat Liver (Stanford Batch 1): What Drives Variance?"
#' author: "Alec Steep and Jiayu Zhang" 
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
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
#' * TODO: Continue Step-by-step goals of analysis
#' * Perform a PCA of Liver
#'     * If Liver clusters by sex, investigate clustering by autosomal and sex genes
#'         * Histogram of p values by autosomal genes and sex genes (BOTH X and Y chromosomes!)
#'         * Match with volcano plots
#'     * Either correct for sex batch, remove sex genes, or subset one sex
#'         * Decision: Subset by sex
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
#BiocManager::install("caret")
#install.packages("tidyverse")

# Load dependencies
pacs...man <- c("tidyverse","GenomicRanges", "DESeq2","devtools","rafalib","GO.db","vsn","hexbin","ggplot2", "GenomicFeatures","Biostrings","BSgenome","AnnotationHub","plyr","dplyr", "org.Rn.eg.db","pheatmap","sva","formula.tools","pathview","biomaRt", "PROPER","SeqGSEA",'purrr','BioInstaller','RColorBrewer','lubridate', "hms","ggpubr", "ggrepel","genefilter","qvalue","ggfortify","som", "vsn")
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
in_file <- paste0(WD,'/data/20200309_rnaseq-countmatrix-pass1a-stanford-sinai_steep.csv')
count_data <- read.table(in_file,sep = ',', header = TRUE,row.names = 1,check.names = FALSE)

# Meatdata table
in_file <- paste0(WD,'/data/20200309_rnaseq-meta-pass1a-stanford-sinai_steep.txt')
col_data <- read.table(in_file, header = TRUE, check.names = FALSE, sep = '\t')
row.names(col_data) <- col_data$sample_key

# Adjust columns
col_data$animal.key.exlt4 <- as.factor(col_data$animal.key.exlt4)

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

#' ## Normalization for Sequencing Depth

#+ Normalization for Sequencing Depth

################################################################################
#####     Normalization for Sequencing Depth      ##############################
################################################################################

#' ##### Liver Samples (80 unique) from 1 batch were filtered prior to normalization
#' One outlier removed: TODO: Document removal of outlier in seperate script
livers <- col_data %>%
  filter(Tissue == 'Liver') %>%
  filter(sample_key != '90042016803_SF1') %>%
  select(sample_key) %>% unlist() %>% as.character()
liver_counts <- count_data[,livers] %>% as.matrix()

# All samples unique?
#col_data %>%
#  filter(Tissue == 'Liver') %>%
#  select(vial_label) %>% unlist() %>% unique() %>% length()

# SUbset liver metadata
liver_cols <- col_data %>%
  filter(Tissue == 'Liver') %>%
  filter(sample_key != '90042016803_SF1')
row.names(liver_cols) <- liver_cols$sample_key

#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(liver_cols) == colnames(liver_counts))

design = ~ 1 # Primary variable needs to be last
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = liver_counts,
                              colData = liver_cols,
                              design = design)

# Perform pre-filtering.
# Filter genes with average count of 10 or less.
# Reasoning from:
#citation("PROPER")
#dds
#' #### We remove genes with an average sequencing depth of 10 or less
#' Before Filtering
dds
keep <- rowSums(counts(dds))/ncol(dds) > 3
dds <- dds[keep,]
#' #### Summary of counts and annotation data in a DESeqDataSet after filtering out genes with low sequencing depth
dds

mypar()
#' #### When we plot counts for genes across 2 samples, we see that genes with high counts will demonstrate high variance in untransformed data.
plot(assay(dds)[,1:2], cex=.3)
#' To see the reads per million for each sample
sort(colSums(assay(dds)))/1e6

# estimateSizeFactors gives us a robust estimate in sequencing depth
dds <- estimateSizeFactors(dds)
#' #### If we plot the total counts for each sample (colSums) vs each sample size factors, we see they are strongly correlated--size factors act as a coefficient for sequencing depth
plot(sizeFactors(dds), colSums(counts(dds)), cex=.3)
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
#' Size facotrs are generally around 1 (scaled) and calculated using the median and are robust to genes with large read counts
sort(sizeFactors(dds))
summary(sizeFactors(dds))

#' How size factor is computed for a sample:
#' Size factors are from the median ratio of samples compared to a pseudo-sample, which is the geometric mean of all samples.
loggeomeans <- rowMeans(log(counts(dds)))
exp(median((log(counts(dds)[,1]) - loggeomeans)[is.finite(loggeomeans)])) == sizeFactors(dds)[1]


#' ### Comparison of transformation techniques

#+ Transformation Techniques
################################################################################
#####     Transformation Techniques      #######################################
################################################################################

#' ##### Log2 normalization (with size factor and pseudo count)"
log_norm_counts <- log2(counts(dds, normalized=TRUE) + 1)
#' ##### Variance stabilizing transformation (VST) uses a subset of genes and divides each column of the counts matrix by our size factors. It also normalizes with respect to library size.
vstd <- vst(dds)
#' Regularized Log (rlog)
#' The purpose of using the "regularized log" transformation is to shrink together the values of the genes that have very low counts ("shrinkage" technique in statistics). Rlog is close to log2 transform.
#Time difference of 2.914409 mins (on Macbook pro)
for( n in 1){
  start_time <- Sys.time()
  rld <- rlog(dds)
  end_time <- Sys.time()
  print(end_time - start_time)
}

# This command is redundent, but included for safety
rs <- rowSums(counts(dds))
#' #### Here we visualize the counts from log2 transformation and from the rlog transform
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] +1), cex=0.2, main = "log2 Transform") # Log2 Transform
boxplot(assay(rld)[rs > 0,], cex=0.2, main = "rlog Transform") # Regularized Log Transform
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] +1), cex=0.2, main = "log2 Transform") # log2 Transform
boxplot(assay(vstd)[rs > 0,], cex=0.2, main = "VST Transform") # VST Transform

#' ##### Returning to our prior plot of comparing expression of 2 samples, we can see the the variance is now more stable in VST and rlog transforms. However, the VST transform demonstrates a lot of variance in genes with low read counts. For that reason, we will proceed with rlog transform.
mypar(1,1)
plot(assay(vstd)[,c(1,2)],
     xlab = "Sample 1 VST Gene Counts", 
     ylab = "Sample 2 VST Gene Counts", 
     cex = 0.3, main = "VST Transform" )
plot(assay(rld)[,c(1,2)],
     xlab = "Sample 1 rlog Gene Counts", 
     ylab = "Sample 2 rlog Gene Counts", 
     cex = 0.3, main = "rlog Transform")

#' ##### Another plot to examine the difference between transformations is a mean standard deviation plot, which calculates--for each gene--the mean over all samples and the standard deviation over all samples. rlog and vst transforms are comparable.
mypar(1,2)
#meanSdPlot(log_norm_counts, ranks=FALSE) 
mds_vst <- meanSdPlot(assay(vstd)[rs > 0,], ranks=FALSE)
mds_vst$gg + ggtitle("VST Transform")
mds_rlog <- meanSdPlot(assay(rld)[rs > 0,], ranks=FALSE)
mds_rlog$gg + ggtitle("rlog Transform")

################################################################################
#' ## PCA of Liver Samples 
#' Samples from Stanford Batch 1, which we suspect demonstrates a strong technical batch effect

#+ PCA of Liver Samples

################################################################################
#####     PCA of Liver Samples       ###########################################
################################################################################

# Extract matrix of normalized counts
#rld_backup <- rld
#rld <- vstd
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
DESeq2::plotPCA(rld, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="Sex"))

# Variables of interest
male_livers <- (liver_cols %>% 
                        filter(animal.registration.sex == 'Male'))$sample_key %>% 
        as.character()
female_livers <- (liver_cols %>% 
                          filter(animal.registration.sex == 'Female'))$sample_key %>% 
        as.character()
ref_livers <- (liver_cols %>%
                       filter(is.na(animal.registration.sex)))$sample_key %>% 
        as.character()
livers <- c(male_livers,female_livers,ref_livers)
Y_genes <- Y_ens_id[Y_ens_id %in% row.names(norm_counts)]
X_genes <- X_ens_id[X_ens_id %in% row.names(norm_counts)]
sex <- liver_cols[livers,"animal.registration.sex"]
group <- liver_cols[livers,"animal.key.anirandgroup"]

#' #### Predict the sex of reference samples (all samples for that matter) by calculating the median expression of genes on the Y chromosome. We should expect a bimodal distribution with males demonstrating significantly higher median expression.
chryexp <- colMeans(norm_counts[Y_genes,])
#chryexp <- colMeans(assay(dds)[Y_genes,])

chryexp_df <- data.frame("counts" = chryexp)
chryexp_df$sample <- row.names(chryexp_df) %>% as.factor()
chryexp_df <- chryexp_df %>%
  mutate(sex = ifelse(sample %in% male_livers, 'Male', 'Female'))

#' ##### If we create a histogram of the median gene expression values on chromosome Y, we should expect to see a bimodal distribution. However, distinct peaks are not detected. This was a surprising result. 
mypar()
hist(chryexp_df$counts, breaks = 200, 
     main = "Histogram: \nY Chromosome Genes across Samples",
     xlab = "Total reads on Y chromosome genes per sample")
ggplot(chryexp_df, aes(counts, colour = sex)) +
  geom_freqpoly(binwidth = 4) +
  xlab("Total reads on Y chromosome genes per sample") +
  ylab("Frequency") +
  ggtitle("Frequency Polygon: \nY Chromosome Genes across Samples")
# summary(chryexp)

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

