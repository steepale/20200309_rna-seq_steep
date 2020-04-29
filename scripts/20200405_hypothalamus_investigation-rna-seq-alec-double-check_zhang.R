#'---
#' title: "PASS1A Rat Hypothalamus: -- Correct for Batch -- (Alec's) Double Check of Jiayu's Script"
#' author: "Jiayu Zhang"
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
auth <- "jiayu"
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


#' #### Retrieve Circadian Genes Associated with Tissue (Kidney)
#' Data from Supplementary Table 2 from 1. Yan, J., Wang, H., Liu, Y. & Shao, C. Analysis of gene regulatory networks in the mammalian circadian rhythm. PLoS Comput. Biol. 4, (2008).
#' Downloaded 20200326 by Alec Steep
#' Data previosuly saved in script 20200413_rnaseq-kidney-female-temporal-clusters_steep.R

# Circadian Genes (Kidney)
in_file <- paste0(WD,'/data/20200409_rnaseq-circadian-kidney-mouse-rat-ortho_steep-yan.txt')
circ_kid <- read.table(file=in_file, sep = '\t', header = TRUE)

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

#' ## Subset Data for Samples of Interest

#+ Subset Data for Samples of Interest
################################################################################
###########     Subset Data for Samples of Interest   ##########################
################################################################################

#' ##### Female Kidney Samples (39 unique) from 1 batch were filtered prior to normalization
fkid <- col_data %>%
        filter(Tissue == 'Hypothalamus') %>%
        filter(!is.na(animal.registration.sex)) %>%
        select(sample_key) %>% unlist() %>% as.character()
fkid_counts <- count_data[,fkid] %>% as.matrix()

# SUbset kidney metadata
fkid_cols <- col_data %>%
        filter(Tissue == 'Hypothalamus') %>%
        filter(!is.na(animal.registration.sex))

#filter(sample_key != '90042016803_SF1')
row.names(fkid_cols) <- fkid_cols$sample_key
#' #### Sanity Check: Ensure that the metadata rownames are identical to count matrix column names
all(rownames(fkid_cols) == colnames(fkid_counts))
# Adjust labels within columns (e.g. 1 and 2 to female and male)
# see ../phenotype/merged/20191013_merged-column-dictionary_steep-edits.xlsx
########################


### PCA Visualization

count_data = fkid_counts
col_data = fkid_cols

design = ~1
title = paste0('Design: ',as.character(design))
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = design)


# Filter genes with average count of 10 or less. (this time I won't perform this)
# Reasoning from:
#citation("PROPER")
#dds
# keep <- rowMeans(counts(dds)) >= 10
# dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
rs <- rowSums(counts(dds))
# Normalize the counts
rld <- vst(dds) #vst and rlog comparable with all samples
#rld <- rlog(dds, blind=FALSE)

# Extract matrix of normalized counts
counts <- assay(rld)

# Annotate normalized counts
counts$ensembl <- rownames(counts)
counts$symbol <- mapIds(org.Rn.eg.db, counts$ensembl, "SYMBOL", "ENSEMBL")
counts$entrez <- mapIds(org.Rn.eg.db, counts$ensembl, "ENTREZID", "ENSEMBL")
counts$genename <- mapIds(org.Rn.eg.db, counts$ensembl, "GENENAME", "ENSEMBL")
counts$go <- mapIds(org.Rn.eg.db, counts$ensembl, "GO", "ENSEMBL")
counts$path <- mapIds(org.Rn.eg.db, counts$ensembl, "PATH", "ENSEMBL")

# Here we examine all samples by batch and tissue
# 
# first PCA with intGroup of tissue
DESeq2::plotPCA(rld, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="sex"))


pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.key.anirandgroup","animal.registration.sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup, shape = animal.registration.sex)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Hypothalamus Samples Across Groups")


# Variables of interest
male_hypothalamus <- (col_data %>% 
                              filter(animal.registration.sex == 'Male'))$vial_label %>% 
        as.character()
female_hypothalamus <- (col_data %>% 
                                filter(animal.registration.sex == 'Female'))$vial_label %>% 
        as.character()
ref_hypothalamus <- (col_data %>%
                             filter(is.na(animal.registration.sex)))$vial_label %>% 
        as.character()
hypothalamus <- c(male_hypothalamus,female_hypothalamus,ref_hypothalamus)
Y_genes <- Y_ens_id[Y_ens_id %in% row.names(assay(rld))]
X_genes <- X_ens_id[X_ens_id %in% row.names(assay(rld))]
sex <- col_data[hypothalamus,"animal.registration.sex"]
group <- col_data[hypothalamus,"animal.key.anirandgroup"]


# t-test on expression associated with sex
tt = rowttests(counts[, hypothalamus], sex)
par(myfrow = c(1,2))
# Histogram of p values associated with ttest
hist(tt$p.value,main="",ylim=c(0,1300), breaks = 100)
plot(tt$dm,-log10(tt$p.value))
points(tt[X_genes,]$dm,-log10(tt[X_genes,]$p.value),col=3,pch=16)
points(tt[Y_genes,]$dm,-log10(tt[Y_genes,]$p.value),col=2,pch=16, xlab="Effect size",ylab="-log10(p-value)")
legend("bottomright",c("X","Y"),col=3:2,pch=16)
p <- tt$p.value
qvals <- qvalue(tt$p.value)$qvalue
index <- which(qvals<=0.05)
abline(h=-log10(max(tt$p.value[index])))
```

### PCA after removing sex-differentially expressed genes 
#### (have I removed sex-related vairance?)
```{r}
# first try removing all sex-significant genes
hypothalamus_counts_rm_sig = hypothalamus_counts[!(rownames(hypothalamus_counts)) %in% rownames(hypothalamus_counts[index,hypothalamus]),]

# then try removing only X or Y related sex-significant genes
# here's a new index
sex_index = c(index[which(rownames(hypothalamus_counts[index,hypothalamus]) %in% X_genes)], index[which(rownames(hypothalamus_counts[index,hypothalamus]) %in% Y_genes)])

hypothalamus_counts_rm_sex = hypothalamus_counts[!(rownames(hypothalamus_counts)) %in% rownames(hypothalamus_counts[sex_index,hypothalamus]),]

# build dds with those 2 new count matrix
hypo_dds_rm_sig <- DESeqDataSetFromMatrix(countData = hypothalamus_counts_rm_sig,
                                          colData = hypothalamus_cols,
                                          design = design)

# remove genes that have avg counts less than 1 (reduced to maintain more sex-relate genes)
keep <- rowSums(counts(hypo_dds_rm_sig))/ncol(hypo_dds_rm_sig) > 1
hypo_dds_rm_sig <- hypo_dds_rm_sig[keep,]
hypo_dds_rm_sig = estimateSizeFactors(hypo_dds_rm_sig)
rlog_norm_hypo_dds_rm_sig = rlog(hypo_dds_rm_sig)


DESeq2::plotPCA(rlog_norm_hypo_dds_rm_sig, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="sex"))


pcaData <- DESeq2::plotPCA(rlog_norm_hypo_dds_rm_sig, intgroup=c("animal.key.anirandgroup","animal.registration.sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup, shape = animal.registration.sex)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Hypothalamus Samples With 22 Genes Removed")
```
```{r}
hypo_dds_rm_sex <- DESeqDataSetFromMatrix(countData = hypothalamus_counts_rm_sex,
                                          colData = hypothalamus_cols,
                                          design = design)

# remove genes that have avg counts less than 1 (reduced to maintain more sex-relate genes)
keep <- rowSums(counts(hypo_dds_rm_sex))/ncol(hypo_dds_rm_sex) > 1
hypo_dds_rm_sex <- hypo_dds_rm_sex[keep,]
hypo_dds_rm_sex = estimateSizeFactors(hypo_dds_rm_sex)
rlog_norm_hypo_dds_rm_sex = rlog(hypo_dds_rm_sex)


DESeq2::plotPCA(rlog_norm_hypo_dds_rm_sex, intgroup ="animal.registration.sex") +
        guides(color=guide_legend(title="sex"))


pcaData <- DESeq2::plotPCA(rlog_norm_hypo_dds_rm_sex, intgroup=c("animal.key.anirandgroup","animal.registration.sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.registration.sex)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Hypothalamus Samples with 16 Genes Removed")
```

## Defferential gene expression
### Modeling raw counts with normalization

now perform *differential gene expression* on the counts, try to find genes with dfference in counts across groups due to the condition of interest rises above the biological and technical variance we observe

first prepare for col data and count data
```{r}
# variables in design formula cannot contain NA, so I have to remove ref samples here
hypothalamus_counts_rm_sex = hypothalamus_counts_rm_sex[,-which(is.na(col_data$animal.key.anirandgroup))]
hypothalamus_cols = hypothalamus_cols[-which(is.na(col_data$animal.key.anirandgroup)),]

# group all samples by acute exercise or not c(0,1)
for (i in 1:nrow(hypothalamus_cols)){
        if (hypothalamus_cols[i,]$animal.key.anirandgroup %in% c("Exercise - 0.5 hr", "Exercise - 1 hr", "Exercise - IPE")){
                hypothalamus_cols[i,]$acute_ex = 1
        }else{
                hypothalamus_cols[i,]$acute_ex = 0
        }
}
hypothalamus_cols$acute_ex = as.factor(hypothalamus_cols$acute_ex)

# one more step to remove batch effect between male and female
f_counts = hypothalamus_counts_rm_sex[, female_hypothalamus[ which(female_hypothalamus %in% colnames(hypothalamus_counts_rm_sex))]]
m_counts = hypothalamus_counts_rm_sex[, male_hypothalamus[ which(male_hypothalamus %in% colnames(hypothalamus_counts_rm_sex))]]
f_median = rowMedians(f_counts)
m_median = rowMedians(m_counts)
f_var = matrixStats::rowVars(f_counts)
m_var = matrixStats::rowVars(m_counts)

# for lots of genes, there's 0 expression in all samples which results in a variance of 0
# they will be removed after building dds object, so I'll just keep them here

new_var = pmin(f_var, m_var)
f_var_scale = ifelse(f_var == 0, 0, sqrt(new_var / f_var))
m_var_scale = ifelse(m_var == 0, 0, sqrt(new_var / m_var))
for (i in 1:nrow(hypothalamus_counts_rm_sex)){
        f_counts[i] = f_counts[i] + max(f_median[i], m_median[i]) - f_median[i]
        m_counts[i] = m_counts[i] + max(f_median[i], m_median[i]) - m_median[i]
        # m_counts[i] = ifelse(m_var_scale == 0, m_counts[i], m_var_scale[i] * m_counts[i] - (m_var_scale[i] - 1) * mean(m_counts[i])) %>% as.integer()
        # f_counts[i] = ifelse(f_var_scale == 0, f_counts[i], f_var_scale[i] * f_counts[i] - (f_var_scale[i] - 1) * mean(f_counts[i])) %>% as.integer()
}

hypothalamus_counts_scaled = cbind(f_counts, m_counts)
hypothalamus_cols_scaled = hypothalamus_cols[c(female_hypothalamus[ which(female_hypothalamus %in% colnames(hypothalamus_counts_rm_sex))], male_hypothalamus[ which(male_hypothalamus %in% colnames(hypothalamus_counts_rm_sex))]),]
```

```{r}
# Create a DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = hypothalamus_counts_scaled,
                              colData = hypothalamus_cols_scaled,
                              design = ~ animal.registration.sex + acute_ex)
# based on the vignette of DESeq2: acute_ex (the last variable in design) effect represents the overall effect controlling for differences due to sex

# Filter genes with average count of 10 or less. (this time I won't perform this)
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowMeans(counts(dds)) >= 1
dds <- dds[keep,]

dds_1 <- estimateSizeFactors(dds)
# Normalize the counts
rld <- rlog(dds_1) #vst and rlog comparable with all samples
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

pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.key.anirandgroup","animal.registration.sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup, shape = animal.registration.sex)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Hypothalamus Samples with sex median-centered") 
```
### to Compare whether I have removed the sex-related variance
* subset by sex == female
* PCA show the data structure

```{r}
hypothalamus_cols_f = hypothalamus_cols[female_hypothalamus,]
hypothalamus_counts_f = hypothalamus_counts[,female_hypothalamus]
# Create a DESeqDataSet Object
dds_f <- DESeqDataSetFromMatrix(countData = hypothalamus_counts_f,
                                colData = hypothalamus_cols_f,
                                design = ~ 1)


# Filter genes with average count of 10 or less. (this time I won't perform this)
# Reasoning from:
#citation("PROPER")
#dds
keep <- rowMeans(counts(dds_f)) >= 1
dds_f <- dds_f[keep,]

dds_f <- estimateSizeFactors(dds_f)
# Normalize the counts
rld <- rlog(dds_f) #vst and rlog comparable 

pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.key.anirandgroup","acute_ex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=animal.key.anirandgroup, shape = acute_ex)) +
        geom_point(size=3) +
        #geom_text(aes(label=Tissue),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Hypothalamus Samples (Female only)")
```






Then run the DESeq2 model:
        ```{r}
dds <- DESeq2::DESeq(dds)
res = results(dds)
head(res)

```
**for myself reference**
        Another way to look at the difference is that a p-value of 0.05 implies that 5% of all tests  
will result in false positives. An FDR adjusted p-value (or q-value) of 0.05 implies that 5%   
of significant tests will result in false positives. The latter will result in fewer false  
positives.  
```{r}
# table(res$padj < 0.1)
# summary(res)
#which()
res2 <- results(dds, alpha=0.05)
table(res2$padj < 0.05)

# a sorted results table:
resSort <- res2[order(res2$padj),]

gene_ls = counts[rownames(res2[which(res2$padj <= 0.05),]),]$symbol
gene_ls = gene_ls[!is.na(gene_ls)]
write.table(gene_ls, file = paste0(WD, "/Data/intermediate/20200423-hypothalamus-DE-gene-list_Zhang.txt"), sep = "\t", col.names = F, row.names = F, quote = FALSE)

```
### Visualizing results
```{r}
plotMA(res2, ylim=c(-4,4), alpha = 0.05)
hist(res2$pvalue[res2$baseMean > 1], 
     col="grey", border="white", xlab="pvalue", ylab="frequency", main="")
```
Examine the counts for the top gene, sorting by FDR:
        ```{r}
#plotCounts(dds, gene=which.min(res2$padj), intgroup="animal.key.anirandgroup")
data <- plotCounts(dds, gene=which.min(res2$padj), intgroup="animal.key.anirandgroup", returnData=TRUE)
ggplot(data, aes(x=animal.key.anirandgroup, y=count, col = hypothalamus_cols$hour_death)) +
        geom_point(position=position_jitter(width=.1,height=0)) +
        #stat_smooth() + 
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 45))

#which gene it is
counts[which.min(res2$padj),c("symbol","genename")]
```
Draw heatmap of the top genes:
        
        ```{r}
topgenes <- head(rownames(resSort),20)
mat <- assay(rld)[topgenes,]

# standarized the count matrix
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)$animal.key.anirandgroup)
rownames(df) = colnames(mat)
pheatmap(mat, annotation_col=df)
```
Generally, exercise-IPE and exercise-0.5hr are clustered together, apart from other groups.  


As we talked in meeting, the control-7 group and exercise-7 group are well matched  
So I try to find genes that are significantly differently expressed in those 2 groups
```{r}


res3 <- results(dds, alpha=0.05, contrast = c("animal.key.anirandgroup", "Control - 7 hr", "Exercise - 7 hr"))
table(res3$padj < 0.05)
hist(res3$pvalue[res3$baseMean > 1], 
     col="grey", border="white", xlab="pvalue", ylab="frequency", main="")
data <- plotCounts(dds, gene=which.min(res3$padj), intgroup="animal.key.anirandgroup", returnData=TRUE)
ggplot(data, aes(x=animal.key.anirandgroup, y=count, col = col_data$hour_death)) +
        geom_point(position=position_jitter(width=.1,height=0)) +
        #scale_y_log10() +
        theme(axis.text.x = element_text(angle = 45))

counts[which.min(res3$padj),c("symbol","genename")]
```


Also draw a heatmap for the top DE genes in the control-7 group vs the exercise-7 group:
        ```{r}
resSort_7 <- res3[order(res3$padj),]
topgenes_7 <- head(rownames(resSort_7),20)
mat <- assay(rld)[topgenes_7,]

# standarized the count matrix
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=df)
```

