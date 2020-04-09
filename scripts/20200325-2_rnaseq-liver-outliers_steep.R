# TODO: Remove an outlier
mypar()
pcaData <- DESeq2::plotPCA(rld, intgroup=c("animal.registration.sex"), 
                           returnData=TRUE, ntop = nrow(assay(rld)))
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color =animal.registration.sex)) +
        geom_point(size=3) +
        #geom_text_repel(aes(label=sample_key),hjust=0, vjust=0) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        #coord_fixed() +
        ggtitle("Male Liver Samples Sequenced ")

"To gain an overview of sample heterogeneity, we calculated sample-sample similarities for each region using pairwise Pearson’s correlation coefficients (r) and calculated the average r of each sample compared with all other samples of the same region. We chose the threshold of average r = 0.85–0.94 (varying by region) to define and remove outlier microarrays. The outliers could result from either technical or biological differences."

# Outlier detection and repetative samples

# Hierarchical clustering based on euclidean distance