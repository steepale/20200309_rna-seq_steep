
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
plot(grad_tods[[5]])
heat_tods[[5]]
melt_plot %>%
        filter(variable %in% circ_tod_df[[5]]$SYMBOL_RAT) %>%
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
all_out <- kmeans_cluster_tod[[1]]
circ_out <- circ_tod[[1]]

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

# Pathway enrichment for all genes, CIRC genes, and NON-CIRC genes in any cluster_tod that resembles any sort of circadian effect
