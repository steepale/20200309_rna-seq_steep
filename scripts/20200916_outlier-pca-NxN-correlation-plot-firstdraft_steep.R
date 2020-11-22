#'---
#' title: "Chapter 3: Data Pre-Processing"
#' author: Alec Steep
#' date: "`r format.Date( Sys.Date(), '%Y%m%d' )`"
#' output: 
#'     html_document:
#'         code_folding: hide
#'         toc: true
#'         highlight: zenburn
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(cache = FALSE)

#' ## Setup the Environment

#+ Setup Environment

################################################################################
##### Resources and Dependencies ###############################################
################################################################################

# Set the working directory
#WD <- '/Users/Alec/Documents/applied_predictive_modeling'
#setwd(WD)

# Load the dependencies
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("ff")
#install.packages("AppliedPredictiveModeling")

# Load dependencies
pacs...man <- c("tidyverse","data.table","R.utils","ggpubr","plyr","ROCR","ff","GenomicRanges","BSgenome","rtracklayer", "caret","pROC","modelr","ggplot2","AppliedPredictiveModeling","e1071")
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

############################################################
##### Functions ############################################
############################################################


# Function to speed up making rows into lists for interation with lapply
################################################################################
f_pmap_aslist <- function(df) {
        purrr::pmap(as.list(df), list)
}
################################################################################

################################################################################
#####     Reference      #######################################################
################################################################################

# Determine where the in-depth scripts are located
AppliedPredictiveModeling::scriptLocation()

#' ## Load & Clean Data
#+ Load the Data

################################################################################
#####     Load & Clean Data      ###############################################
################################################################################

# The raw segmentation data
########################
data("segmentationOriginal")

df <- as_tibble(segmentationOriginal)

# Columns:
# Cell: Each cell
# Class: Which cells were segmented
# Case: Which cells were used for training and test set

# Filter data to contain only training data
train_df <- df %>% filter(Case == 'Train')

# Save the Class and Cell fileds into vectors and remove thme from the tibble
cellID <- train_df$Cell
class <- train_df$Class
case <- train_df$Case
# Remove these fileds from dataset
train_df <- train_df %>% dplyr::select(-Cell,-Class,-Case)

# The original data contained "status" columns that contained binary versions of some predictors. Find them with "status" in the name and remove them.
status_col <- grep("Status", names(train_df))
train_df <- train_df[, -status_col]

################################################################################
#####     Transformations      #################################################
################################################################################

PC_cor

# Load the metadata information
rnaseq_meta_annotation_file <- '/Volumes/Frishman_4TB/motrpac/20200309_rna-seq_steep/merged/merged_column_dictionary2019-10-13.txt'
meta_col_df <- read.table(file = rnaseq_meta_annotation_file, header = T, sep = '\t')
meta_col_df <- meta_col_df %>% select(Field.Name, Data.Type) %>%
        mutate(Data.Type = ifelse(Data.Type == "varchar(MAX)", "varchar", Data.Type)) %>%
        mutate(Data.Type = ifelse(Data.Type == "smallint", 'int', Data.Type))

#table(meta_col_df$Data.Type)
# Collect the numeric coulmns
# quant_cols <- meta_col_df %>%
#         filter(Data.Type %in% c('int', 'numeric')) %>%
#         select(Field.Name) %>% unlist() %>% as.character()
# qual_cols <- meta_col_df %>%
#         filter(Data.Type %in% c('varchar', 'datetime')) %>%
#         select(Field.Name) %>% unlist() %>% as.character()
# Columns to add
nums <- unlist(lapply(PC_cor, is.numeric))
quant_cols <- names(PC_cor[, nums]) %>% unique()
#quant_cols <- quant_cols[quant_cols %!in% quant_cols[grepl("PC", quant_cols)]]
qual_cols <- names(PC_cor[, !nums]) %>% unique()

# # Columns of interest (numeric)
# coi_num <- c('PC1', 'PC2', 'RNA_extr_conc..ng.uL.', 'RIN','r_260_280','r_260_230',
#              'Lib_RNA_conc..ng.ul.','Lib_RNA_vol..ul.','Lib_UMI_cycle_num..bp.',
#              'Lib_adapter_size..bp.','Lib_frag_size..bp.','Lib_DNA_conc..ng.uL.',
#              'Lib_molarity..nM.','RNA_extr_conc','Lib_RNA_conc','Lib_RNA_vol',
#              'Lib_frag_size', 'Lib_adapter_size', 'Lib_DNA_conc','Lib_molarity',
#              'acute.test.days_visit')
# 
# # Columns of interest (categorical)
# coi_cat <- c('GET_site','RNA_extr_plate_ID','RNA_extr_date','Lib_prep_date',
#              'Lib_robot','Lib_vendor','Lib_type','Lib_kit_id','Lib_batch_ID',
#              'Lib_barcode_well','Lib_index_1','Lib_index_2','Lib_adapter_1',
#              'Lib_adapter_2','Seq_platform','Seq_date','Seq_machine_ID',
#              'Seq_flowcell_ID','Seq_flowcell_run','Seq_flowcell_lane',
#              'Seq_flowcell_type','Seq_length','Seq_end_type','Lib_UMI_cycle_num',
#              'acute.test.d_visit', '')


################################################################################
####### Process the Metadata ###################################################
################################################################################

# -1. Qualitative Columns to keep
################################################################################
qual_keep <- c('GET_site','RNA_extr_plate_ID','RNA_extr_date','Lib_prep_date','Lib_robot',
               'Lib_kit_id','Lib_batch_ID','Seq_date','Seq_machine_ID','Seq_flowcell_ID',
               'Seq_flowcell_run','Seq_flowcell_lane','acute.test.d_visit',
               'acute.test.d_start','acute.test.t_start','acute.test.contactshock',
               'acute.test.comments','animal.familiarization.d_visit',
               'animal.familiarization.d_treadmillbegin','animal.familiarization.d_treadmillcomplete',
               'animal.key.batch','animal.key.intervention',
               'animal.key.anirandgroup','animal.registration.d_visit','animal.registration.ratid',
               'animal.registration.d_arrive','animal.registration.batchnumber',
               'animal.registration.d_reverselight','animal.registration.d_birth',
               'animal.registration.sex','animal.registration.cagenumber',
               'specimen.collection.d_visit','specimen.collection.t_anesthesia',
               'specimen.collection.bloodcomplete','specimen.collection.t_bloodstart',
               'specimen.collection.t_bloodstop','specimen.collection.t_edtafill',
               'specimen.collection.t_death','specimen.collection.deathtype','bid',
               'specimen.processing.samplenumber','specimen.processing.t_collection',
               'specimen.processing.t_freeze','shiptositeid','animal.key.is_control',
               'Seq_batch','animal.key.exvsctrl','animal.key.exlt4',
               'animal.key.anirandgroup.bins.1','specimen.collection.t_death_bins.type')

# 0. Create a list of columns that are removed in process
################################################################################
# as.numeric(PC_cor2$acute.test.d_start)
# length(qual_cols)
# for(i in 141:150){
#         print(paste0('COLUMN ',i))
#         print(qual_cols[i])
#         print(table(as.character(PC_cor2[,qual_cols[i]])))
#         print(table(as.character(col_data[,qual_cols[i]])))
# }
# summary(col_data$animal.key.intervention)
# table(col_data$animal.familiarization.versionguid)


# Reason 1: Columns are not of interest, all value unique
remove_1 <- names(col_data)[apply(col_data, 2, function(x) length(unique(x)) == nrow(col_data))]

# Reason 2: All values are NA
remove_2 <- names(PC_cor2)[apply(PC_cor2, 2, function(x)all(is.na(x)))]

# Reason 3: Manual investigation
remove_3 <- c('sample_key','vial_label','X2D_barcode','Species','BID','PID','Tissue',
             'Sample_category','Lib_barcode_well','Lib_index_1','Lib_index_2',
             'labelid','pid','acute.test.participantguid','acute.test.formname',
             'animal.familiarization.participantguid','animal.familiarization.formname',
             'Lib_UMI_cycle_num','animal.familiarization.fat','animal.familiarization.lean',
             'animal.familiarization.comments','animal.key.participantguid','animal.key.protocol',
             'animal.key.agegroup','animal.key.sitename','animal.registration.participantguid',
             'animal.registration.staffid','animal.registration.siteguid',
             'animal.registration.siteid','animal.registration.formguid','acute.test.formguid',
             'acute.test.siteid','Lib_vendor','Lib_type','Lib_adapter_1','Lib_adapter_2',
             'Seq_platform','Seq_flowcell_type','Seq_end_type','acute.test.staffid',
             'acute.test.siteguid','acute.test.versionguid','acute.test.versionnbr',
             'animal.familiarization.staffid','animal.familiarization.siteguid',
             'animal.familiarization.siteid','animal.familiarization.versionguid',
             'animal.familiarization.versionnbr','animal.familiarization.compliant',
             'animal.registration.staffid','animal.registration.siteguid',
             'animal.registration.siteid','animal.registration.formguid',
             'animal.registration.formname','animal.registration.versionguid',
             'animal.registration.versionnbr','animal.registration.comments',
             'specimen.collection.participantguid','specimen.collection.siteguid',
             'specimen.collection.formguid','specimen.collection.formname',
             'specimen.collection.versionguid','specimen.collection.anesthesiacomments',
             'specimen.collection.bloodtype','specimen.collection.bloodtube',
             'specimen.collection.bloodtechid','specimen.collection.bloodcomments',
             'specimen.collection.uterustype','specimen.collection.uterustechid',
             'specimen.collection.uteruscomplete','specimen.collection.t_uterusstart',
             'specimen.collection.t_uterusstop','specimen.collection.uterusweight',
             'specimen.collection.uteruscomments','specimen.processing.formguid',
             'specimen.processing.versionguid','specimen.processing.versionnbr',
             'specimen.processing.formname','specimen.processing.siteguid',
             'specimen.processing.siteid','specimen.processing.participantguid',
             'specimen.processing.labelguid','specimen.processing.sampletypedescription',
             'specimen.processing.sampletypeguid','specimen.processing.aliquotdescription',
             'specimen.processing.aliquotguid','specimen.processing.volume',
             'specimen.processing.partialamt','specimen.processing.hemolyzed',
             'specimen.processing.comments','specimen.processing.t_edtaspin',
             'specimen.processing.techid','barcode','receivedcas','receivestatuscas')

# Determine which columns ahve no variance
remove_4 <- names(which(apply(PC_cor2, 2, var) == 0))
remove_4 <- c(remove_4, 'acute.test.intenseinitial')

PC_cor2 <- PC_cor %>% select(-all_of(c(remove_1,remove_2,remove_3,remove_4)))

# 2. Process columns with NA values
################################################################################
# Find columns with missing values
NA_cols <- colnames(PC_cor2)[ apply(PC_cor2, 2, anyNA) ]
# Most NAs associated with controls (confounded) or RNA quality measures

# 3. Asjust the date & time columns to be both numeric and leave the originals to be one hot encoded
################################################################################
# Determine date and time columns that should be converted to both numeric and categorical
date_to_num <- c("acute.test.d_start","animal.familiarization.d_treadmillbegin",
                 "animal.registration.d_visit","animal.registration.d_reverselight",
                 "specimen.collection.d_visit","acute.test.d_visit","animal.familiarization.d_visit",
                 "animal.familiarization.d_treadmillcomplete","animal.registration.d_arrive",
                 "animal.registration.d_birth")
time_to_num <- c("acute.test.t_complete","specimen.collection.t_bloodstart",
                 "specimen.collection.t_edtafill","specimen.processing.t_collection",
                 "specimen.processing.t_freeze","calculated.variables.time_death_to_collect_min",
                 "acute.test.t_start","specimen.collection.t_anesthesia",
                 "specimen.collection.t_bloodstop","specimen.collection.t_death",
                 "calculated.variables.time_collect_to_freeze_min","animal_time_last_fed")
# Create numeric versions of all date and time columns
for(C in c(date_to_num,time_to_num)){
        ( C_num <- as.symbol(paste0(C,'_numeric')) )
        PC_cor2[[sym(C_num)]] <- as.numeric(PC_cor2[[sym(C)]])
}

# 4. Impute additional numeric columns to have Gaussian distributions
################################################################################
nums <- unlist(lapply(PC_cor2, is.numeric))
quant_cols <- names(PC_cor2[, nums]) %>% unique()

# Examine each predictor and determine if it demonstrates skewness. If so, correct for skewness via transformation. (e1071 package)

# All columns are numeric so we can apply the skewness function to all columns
skew_vals <- lapply(PC_cor2[,quant_cols], e1071::skewness) %>% unlist()

r_skew_names <- names(skew_vals)[skew_vals > 1]
r_skewed <- r_skew_names[!is.na(r_skew_names)]
mypar(2,2)
for(i in 1:4){
        hist(PC_cor2[[r_skewed[i]]], breaks = 50)
}
l_skew_names <- names(skew_vals)[skew_vals < -1]
l_skewed <- l_skew_names[!is.na(l_skew_names)]
mypar(2,2)
for(i in 1:4){
        hist(PC_cor2[[l_skewed[i]]], breaks = 50)
}

# The BoxCoxTrans function in the caret package can be used to determine which sort of transformation to use (cite Box & Cox, 1964).
trans <- preProcess(PC_cor2[, quant_cols],
                    method = c("BoxCox", "center", "scale"))
# Visualize the results
trans
# Apply the transformations:
PC_cor3 <- predict(trans, PC_cor2)
# We can revisualize the columns
mypar(2,2)
for(i in 1:4){
        hist(PC_cor3[[r_skewed[i]]], breaks = 50)
}
for(i in 1:4){
        hist(PC_cor3[[l_skewed[i]]], breaks = 50)
}

# We can add suffixes to columns
colnames(PC_cor3) <- paste(colnames(PC_cor3), "bxcs", sep = "_")

# Add columns 
PC_cor4 <- cbind(PC_cor2, PC_cor3)

# 5. One hot encode the qualitative variables (dummify)
################################################################################
# Update the columns we have choosen to keep (we need to filter out some zerovar columns)
qual_keep <- qual_keep[qual_keep %in% names(PC_cor4)]
# Convert all quantitative variables to factors
qual_keep2 <- c()
for(C in qual_keep){
        PC_cor4[[C]] <- as.character(PC_cor4[[C]])
        PC_cor4[[C]] <- as.factor(PC_cor4[[C]])
        if(nlevels(PC_cor4[[C]]) >= 2){
                qual_keep2 <- c(qual_keep2,C)
        }
}
# Dummify the qualitative variables with more than 2 levels
dmy <- dummyVars(" ~ .", data = PC_cor4[, qual_keep2])
PC_cor5 <- data.frame(predict(dmy, newdata = PC_cor4[, qual_keep2]))
PC_cor4 <- PC_cor4 %>% select(-all_of(qual_keep))
PC_cor4 <- cbind(PC_cor4, PC_cor5)

# Remove columns without at least 2 unique values
PC_cor4 <- PC_cor4[, sapply(PC_cor4, function(col) length(unique(col))) > 1]

################################################################################
####### Perform Linear Models on the Processed Metadata ########################
################################################################################

# Get the names of the PC variables
VARS <- names(PC_cor4)[names(PC_cor4) %!in% names(PC_cor4)[grepl('PC',names(PC_cor4))]]

PC_mod <- data.frame()
for(PC in c('PC1','PC2','PC1_bxcs','PC2_bxcs')){
        for(VAR in VARS){
                design <- as.formula(paste0(PC,' ~ ',VAR))
                mod_pc1 <- lm(design, data = PC_cor4)
                if(!is.na(summary(mod_pc1)$adj.r.squared) & summary(mod_pc1)$adj.r.squared > 0.3){
                        PC_append <- data.frame('Var1' = PC, 'Var2' = VAR, 
                                                'R2' = summary(mod_pc1)$adj.r.squared)
                        PC_mod <- rbind(PC_mod, PC_append)
                }
        }
}

PC_mod$Var1 <- factor(PC_mod$Var1, levels = c("PC2","PC2_bxcs","PC1","PC1_bxcs"))

################################################################################
####### Plot a heatmap of the Correlated Values ########################
################################################################################

min_corr <- min(PC_mod$R2)
max_corr <- max(PC_mod$R2)
# Set the breaks
brks <- round(seq(min_corr, max_corr, length.out=10), digits = 3)
# Set the colors
col <- colorRampPalette(c("orange","red"))((max_corr-min_corr)*1000)
# Plot the heatmaps
quantmap_list <- list()
quantmap_list[[i]] <- PC_mod %>%
        arrange(Var1) %>%
        ggplot(aes(x = Var2,
                   y = Var1,
                   fill = R2)) +
        geom_tile(color = "gray") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                         size = 5, hjust = 1)) +
        theme(axis.text.y = element_text(size = 5)) +
        theme(axis.title.x=element_blank(),
              axis.title.y=element_blank()) +
        #guides(fill = guide_legend(override.aes = list(size = 1))) +
        theme(legend.title = element_text(size = 4),
              legend.text = element_text(size = 4)) +
        theme(legend.key.size = unit(0.5, 'lines')) +
        theme(legend.key.width = unit(0.5, 'lines')) +
        scale_fill_gradientn(colors = col,
                             name="Correlation",
                             limits=c(min_corr,1),
                             breaks = brks) +
        ggtitle(paste0(TIS,'-',ZeroGenes,'-',TForm,'-',Corr,'-',TForder,'-',BatchOrder)) +
        theme(plot.title = element_text(size=6)) +
        coord_fixed() +
        theme(aspect.ratio=1)

PC_mod %>% arrange(desc(R2)) %>% filter(Var1 %in% c('PC1', 'PC1_bxcs')) %>% head(n=8)





# Melt the correlation matrix
melted_quantmat <- melt(quant_mat, na.rm = TRUE) %>%
        filter(Var1 %in% c('PC1','PC2','PC1_')) %>%
        filter(Var2 %!in% c('PC1','PC2','PC3','PC4','PC5')) %>%
        mutate(Var1 = factor(Var1, levels = c('PC2','PC1'))) %>%
        arrange(rev(Var1))
melted_quantmat$Tissue <- TISSUE
melted_quantmat$Transform <- TForm
melted_quantmat$Correlation <- Corr
melted_quantmat$Gene_Filter <- ZeroGenes
melted_quantmat$Filter_Order <- TForder
melted_quantmat$value <- round(melted_quantmat$value, digits = 3)











