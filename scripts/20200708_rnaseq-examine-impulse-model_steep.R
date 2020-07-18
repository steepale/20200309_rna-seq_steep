vignette("ImpulseDE2_Tutorial")

library(ImpulseDE2)
lsSimulatedData <- simulateDataSetImpulseDE2(
        vecTimePointsA   = rep(seq(1,8),3),
        vecTimePointsB   = NULL,
        vecBatchesA      = NULL,
        vecBatchesB      = NULL,
        scaNConst        = 30,
        scaNImp          = 10,
        scaNLin          = 10,
        scaNSig          = 10,
        scaMuBatchEffect = NULL,
        scaSDBatchEffect = NULL,
        dirOutSimulation = NULL)


objectImpulseDE2 <- runImpulseDE2(
        matCountData    = lsSimulatedData$matObservedCounts, 
        dfAnnotation    = lsSimulatedData$dfAnnotation,
        boolCaseCtrl    = FALSE,
        vecConfounders  = NULL,
        scaNProc        = 1 )

ImpulseDE2:::processData
################################################################################
processData <- function (dfAnnotation, matCountData, boolCaseCtrl, vecConfounders, 
          vecDispersionsExternal, vecSizeFactorsExternal) 
{
        checkNull <- function(objectInput, strObjectInput) {
                if (is.null(objectInput)) {
                        stop(paste0("ERROR: ", strObjectInput, " was not given as input."))
                }
        }
        checkNA <- function(objectInput, strObjectInput) {
                if (is.na(objectInput)) {
                        stop(paste0("ERROR: ", strObjectInput, " is NA and needs to be specifief."))
                }
        }
        checkDimMatch <- function(matInput1, matInput2, strMatInput1, 
                                  strMatInput2) {
                if (any(dim(matInput1) != dim(matInput2))) {
                        stop(paste0("ERROR: ", strMatInput1, " does not have the dimensions as ", 
                                    strMatInput2, "."))
                }
        }
        checkElementMatch <- function(vec1, vec2, strVec1, strVec2) {
                if (!any(vec1 == vec2)) {
                        stop(paste0("ERROR: ", strVec1, " do not agree with ", 
                                    strVec2, "."))
                }
        }
        checkNumeric <- function(matInput, strMatInput) {
                if (any(!is.numeric(matInput))) {
                        stop(paste0("ERROR: ", strMatInput, " contains non-numeric elements. ", 
                                    "Requires numeric data."))
                }
        }
        checkProbability <- function(matInput, strMatInput) {
                checkNumeric(matInput, strMatInput)
                if (any(matInput < 0 | matInput > 1 | is.na(matInput))) {
                        stop(paste0("ERROR: ", strMatInput, " contains elements outside", 
                                    " of interval [0,1]."))
                }
        }
        checkCounts <- function(matInput, strMatInput) {
                checkNumeric(matInput, strMatInput)
                if (any(matInput[!is.na(matInput)]%%1 != 0)) {
                        stop(paste0("ERROR: ", strMatInput, " contains non-integer elements.", 
                                    " Requires count data."))
                }
                if (any(!is.finite(matInput[!is.na(matInput)]))) {
                        stop(paste0("ERROR: ", strMatInput, " contains infinite elements.", 
                                    " Requires count data."))
                }
                if (any(matInput[!is.na(matInput)] < 0)) {
                        stop(paste0("ERROR: ", strMatInput, " contains negative elements.", 
                                    " Requires count data."))
                }
        }
        checkData <- function(dfAnnotation, matCountData, boolCaseCtrl, 
                              vecConfounders, vecDispersionsExternal, vecSizeFactorsExternal, 
                              strReportProcessing) {
                strReportProcessing <- "Processing Details:"
                checkNull(dfAnnotation, "dfAnnotation")
                checkNull(matCountData, "matCountData")
                checkNull(boolCaseCtrl, "boolCaseCtrl")
                vecColNamesRequired <- c("Sample", "Condition", "Time", 
                                         vecConfounders)
                if (!all(vecColNamesRequired %in% colnames(dfAnnotation))) {
                        stop(paste0("Could not find column ", 
                                    vecColNamesRequired[!(vecColNamesRequired %in% 
                                        colnames(dfAnnotation))], " in annotation table."))
                }
                if (any(duplicated(dfAnnotation$Sample))) {
                        stop(paste0("ERROR: [Annotation table] ", "Sample names must be unique: Sample(s) ", paste0((dfAnnotation$Sample)[duplicated(dfAnnotation$Sample)], collapse = ","), " is/are duplicated."))
                }
                vecTimepoints <- unique(as.vector(dfAnnotation$Time))
                checkNumeric(dfAnnotation$Time, "dfAnnotation$Time")
                lsConditions <- unique(dfAnnotation$Condition)
                if (!("case" %in% lsConditions)) {
                        stop(paste0("ERROR: Condition \"case\" does", " not occur in annotation table condition column."))
                }
                if (boolCaseCtrl) {
                        if (!("control" %in% lsConditions)) {
                                stop(paste0("ERROR: Condition \"control\" does", 
                                            " not occur in annotation table condition column."))
                        }
                }
                if (!is.null(vecConfounders)) {
                        for (confounder in vecConfounders) {
                                if (length(unique(dfAnnotation[, confounder])) == 
                                    1) {
                                        stop(paste0("ERROR: Model matrix based on confounding ", 
                                                    "variables {", vecConfounders, "} is not full rank: Only one batch", 
                                                    " given for confounder ", confounder, ". Remove from vecConfounders or correct", 
                                                    " dfAnnotation."))
                                }
                        }
                        matModelMatrix <- do.call(cbind, lapply(vecConfounders, 
                                                                function(confounder) {
                                                                        match(dfAnnotation[, confounder], unique(dfAnnotation[, 
                                                                                                                              confounder]))
                                                                }))
                        if (rankMatrix(matModelMatrix)[1] != dim(matModelMatrix)[2]) {
                                stop(paste0("Model matrix based on confounding variables {", 
                                            vecConfounders, "} is not full rank. ", "Correct the confounding variables. ", 
                                            " Note that it is not possible to model", " NESTED confounding variables: ", 
                                            " Any confounding variables cannot", " be a linear combination of the other ", 
                                            " confounding variables."))
                        }
                }
                if (any(!(colnames(matCountData) %in% dfAnnotation$Sample))) {
                        strReportProcessing <- paste0(strReportProcessing, 
                                                      "\n", "WARNING: The column(s) ", paste0(as.character(colnames(matCountData)[!(colnames(matCountData) %in% 
                                                                                                                                            dfAnnotation$Sample)]), collapse = ","), " in the count data table do(es) not occur", 
                                                      " in annotation table and will be ignored.")
                }
                checkNull(rownames(matCountData), "[Rownames of matCountData]")
                checkCounts(matCountData, "matCountData")
                if (is.null(vecDispersionsExternal) & length(vecConfounders) > 
                    1) {
                        stop(paste0("DESeq2 cannot be run in an automated fashion for multiple ", 
                                    "confounding variables. Run DESeq separately and supply ", 
                                    "dispersion parameters through vecDispersionsExternal."))
                }
                if (!is.null(vecDispersionsExternal)) {
                        if (is.null(names(vecDispersionsExternal))) {
                                stop(paste0("ERROR: vecDispersionsExternal was not named.", 
                                            " Name according to rownames of matCountData."))
                        }
                        if (any(!(names(vecDispersionsExternal) %in% rownames(matCountData))) | 
                            any(!(rownames(matCountData) %in% names(vecDispersionsExternal)))) {
                                stop(paste0("ERROR: vecDispersionsExternal supplied but names", 
                                            " do not agree with rownames of matCountData."))
                        }
                        checkNumeric(vecDispersionsExternal, "vecDispersionsExternal")
                        if (any(vecDispersionsExternal <= 0)) {
                                stop(paste0("WARNING: vecDispersionsExternal contains negative", 
                                            " or zero elements which leads.", "Dispersion parameters must be positive."))
                        }
                }
                if (!is.null(vecSizeFactorsExternal)) {
                        if (is.null(names(vecSizeFactorsExternal))) {
                                stop(paste0("ERROR: vecSizeFactorsExternal was not named.", 
                                            " Name according to colnames of matCountData."))
                        }
                        if (any(!(names(vecSizeFactorsExternal) %in% colnames(matCountData))) | 
                            any(!(colnames(matCountData) %in% names(vecSizeFactorsExternal)))) {
                                stop(paste0("ERROR: vecSizeFactorsExternal supplied but", 
                                            " names do not agree with colnames of matCountData."))
                        }
                        checkNumeric(vecSizeFactorsExternal, "vecSizeFactorsExternal")
                        if (any(vecSizeFactorsExternal <= 0)) {
                                stop(paste0("WARNING: vecSizeFactorsExternal contains negative ", 
                                            "or zero elements which leads.", "Size factors must be positive, remove samples if size ", 
                                            "factor is supposed to be zero."))
                        }
                }
                if (boolCaseCtrl) {
                        strReportProcessing <- paste0(strReportProcessing, 
                                                      "\n", "ImpulseDE2 runs in case-ctrl mode.")
                }
                else {
                        strReportProcessing <- paste0(strReportProcessing, 
                                                      "\n", "ImpulseDE2 runs in case-only mode.")
                }
                strReportProcessing <- paste0(strReportProcessing, "\n", 
                                              paste0("Found time points: ", paste(vecTimepoints, 
                                                                                  collapse = ",")))
                for (tp in vecTimepoints) {
                        strReportProcessing <- paste0(strReportProcessing, 
                                                      "\n", paste0("Case: Found the samples at time point ", 
                                                                   tp, ": ", paste0(dfAnnotation[(dfAnnotation$Time %in% 
                                                                                                          tp) & (dfAnnotation$Condition %in% "case") & 
                                                                                                         (dfAnnotation$Sample %in% colnames(matCountData)), 
                                                                                                 ]$Sample, collapse = ","), collapse = ","))
                }
                if (boolCaseCtrl) {
                        for (tp in vecTimepoints) {
                                strReportProcessing <- paste0(strReportProcessing, 
                                                              "\n", paste0("Control: Found the following samples at time point ", 
                                                                           tp, ":", paste0(dfAnnotation[dfAnnotation$Time %in% 
                                                                                                                tp & dfAnnotation$Condition %in% "control" & 
                                                                                                                dfAnnotation$Sample %in% colnames(matCountData), 
                                                                                                        ]$Sample, collapse = ","), collapse = ","))
                        }
                }
                if (!is.null(vecConfounders)) {
                        for (confounder in vecConfounders) {
                                for (batch in unique(dfAnnotation[, confounder])) {
                                        strReportProcessing <- paste0(strReportProcessing, 
                                                                      "\n", "Found the following samples for confounder ", 
                                                                      confounder, " and batch ", batch, ": ", paste0(dfAnnotation[dfAnnotation[, 
                                                                                                                                               confounder] %in% batch & dfAnnotation$Sample %in% 
                                                                                                                                          colnames(matCountData), ]$Sample, collapse = ","))
                                }
                        }
                }
                return(strReportProcessing)
        }
        procAnnotation <- function(dfAnnotation, matCountData, boolCaseCtrl, 
                                   vecConfounders) {
                for (col in seq(1, dim(dfAnnotation)[2])) dfAnnotation[, 
                                                                       col] <- as.vector(dfAnnotation[, col])
                if (!boolCaseCtrl) {
                        dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition == 
                                                                 "case", ]
                }
                else {
                        dfAnnotationProc <- dfAnnotation[dfAnnotation$Condition == 
                                                                 "case" | dfAnnotation$Condition == "control", 
                                                         ]
                }
                dfAnnotationProc <- dfAnnotationProc[dfAnnotationProc$Sample %in% 
                                                             colnames(matCountData), ]
                dfAnnotationProc <- dfAnnotationProc[, c("Sample", "Time", 
                                                         "Condition", vecConfounders)]
                dfAnnotationProc$TimeCateg <- paste0("_", dfAnnotationProc$Time)
                if (boolCaseCtrl) {
                        for (confounder in vecConfounders) {
                                vecBatchesCase <- dfAnnotationProc[dfAnnotationProc$Condition == 
                                                                           "case", confounder]
                                vecNestedBatchesCase <- paste0("NestedBatch", 
                                                               match(vecBatchesCase, unique(vecBatchesCase)))
                                vecBatchesCtrl <- dfAnnotationProc[dfAnnotationProc$Condition == 
                                                                           "control", confounder]
                                vecNestedBatchesCtrl <- paste0("NestedBatch", 
                                                               match(vecBatchesCtrl, unique(vecBatchesCtrl)))
                                vecBatchesNested <- array(NA, dim(dfAnnotationProc)[1])
                                vecBatchesNested[dfAnnotationProc$Condition == 
                                                         "case"] <- vecNestedBatchesCase
                                vecBatchesNested[dfAnnotationProc$Condition == 
                                                         "control"] <- vecNestedBatchesCtrl
                                strNameConfounderNested <- paste0(confounder, 
                                                                  "Nested")
                                dfAnnotationProc[strNameConfounderNested] <- vecBatchesNested
                        }
                }
                return(dfAnnotationProc)
        }
        nameGenes <- function(matCountDataProc) {
                if (is.null(rownames(matCountDataProc))) {
                        rownames(matCountDataProc) <- paste0("Gene_", seq(1, 
                                                                          nrow(matCountDataProc)))
                }
                return(matCountDataProc)
        }
        reduceCountData <- function(dfAnnotation, matCountDataProc) {
                strReportCountRed <- paste0("Input contained ", dim(matCountDataProc)[1], 
                                            " genes/regions.")
                if (boolCaseCtrl) {
                        vecSampleNames_Case <- as.character(as.vector(dfAnnotation[dfAnnotation$Condition %in% 
                                                                                           "case" & dfAnnotation$Sample %in% colnames(matCountDataProc), 
                                                                                   ]$Sample))
                        vecSampleNames_Ctrl <- as.character(as.vector(dfAnnotation[dfAnnotation$Condition %in% 
                                                                                           "control" & dfAnnotation$Sample %in% colnames(matCountDataProc), 
                                                                                   ]$Sample))
                        matCountDataProc <- matCountDataProc[, c(vecSampleNames_Case, 
                                                                 vecSampleNames_Ctrl)]
                }
                else {
                        vecSampleNames_Case <- as.character(as.vector(dfAnnotation[dfAnnotation$Condition %in% 
                                                                                           "case" & dfAnnotation$Sample %in% colnames(matCountDataProc), 
                                                                                   ]$Sample))
                        matCountDataProc <- matCountDataProc[, vecSampleNames_Case]
                }
                vecboolNASample <- apply(matCountDataProc, 2, function(sample) all(is.na(sample)))
                if (any(vecboolNASample)) {
                        strReportCountRed <- paste0(strReportCountRed, "\n", 
                                                    "WARNING: Sample(s) ", paste0(colnames(matCountDataProc)[vecboolNASample], 
                                                                                  collapse = ","), " only contain(s) NA values and will be ", 
                                                    "removed from the analysis.")
                        matCountDataProc <- matCountDataProc[, !vecboolNASample]
                }
                vecboolAllZeroSample <- apply(matCountDataProc, 2, function(sample) all(sample[!is.na(sample)] == 
                                                                                                0))
                if (any(vecboolAllZeroSample)) {
                        strReportCountRed <- paste0(strReportCountRed, "\n", 
                                                    "WARNING: Sample(s) ", paste0(colnames(matCountDataProc)[vecboolAllZeroSample], 
                                                                                  collapse = ","), " only contain(s) zeros (and NAs).", 
                                                    " These samples are kept for analysis.")
                }
                vecboolNonZeroGene <- apply(matCountDataProc, 1, function(gene) any(!is.na(gene) & 
                                                                                            gene > 0))
                matCountDataProc <- matCountDataProc[vecboolNonZeroGene, 
                                                     ]
                if (sum(!vecboolNonZeroGene) > 0) {
                        strReportCountRed <- paste0(strReportCountRed, "\n", 
                                                    "WARNING: ", sum(!vecboolNonZeroGene), " out of ", 
                                                    length(vecboolNonZeroGene), " genes do not have obserserved non-zero counts", 
                                                    " and are excluded.")
                }
                matCountDataProc <- matCountDataProc[, match(as.vector(dfAnnotation$Sample), 
                                                             colnames(matCountDataProc))]
                strReportCountRed <- paste0(strReportCountRed, "\n", 
                                            "Selected ", dim(matCountDataProc)[1], " genes/regions for analysis.")
                return(list(matCountDataProc = matCountDataProc, strReportCountRed = strReportCountRed))
        }
        strReportProcessing <- checkData(dfAnnotation = dfAnnotation, 
                                         matCountData = matCountData, boolCaseCtrl = boolCaseCtrl, 
                                         vecConfounders = vecConfounders, vecDispersionsExternal = vecDispersionsExternal, 
                                         vecSizeFactorsExternal = vecSizeFactorsExternal)
        dfAnnotationProc <- procAnnotation(dfAnnotation = dfAnnotation, 
                                           matCountData = matCountData, boolCaseCtrl = boolCaseCtrl, 
                                           vecConfounders = vecConfounders)
        matCountDataProc <- nameGenes(matCountDataProc = matCountData)
        lsReduceCounts <- reduceCountData(dfAnnotation = dfAnnotationProc, 
                                          matCountDataProc = matCountDataProc)
        matCountDataProc <- lsReduceCounts$matCountDataProc
        strReportProcessing <- paste0(strReportProcessing, "\n", 
                                      lsReduceCounts$strReportCountRed)
        if (!is.null(vecSizeFactorsExternal)) {
                vecSizeFactorsExternalProc <- vecSizeFactorsExternal[colnames(matCountDataProc)]
        }
        else {
                vecSizeFactorsExternalProc <- NULL
        }
        if (!is.null(vecDispersionsExternal)) {
                vecDispersionsExternalProc <- vecDispersionsExternal[rownames(matCountDataProc)]
        }
        else {
                vecDispersionsExternalProc <- NULL
        }
        return(list(matCountDataProc = matCountDataProc, dfAnnotationProc = dfAnnotationProc, 
                    vecSizeFactorsExternalProc = vecSizeFactorsExternalProc, 
                    vecDispersionsExternalProc = vecDispersionsExternalProc, 
                    strReportProcessing = strReportProcessing))
}
################################################################################

ImpulseDE2:::runDESeq2
################################################################################
runDESeq2 <- function (dfAnnotationProc, matCountDataProc, boolCaseCtrl, vecConfounders) 
{
        dds <- NULL
        if (!boolCaseCtrl) {
                if (is.null(vecConfounders)) {
                        dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                       colData = dfAnnotationProc, design = ~TimeCateg))
                        dds <- estimateSizeFactors(dds)
                        dds <- estimateDispersions(dds)
                }
                else {
                        tryCatch({
                                dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                               colData = dfAnnotationProc, design = ~TimeCateg + 
                                                                                       Batch))
                                dds <- estimateSizeFactors(dds)
                                dds <- estimateDispersions(dds)
                        }, error = function(strErrorMsg) {
                                print(strErrorMsg)
                                print(paste0("WARNING: DESeq2 failed on full model ", 
                                             "- dispersions may be inaccurate.", "Estimating dispersions on reduced model ", 
                                             "formulation [full = ~Batch", " reduced = ~1]. Supply externally generated", 
                                             " dispersion parameters via ", "vecDispersionsExternal if there is a more ", 
                                             "accurate model for your data set."))
                                warning("Warning generated in dispersion factor", 
                                        " estimation, read stdout.")
                        }, finally = {
                                if (is.null(dds)) {
                                        dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                                       colData = dfAnnotationProc, design = ~Batch))
                                        dds <- estimateSizeFactors(dds)
                                        dds <- estimateDispersions(dds)
                                }
                        })
                }
        }
        else {
                if (is.null(vecConfounders)) {
                        tryCatch({
                                dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                               colData = dfAnnotationProc, design = ~Condition + 
                                                                                       Condition:TimeCateg))
                                dds <- estimateSizeFactors(dds)
                                dds <- estimateDispersions(dds)
                        }, error = function(strErrorMsg) {
                                print(strErrorMsg)
                                print(paste0("WARNING: DESeq2 failed on full model ", 
                                             "- dispersions may be inaccurate.", "Estimating dispersions on reduced model", 
                                             " formulation [full = ~Condition", " reduced = ~1]. Supply externally ", 
                                             "generated dispersion parameters via ", "vecDispersionsExternal if there is a more", 
                                             " accurate model for your data set."))
                                warning("Warning generated in dispersion factor", 
                                        " estimation, read stdout.")
                        }, finally = {
                                if (is.null(dds)) {
                                        dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                                       colData = dfAnnotationProc, design = ~Condition))
                                        dds <- estimateSizeFactors(dds)
                                        dds <- estimateDispersions(dds)
                                }
                        })
                }
                else {
                        tryCatch({
                                dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                               colData = dfAnnotationProc, design = ~Condition + 
                                                                                       Condition:TimeCateg + Condition:BatchNested))
                                dds <- estimateSizeFactors(dds)
                                dds <- estimateDispersions(dds)
                        }, error = function(strErrorMsg) {
                                print(strErrorMsg)
                                print(paste0("WARNING: DESeq2 failed on full model ", 
                                             "- dispersions may be inaccurate.", "Estimating dispersions on reduced model.", 
                                             " Supply externally generated dispersion parameters via ", 
                                             "vecDispersionsExternal if there is a more accurate ", 
                                             "model for your data set."))
                                warning("Warning generated in dispersion factor estimation.", 
                                        " Eead stdout.")
                        }, finally = {
                                if (is.null(dds)) {
                                        matModelMatrixBatches <- do.call(cbind, lapply(vecConfounders, 
                                                                                       function(confounder) {
                                                                                               match(dfAnnotationProc[, confounder], unique(dfAnnotationProc[, 
                                                                                                                                                             confounder]))
                                                                                       }))
                                        matModelMatrixBatchesTime <- cbind(matModelMatrixBatches, 
                                                                           match(dfAnnotationProc$Time, unique(dfAnnotationProc$Time)))
                                        boolFullRankBatchTime <- rankMatrix(matModelMatrixBatchesTime)[1] == 
                                                dim(matModelMatrixBatchesTime)[2]
                                        if (!boolFullRankBatchTime) {
                                                paste0("Model matrix based on confounding variables {", 
                                                       vecConfounders, "} and time is not full rank: ", 
                                                       "There are confounding variables with ", 
                                                       " batch structures which are linear", " combinations of the time points.")
                                                paste0("Using reduced model formulation ", 
                                                       "[full= ~Condition+Condition:TimeCateg, ", 
                                                       "reduced= ~TimeCateg].")
                                                dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                                               colData = dfAnnotationProc, design = ~Condition + 
                                                                                                       Condition:TimeCateg))
                                                dds <- estimateSizeFactors(dds)
                                                dds <- estimateDispersions(dds)
                                        }
                                        else {
                                                paste0("Found 1. or {1. and 2.}:")
                                                paste0("1. Model matrix based on confounding ", 
                                                       "variables {", vecConfounders, "} and time is not full rank:", 
                                                       " There are confounding variables with ", 
                                                       " batch structures which are ", "linear combinations of the time points.")
                                                paste0(" 2. Model matrix based on condition", 
                                                       " and time is not full rank: ", "Conditions case and control have", 
                                                       " mutually exclusive sets of timepoints.")
                                                paste0("Using reduced model formulation ", 
                                                       "[full= ~Condition+Condition:BatchNested, ", 
                                                       "reduced= ~1].")
                                                dds <- suppressWarnings(DESeqDataSetFromMatrix(countData = matCountDataProc, 
                                                                                               colData = dfAnnotationProc, design = ~Condition + 
                                                                                                       Condition:BatchNested))
                                                dds <- estimateSizeFactors(dds)
                                                dds <- estimateDispersions(dds)
                                        }
                                }
                        })
                }
        }
        vecDispersionsInv <- mcols(dds)$dispersion
        vecindDESeq2HighOutliesFailure <- !is.na(mcols(dds)$dispOutlier) & 
                mcols(dds)$dispOutlier == TRUE & mcols(dds)$dispersion > 
                (20 - 10^(-5)) & apply(matCountDataProc, 1, function(gene) any(gene == 0))
        vecDispersionsInv[vecindDESeq2HighOutliesFailure] <- mcols(dds)$dispMAP[vecindDESeq2HighOutliesFailure]
        if (sum(vecindDESeq2HighOutliesFailure) > 0) {
                print(paste0("Corrected ", sum(vecindDESeq2HighOutliesFailure), 
                             " DESEq2 dispersion estimates which ", "to avoid variance overestimation and loss of ", 
                             "discriminatory power for model selection."))
        }
        vecDispersions <- 1/vecDispersionsInv
        names(vecDispersions) <- rownames(dds)
        return(vecDispersions)
}
################################################################################


ImpulseDE2:::fitConstImpulse
################################################################################
fitConstImpulse <- function (matCountDataProcCondition, 
                             vecDispersions, 
                             vecSizeFactors, 
                             vecTimepoints, 
                             lsvecBatches, 
                             boolFitConst){
        MAXIT <- 1000
        vecTimepointsUnique <- sort(unique(vecTimepoints))
        vecidxTimepoint <- match(vecTimepoints, vecTimepointsUnique)
        if (length(lsvecBatches) > 0) {
                lsvecBatchUnique <- list()
                lsvecidxBatch <- list()
                for (batchfactor in lsvecBatches) {
                        lsvecBatchUnique[[length(lsvecBatchUnique) + 1]] <- unique(batchfactor)
                        lsvecidxBatch[[length(lsvecidxBatch) + 1]] <- match(batchfactor, 
                                                                            unique(batchfactor))
                }
                names(lsvecBatchUnique) <- names(lsvecBatches)
                names(lsvecidxBatch) <- names(lsvecBatches)
        }
        else {
                lsvecBatchUnique <- NULL
                lsvecidxBatch <- NULL
        }
        lsFits <- bplapply(rownames(matCountDataProcCondition), function(x) {
                fitConstImpulseGene(
                        vecCounts = matCountDataProcCondition[x,], 
                        scaDisp = vecDispersions[x], 
                        vecSizeFactors = vecSizeFactors, 
                        vecTimepointsUnique = vecTimepointsUnique, 
                        vecidxTimepoint = vecidxTimepoint, 
                        lsvecidxBatch = lsvecidxBatch, 
                        boolFitConst = boolFitConst, 
                        MAXIT = MAXIT)
        })
        names(lsFits) <- rownames(matCountDataProcCondition)
        return(list(lsFits = lsFits, vecTimepointsUnique = vecTimepointsUnique, 
                    vecidxTimepoint = vecidxTimepoint, lsvecBatchUnique = lsvecBatchUnique, 
                    lsvecidxBatch = lsvecidxBatch))
}
################################################################################

ImpulseDE2:::evalImpulse_comp
################################################################################
evalImpulse_comp <- function (vecImpulseParam, vecTimepoints) 
{
        vecImpulseValue <- sapply(vecTimepoints, function(t) {
                (1/vecImpulseParam[3]) * (vecImpulseParam[2] + (vecImpulseParam[3] - 
                                                                        vecImpulseParam[2]) * (1/(1 + exp(-vecImpulseParam[1] * 
                                                                                                                  (t - vecImpulseParam[5]))))) * (vecImpulseParam[4] + 
                                                                                                                                                          (vecImpulseParam[3] - vecImpulseParam[4]) * (1/(1 + 
                                                                                                                                                                                                                  exp(vecImpulseParam[1] * (t - vecImpulseParam[6])))))
        })
        vecImpulseValue[vecImpulseValue < 10^(-10)] <- 10^(-10)
        return(vecImpulseValue)
}
################################################################################

ImpulseDE2:::evalLogLikImpulse_comp
################################################################################
evalLogLikImpulse_comp <- function (vecTheta, vecCounts, scaDisp, vecSizeFactors, vecTimepointsUnique, 
                                    vecidxTimepoint, lsvecidxBatch, vecboolObserved) 
{
        vecImpulseParam <- vecTheta[1:6]
        vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
        vecImpulseValue <- evalImpulse_comp(vecImpulseParam = vecImpulseParam, 
                                            vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
        scaNParamUsed <- 6
        vecBatchFactors <- array(1, length(vecCounts))
        if (!is.null(lsvecidxBatch)) {
                for (vecidxConfounder in lsvecidxBatch) {
                        scaNBatchFactors <- max(vecidxConfounder) - 1
                        vecBatchFacConf <- c(1, exp(vecTheta[(scaNParamUsed + 
                                                                      1):(scaNParamUsed + scaNBatchFactors)]))[vecidxConfounder]
                        scaNParamUsed <- scaNParamUsed + scaNBatchFactors
                        vecBatchFacConf[vecBatchFacConf < 10^(-10)] <- 10^(-10)
                        vecBatchFacConf[vecBatchFacConf > 10^(10)] <- 10^(10)
                        vecBatchFactors <- vecBatchFactors * vecBatchFacConf
                }
        }
        scaLogLik <- sum(dnbinom(vecCounts[vecboolObserved], mu = vecImpulseValue[vecboolObserved] * 
                                         vecBatchFactors[vecboolObserved] * vecSizeFactors[vecboolObserved], 
                                 size = scaDisp, log = TRUE))
        return(scaLogLik)
}
################################################################################

ImpulseDE2:::fitConstModel
################################################################################
fitConstModel <- function (vecCounts, scaDisp, vecSizeFactors, lsvecidxBatch, 
          MAXIT = 1000, RELTOL = 10^(-8), trace = 0, REPORT = 10) 
{
        vecParamGuess <- log(mean(vecCounts, na.rm = TRUE) + 1)
        if (!is.null(lsvecidxBatch)) {
                for (vecidxConfounder in lsvecidxBatch) {
                        vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder)) - 
                                                                      1))
                }
        }
        lsFit <- tryCatch({
                optim(par = vecParamGuess, fn = evalLogLikMu_comp, vecCounts = vecCounts, 
                      scaDisp = scaDisp, vecSizeFactors = vecSizeFactors, 
                      lsvecidxBatch = lsvecidxBatch, vecboolObserved = !is.na(vecCounts), 
                      method = "BFGS", control = list(maxit = MAXIT, reltol = RELTOL, 
                                                      fnscale = -1))[c("par", "value", "convergence")]
        }, error = function(strErrorMsg) {
                print(paste0("ERROR: Fitting null model: fitConstModel().", 
                             " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
                print(paste0("vecParamGuess ", paste(vecParamGuess, collapse = " ")))
                print(paste0("vecCounts ", paste(vecCounts, collapse = " ")))
                print(paste0("scaDisp ", paste(scaDisp, collapse = " ")))
                print(paste0("vecSizeFactors ", paste(vecSizeFactors, 
                                                      collapse = " ")))
                print(paste0("lsvecidxBatch ", paste(lsvecidxBatch, collapse = " ")))
                print(paste0("MAXIT ", MAXIT))
                print(strErrorMsg)
                stop(strErrorMsg)
        })
        scaMu <- exp(lsFit$par[1])
        if (scaMu < 10^(-10)) {
                scaMu <- 10^(-10)
        }
        scaNParamUsed <- 1
        if (!is.null(lsvecidxBatch)) {
                lsvecBatchFactors <- list()
                for (i in seq(1, length(lsvecidxBatch))) {
                        vecidxConfounder <- lsvecidxBatch[[i]]
                        scaNBatchFactors <- max(vecidxConfounder) - 1
                        vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed + 
                                                                       1):(scaNParamUsed + scaNBatchFactors)]))
                        scaNParamUsed <- scaNParamUsed + scaNBatchFactors
                        vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                        vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                        lsvecBatchFactors[[i]] <- vecBatchFactors
                }
        }
        else {
                lsvecBatchFactors <- NULL
        }
        return(list(scaMu = scaMu, lsvecBatchFactors = lsvecBatchFactors, 
                    scaDispParam = scaDisp, scaLL = lsFit$value, scaConvergence = lsFit$convergence))
}
################################################################################


ImpulseDE2:::fitConstImpulseGene
################################################################################
fitConstImpulseGene <- function (vecCounts, 
                                 scaDisp, 
                                 vecSizeFactors, 
                                 vecTimepointsUnique, 
                                 vecidxTimepoint,
                                 lsvecidxBatch, 
                                 boolFitConst, 
                                 MAXIT = 1000) {
        lsParamGuesses <- estimateImpulseParam(vecCounts = vecCounts, 
                                               vecTimepoints = vecTimepointsUnique[vecidxTimepoint], 
                                               lsvecidxBatch = lsvecidxBatch, vecSizeFactors = vecSizeFactors)
        vecParamGuessPeak <- lsParamGuesses$peak
        vecParamGuessValley <- lsParamGuesses$valley
        lsFitPeak <- fitImpulseModel(vecImpulseParamGuess = vecParamGuessPeak, 
                                     vecCounts = vecCounts, 
                                     scaDisp = scaDisp, 
                                     vecSizeFactors = vecSizeFactors, 
                                     vecTimepointsUnique = vecTimepointsUnique, 
                                     vecidxTimepoint = vecidxTimepoint, 
                                     lsvecidxBatch = lsvecidxBatch, 
                                     MAXIT = MAXIT)
        lsFitValley <- fitImpulseModel(vecImpulseParamGuess = vecParamGuessValley, 
                                       vecCounts = vecCounts, 
                                       scaDisp = scaDisp, 
                                       vecSizeFactors = vecSizeFactors, 
                                       vecTimepointsUnique = vecTimepointsUnique, 
                                       vecidxTimepoint = vecidxTimepoint, 
                                       lsvecidxBatch = lsvecidxBatch, 
                                       MAXIT = MAXIT)
        if (lsFitValley$scaLL > lsFitPeak$scaLL) {
                lsBestImpulseFit <- lsFitValley
        }
        else {
                lsBestImpulseFit <- lsFitPeak
        }
        if (boolFitConst) {
                lsConstFit <- fitConstModel(vecCounts = vecCounts, scaDisp = scaDisp, 
                                            vecSizeFactors = vecSizeFactors, lsvecidxBatch = lsvecidxBatch)
        }
        else {
                lsConstFit <- NULL
        }
        return(list(lsImpulseFit = lsBestImpulseFit, lsConstFit = lsConstFit))
}
################################################################################

ImpulseDE2:::estimateImpulseParam
################################################################################
estimateImpulseParam <- function (vecCounts, vecTimepoints, vecSizeFactors, lsvecidxBatch) 
{
        vecTimepointsUnique <- unique(vecTimepoints)
        scaMeanCount <- mean(vecCounts, na.rm = TRUE)
        vecCountsSFcorrected <- vecCounts/vecSizeFactors
        vecCountsSFcorrectedNorm <- vecCountsSFcorrected/scaMeanCount
        if (!is.null(lsvecidxBatch)) {
                vecBatchFactors <- array(1, length(vecCounts))
                for (vecidxBatch in lsvecidxBatch) {
                        vecBatchFactorsConfounder <- tapply(vecCountsSFcorrectedNorm, 
                                                            vecidxBatch, mean, na.rm = TRUE)
                        vecBatchFactorsConfounder[is.na(vecBatchFactorsConfounder) | 
                                                          vecBatchFactorsConfounder == 0] <- 1
                        vecBatchFactors <- vecBatchFactors * vecBatchFactorsConfounder[vecidxBatch]
                }
                vecCountsSFBatchcorrected <- vecCountsSFcorrected/vecBatchFactors
                vecExpressionMeans <- tapply(vecCountsSFBatchcorrected, 
                                             match(vecTimepoints, vecTimepointsUnique), mean, 
                                             na.rm = TRUE)
        }
        else {
                vecExpressionMeans <- tapply(vecCountsSFcorrected, match(vecTimepoints, 
                                                                         vecTimepointsUnique), mean, na.rm = TRUE)
        }
        scaNTimepoints <- length(vecTimepointsUnique)
        vecidxMiddle <- seq(2, scaNTimepoints - 1, by = 1)
        scaMaxMiddleMean <- max(vecExpressionMeans[vecidxMiddle], 
                                na.rm = TRUE)
        scaMinMiddleMean <- min(vecExpressionMeans[vecidxMiddle], 
                                na.rm = TRUE)
        indMaxMiddleMean <- match(scaMaxMiddleMean, vecExpressionMeans[vecidxMiddle]) + 
                1
        indMinMiddleMean <- match(scaMinMiddleMean, vecExpressionMeans[vecidxMiddle]) + 
                1
        vecDiffExpr <- diff(vecExpressionMeans)
        vecDiffTime <- diff(vecTimepointsUnique)
        vecGradients <- vecDiffExpr/vecDiffTime
        vecGradients[is.na(vecGradients) | !is.finite(vecGradients)] <- 0
        vecdidxFirstPart <- seq(1, indMaxMiddleMean - 1, by = 1)
        vecdidxSecndPart <- seq(indMaxMiddleMean, length(vecGradients), 
                                by = 1)
        indLowerInflexionPoint <- match(max(vecGradients[vecdidxFirstPart], 
                                            na.rm = TRUE), vecGradients[vecdidxFirstPart])
        indUpperInflexionPoint <- indMaxMiddleMean - 1 + match(min(vecGradients[vecdidxSecndPart], 
                                                                   na.rm = TRUE), vecGradients[vecdidxSecndPart])
        vecParamGuessPeak <- c(1, log(vecExpressionMeans[1] + 1), 
                               log(scaMaxMiddleMean + 1), log(vecExpressionMeans[scaNTimepoints] + 
                                                                      1), (vecTimepointsUnique[indLowerInflexionPoint] + 
                                                                                   vecTimepointsUnique[indLowerInflexionPoint + 1])/2, 
                               (vecTimepointsUnique[indUpperInflexionPoint] + vecTimepointsUnique[indUpperInflexionPoint + 
                                                                                                          1])/2)
        indLowerInflexionPoint <- match(min(vecGradients[1:(indMinMiddleMean - 
                                                                    1)], na.rm = TRUE), vecGradients[1:(indMinMiddleMean - 
                                                                                                                1)])
        indUpperInflexionPoint <- indMinMiddleMean - 1 + match(max(vecGradients[indMinMiddleMean:(scaNTimepoints - 
                                                                                                          1)], na.rm = TRUE), vecGradients[indMinMiddleMean:(scaNTimepoints - 
                                                                                                                                                                     1)])
        vecParamGuessValley <- c(1, log(vecExpressionMeans[1] + 1), 
                                 log(scaMinMiddleMean + 1), log(vecExpressionMeans[scaNTimepoints] + 
                                                                        1), (vecTimepointsUnique[indLowerInflexionPoint] + 
                                                                                     vecTimepointsUnique[indLowerInflexionPoint + 1])/2, 
                                 (vecTimepointsUnique[indUpperInflexionPoint] + vecTimepointsUnique[indUpperInflexionPoint + 
                                                                                                            1])/2)
        return(list(peak = vecParamGuessPeak, valley = vecParamGuessValley))
}
################################################################################

ImpulseDE2:::fitImpulseModel
################################################################################
fitImpulseModel <- function (vecImpulseParamGuess, 
                             vecCounts, 
                             scaDisp, 
                             vecSizeFactors, 
                             lsvecidxBatch, 
                             vecTimepointsUnique, 
                             vecidxTimepoint, 
                             MAXIT = 1000, 
                             RELTOL = 10^(-8), 
                             trace = 0, 
                             REPORT = 10) {
        vecParamGuess <- vecImpulseParamGuess
        if (!is.null(lsvecidxBatch)) {
                for (vecidxConfounder in lsvecidxBatch) {
                        vecParamGuess <- c(vecParamGuess, rep(0, length(unique(vecidxConfounder)) - 
                                                                      1))
                }
        }
        lsFit <- tryCatch({
                optim(par = vecParamGuess, 
                      fn = evalLogLikImpulse_comp, 
                      vecCounts = vecCounts, 
                      scaDisp = scaDisp, 
                      vecSizeFactors = vecSizeFactors, 
                      vecTimepointsUnique = vecTimepointsUnique, 
                      vecidxTimepoint = vecidxTimepoint, 
                      lsvecidxBatch = lsvecidxBatch, 
                      vecboolObserved = !is.na(vecCounts), 
                      method = "BFGS", 
                      control = list(maxit = MAXIT, reltol = RELTOL, fnscale = -1))[c("par", "value", "convergence")]
        }, error = function(strErrorMsg) {
                print(paste0("ERROR: Fitting impulse model: fitImpulseModel().", 
                             " Wrote report into ImpulseDE2_lsErrorCausingGene.RData"))
                print(paste0("vecParamGuess ", paste(vecParamGuess, collapse = " ")))
                print(paste0("vecCounts ", paste(vecCounts, collapse = " ")))
                print(paste0("scaDisp ", paste(scaDisp, collapse = " ")))
                print(paste0("vecSizeFactors ", paste(vecSizeFactors, 
                                                      collapse = " ")))
                print(paste0("vecTimepointsUnique ", paste(vecTimepointsUnique, 
                                                           collapse = " ")))
                print(paste0("vecidxTimepoint ", paste(vecidxTimepoint, 
                                                       collapse = " ")))
                print(paste0("lsvecidxBatch ", paste(lsvecidxBatch, collapse = " ")))
                print(paste0("MAXIT ", MAXIT))
                print(strErrorMsg)
                stop(strErrorMsg)
        })
        vecImpulseParam <- lsFit$par[1:6]
        vecImpulseParam[2:4] <- exp(vecImpulseParam[2:4])
        names(vecImpulseParam) <- c("beta", "h0", "h1", "h2", "t1", "t2")
        vecImpulseValue <- evalImpulse_comp(vecImpulseParam = vecImpulseParam, 
                                            vecTimepoints = vecTimepointsUnique)[vecidxTimepoint]
        names(vecImpulseValue) <- names(vecCounts)
        scaNParamUsed <- 6
        if (!is.null(lsvecidxBatch)) {
                lsvecBatchFactors <- list()
                for (i in seq(1, length(lsvecidxBatch))) {
                        vecidxConfounder <- lsvecidxBatch[[i]]
                        scaNBatchFactors <- max(vecidxConfounder) - 1
                        vecBatchFactors <- c(1, exp(lsFit$par[(scaNParamUsed + 
                                                                       1):(scaNParamUsed + scaNBatchFactors)]))
                        scaNParamUsed <- scaNParamUsed + scaNBatchFactors
                        vecBatchFactors[vecBatchFactors < 10^(-10)] <- 10^(-10)
                        vecBatchFactors[vecBatchFactors > 10^(10)] <- 10^(10)
                        lsvecBatchFactors[[i]] <- vecBatchFactors
                }
        }
        else {
                lsvecBatchFactors <- NULL
        }
        return(list(vecImpulseParam = vecImpulseParam, 
                    vecImpulseValue = vecImpulseValue, 
                    lsvecBatchFactors = lsvecBatchFactors, 
                    scaDispParam = scaDisp, 
                    scaLL = lsFit$value, 
                    scaConvergence = lsFit$convergence))
}
################################################################################


ImpulseDE2:::fitModels
################################################################################
fitModels <- function (objectImpulseDE2, vecConfounders, boolCaseCtrl) 
{
        lsFitResults_all = list()
        dfAnnot <- get_dfAnnotationProc(obj = objectImpulseDE2)
        if (boolCaseCtrl) {
                vecLabels <- c("combined", "case", "control")
        }
        else {
                vecLabels <- c("case")
        }
        if (boolCaseCtrl) {
                lsSamplesByCond <- list(combined = dfAnnot$Sample, 
                                        case = dfAnnot[dfAnnot$Condition == "case", ]$Sample, 
                                        control = dfAnnot[dfAnnot$Condition == "control", ]$Sample)
        }
        else {
                lsSamplesByCond <- list(case = dfAnnot[dfAnnot$Condition == "case", ]$Sample)
        }
        if (!is.null(vecConfounders)) {
                lsvecBatches <- lapply(vecConfounders, function(confounder) {
                        vecBatches <- dfAnnot[, confounder]
                        names(vecBatches) <- dfAnnot$Sample
                        return(vecBatches)
                })
                names(lsvecBatches) <- vecConfounders
        }
        else {
                lsvecBatches <- NULL
        }
        vecTimepoints <- dfAnnot$Time
        names(vecTimepoints) <- colnames(get_matCountDataProc(obj = objectImpulseDE2))
        lsFitResultsByCond <- lapply(vecLabels, function(label) {
                lsFitResults <- fitConstImpulse(
                        matCountDataProcCondition = get_matCountDataProc(obj = objectImpulseDE2)[, lsSamplesByCond[[label]]], 
                        vecSizeFactors = get_vecSizeFactors(obj = objectImpulseDE2)[lsSamplesByCond[[label]]], 
                        vecDispersions = get_vecDispersions(obj = objectImpulseDE2), 
                        vecTimepoints = vecTimepoints[lsSamplesByCond[[label]]], 
                        lsvecBatches = lapply(lsvecBatches, function(confounder) confounder[lsSamplesByCond[[label]]]), 
                        boolFitConst = TRUE)
                return(lsFitResults)
        })
        lsModelFitsByCondFormat <- lapply(lsFitResultsByCond, function(res) res$lsFits)
        names(lsModelFitsByCondFormat) <- vecLabels
        lsModelFitsByCondFormat$IdxGroups <- lapply(lsFitResultsByCond, 
                                                    function(res) {
                                                            list(vecTimepointsUnique = res$vecTimepointsUnique, 
                                                                 vecidxTimepoint = res$vecidxTimepoint, lsvecBatchUnique = res$lsvecBatchUnique, 
                                                                 lsvecidxBatch = res$lsvecidxBatch)
                                                    })
        names(lsModelFitsByCondFormat$IdxGroups) <- vecLabels
        for (label in vecLabels) {
                lsModelFitsByCondFormat$IdxGroups[[label]]$Samples <- lsSamplesByCond[[label]]
        }
        objectImpulseDE2 <- set_lsModelFits(obj = objectImpulseDE2, 
                                            element = lsModelFitsByCondFormat)
        return(objectImpulseDE2)
}
################################################################################


ImpulseDE2:::runImpulseDE2
################################################################################
function (matCountData = lsSimulatedData$matObservedCounts, 
          dfAnnotation = lsSimulatedData$dfAnnotation, 
          boolCaseCtrl = FALSE, 
          vecConfounders = NULL, 
          scaNProc = 1, 
          scaQThres = NULL, 
          vecDispersionsExternal = NULL, 
          vecSizeFactorsExternal = NULL, 
          boolIdentifyTransients = FALSE, 
          boolVerbose = TRUE) 
{
        strMessage <- paste0("ImpulseDE2 for count data, v", packageDescription("ImpulseDE2", 
                                                                                fields = "Version"))
        if (boolVerbose) {
                message(strMessage)
        }
        strReport <- strMessage
        tm_runImpulseDE2 <- system.time({
                strMessage <- "# Process input"
                if (boolVerbose) {
                        message(strMessage)
                }
                strReport <- paste0(strReport, "\n", strMessage)
                if (class(matCountData) == "SummarizedExperiment") {
                        matCountData <- assay(matCountData)
                }
                lsProcessedData <- processData(dfAnnotation = dfAnnotation, 
                                               matCountData = matCountData, 
                                               boolCaseCtrl = boolCaseCtrl, 
                                               vecConfounders = vecConfounders, 
                                               vecDispersionsExternal = vecDispersionsExternal, 
                                               vecSizeFactorsExternal = vecSizeFactorsExternal)
                matCountDataProc <- lsProcessedData$matCountDataProc
                dfAnnotationProc <- lsProcessedData$dfAnnotationProc
                vecSizeFactorsExternalProc <- lsProcessedData$vecSizeFactorsExternalProc
                vecDispersionsExternalProc <- lsProcessedData$vecDispersionsExternalProc
                if (boolVerbose) {
                        write(lsProcessedData$strReportProcessing, file = "", 
                              ncolumns = 1)
                }
                strReport <- paste0(strReport, lsProcessedData$strReportProcessing)
                if (scaNProc > 1) {
                        register(MulticoreParam(workers = scaNProc))
                }
                else {
                        register(SerialParam())
                }
                if (is.null(vecDispersionsExternal)) {
                        strMessage <- paste0("# Run DESeq2: Using dispersion factors", 
                                             "computed by DESeq2.")
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        strReport <- paste0(strReport, "\n", strMessage)
                        tm_runDESeq2 <- system.time({
                                vecDispersions <- runDESeq2(dfAnnotationProc = dfAnnotationProc, 
                                                            matCountDataProc = matCountDataProc, 
                                                            boolCaseCtrl = boolCaseCtrl, 
                                                            vecConfounders = vecConfounders)
                        })
                        strMessage <- paste0("Consumed time: ", round(tm_runDESeq2["elapsed"]/60, 
                                                                      2), " min.")
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        strReport <- paste0(strReport, "\n", strMessage)
                }
                else {
                        strMessage <- "# Using externally supplied dispersion factors."
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        strReport <- paste0(strReport, "\n", strMessage)
                        vecDispersions <- vecDispersionsExternalProc
                }
                strMessage <- "# Compute size factors"
                if (boolVerbose) {
                        message(strMessage)
                }
                strReport <- paste0(strReport, "\n", strMessage)
                vecSizeFactors <- computeNormConst(matCountDataProc = matCountDataProc, 
                                                   vecSizeFactorsExternal = vecSizeFactorsExternalProc)
                objectImpulseDE2 <- new("ImpulseDE2Object", 
                                        dfImpulseDE2Results = NULL, 
                                        vecDEGenes = NULL, 
                                        lsModelFits = NULL, 
                                        matCountDataProc = matCountDataProc, 
                                        vecAllIDs = rownames(matCountData), 
                                        dfAnnotationProc = dfAnnotationProc, 
                                        vecSizeFactors = vecSizeFactors, 
                                        vecDispersions = vecDispersions, 
                                        boolCaseCtrl = boolCaseCtrl, 
                                        vecConfounders = vecConfounders, 
                                        scaNProc = scaNProc, 
                                        scaQThres = scaQThres, 
                                        strReport = strReport)
                strMessage <- "# Fitting null and alternative model to the genes"
                if (boolVerbose) {
                        message(strMessage)
                }
                objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                     s = strMessage)
                tm_fitImpulse <- system.time({
                        objectImpulseDE2 <- fitModels(objectImpulseDE2 = objectImpulseDE2, 
                                                      vecConfounders = vecConfounders, 
                                                      boolCaseCtrl = boolCaseCtrl)
                })
                strMessage <- paste0("Consumed time: ", round(tm_fitImpulse["elapsed"]/60, 
                                                              2), " min.")
                if (boolVerbose) {
                        message(strMessage)
                }
                objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                     s = strMessage)
                if (boolIdentifyTransients) {
                        strMessage <- "# Fitting sigmoid model to case condition"
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                             s = strMessage)
                        tm_fitSigmoid <- system.time({
                                objectImpulseDE2 <- fitSigmoidModels(objectImpulseDE2 = objectImpulseDE2, 
                                                                     vecConfounders = vecConfounders, strCondition = "case")
                        })
                        strMessage <- paste0("Consumed time: ", round(tm_fitSigmoid["elapsed"]/60, 
                                                                      2), " min.")
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                             s = strMessage)
                }
                strMessage <- "# Differentially expression analysis based on model fits"
                if (boolVerbose) {
                        message(strMessage)
                }
                objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                     s = strMessage)
                objectImpulseDE2 <- runDEAnalysis(objectImpulseDE2 = objectImpulseDE2, 
                                                  boolCaseCtrl = get_boolCaseCtrl(obj = objectImpulseDE2), 
                                                  boolIdentifyTransients = boolIdentifyTransients)
                if (!is.null(scaQThres)) {
                        vecDEGenes <- as.vector(objectImpulseDE2$dfImpulseDE2Results[as.numeric(objectImpulseDE2$dfImpulseDE2Results$padj) <= scaQThres, "Gene"])
                        strMessage <- paste0("Found ", length(vecDEGenes), 
                                             " DE genes", " at a FDR corrected p-value cut off of ", 
                                             scaQThres, ".")
                        if (boolVerbose) {
                                message(strMessage)
                        }
                        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                                             s = strMessage)
                }
                else {
                        vecDEGenes <- NULL
                }
                objectImpulseDE2 <- set_vecDEGenes(obj = objectImpulseDE2, 
                                                   element = vecDEGenes)
        })
        strMessage <- "Finished running ImpulseDE2."
        if (boolVerbose) {
                message(strMessage)
        }
        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                             s = strMessage)
        strMessage <- paste0("TOTAL consumed time: ", round(tm_runImpulseDE2["elapsed"]/60, 
                                                            2), " min.")
        if (boolVerbose) {
                message(strMessage)
        }
        objectImpulseDE2 <- append_strReport(obj = objectImpulseDE2, 
                                             s = strMessage)
        return(objectImpulseDE2)
}
################################################################################

        

