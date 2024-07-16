# load data frames
library(WGCNA)

male_lnames=load(file="/.../dtiM.RData")
female_lnames=load(file="/.../dtiF.RData")
MDDinf_lnames=load(file="/.../MDDinfMadrs.RData")
MDDwgcna_lnames=load(file="/.../02-MDD-networkConstruction.RData")

new_row_names <- gsub("^CBN01_", "", rownames(MDDinfData5))
rownames(MDDinfData5) <- new_row_names

# MDD_cdi_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_cdi_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_cdi_M <- MDD_cdi_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_cdi_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_cdi_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_dti_FA_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_dti_FA_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_dti_FA_M <- MDD_dti_FA_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_dti_FA_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_dti_FA_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_dti_MD_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_dti_MD_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_dti_MD_M <- MDD_dti_MD_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_dti_MD_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_dti_MD_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_F_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_F_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_F_M <- MDD_F_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_F_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_F_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_FAt_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_FAt_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_FAt_M <- MDD_FAt_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_FAt_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_FAt_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_MD_t_M
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_MD_t_M))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_MD_t_M <- MDD_MD_t_M[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
new_order <- c("MEbrown", "MEgreen", "MEyellow", "MEblue", "MEturquoise", "MEgrey")
MEs <- MEs[, new_order]
moduleTraitCor <- cor(MEs, MDD_MD_t_M, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_cdi_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_cdi_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_cdi_F <- MDD_cdi_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_cdi_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_dti_FA_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_dti_FA_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_dti_FA_F <- MDD_dti_FA_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_dti_FA_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_dti_MD_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_dti_MD_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_dti_MD_F <- MDD_dti_MD_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_dti_MD_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_F_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_F_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_F_F <- MDD_F_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_F_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_FAt_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_FAt_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_FAt_F <- MDD_FAt_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_FAt_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)

# MDD_MD_t_F
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_MD_t_F))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_MD_t_F <- MDD_MD_t_F[common_subject_ids, ]

nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_MD_t_F, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10
significant_pvalue_mask <- moduleTraitPvalue < 0.05
combined_mask <- correlation_mask & significant_pvalue_mask
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
textMatrix[!combined_mask] <- ""
moduleTraitCor[!combined_mask] <- 0
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_M),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               xLabelsAngle = 80,
               plotLegend = FALSE)
