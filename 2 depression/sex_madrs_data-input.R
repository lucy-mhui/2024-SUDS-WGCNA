# the following code is for MDD - adjust as needed for Control 

# load
library(WGCNA)
lnames=load(file="/.../MDDinfMadrs.RData")
lnames=load(file="/.../CntrlinfMadrs.RData")
lnames = load(file = "/.../sex.RData")

# MDDinfM
lnames=load(file="/Users/lucyhui/Downloads/SUDS/Individual/2 module detection/02-MDDinfM-networkConstruction.RData")

common_row_names <- intersect(rownames(MDDinfM), rownames(MDDmadrs1))
MDDmadrsM <- MDDmadrs1[common_row_names, , drop = FALSE]

nMarkers=ncol(MDDinfM)
nSamples=nrow(MDDinfM)

MEs0=moduleEigengenes(MDDinfM,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,MDDmadrsM,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDDmadrsM),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

library(dplyr)
ME_MDD_M = MEs %>% mutate(Gender = "M")

# MDDinfF
lnames=load(file="/Users/lucyhui/Downloads/SUDS/Individual/2 module detection/02-MDDinfF-networkConstruction.RData")

common_row_names <- intersect(rownames(MDDinfF), rownames(MDDmadrs1))
MDDmadrsF <- MDDmadrs1[common_row_names, , drop = FALSE]

nMarkers=ncol(MDDinfF)
nSamples=nrow(MDDinfF)

MEs0=moduleEigengenes(MDDinfF,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,MDDmadrsF,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDDmadrsF),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

library(dplyr)
ME_MDD_F = MEs %>% mutate(Gender = "F")

# CntrlinfM
lnames=load(file="/Users/lucyhui/Downloads/SUDS/Individual/2 module detection/02-CntrlInfM-networkConstruction.RData")

common_row_names <- intersect(rownames(CntrlinfM), rownames(Cntrlmadrs1))
CntrlmadrsM <- Cntrlmadrs1[common_row_names, , drop = FALSE]
CntrlinfM <- CntrlinfM[common_row_names, , drop = FALSE]

nMarkers=ncol(CntrlinfM)
nSamples=nrow(CntrlinfM)

MEs0=moduleEigengenes(CntrlinfM,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,CntrlmadrsM,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(CntrlmadrsM),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# CntrlinfF
lnames=load(file="/Users/lucyhui/Downloads/SUDS/Individual/2 module detection/02-CntrlInfF-networkConstruction.RData")

common_row_names <- intersect(rownames(CntrlinfF), rownames(Cntrlmadrs1))
CntrlmadrsF <- Cntrlmadrs1[common_row_names, , drop = FALSE]
CntrlinfF <- CntrlinfF[common_row_names, , drop = FALSE]

nMarkers=ncol(CntrlinfF)
nSamples=nrow(CntrlinfF)

MEs0=moduleEigengenes(CntrlinfF,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,CntrlmadrsF,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(CntrlmadrsF),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Long format
library(dplyr)
library(tidyr)

MDDinfM = MDDinfM %>% mutate(Gender = "M")
MDDinfF = MDDinfF %>% mutate(Gender = "F")
MDDmadrsM = MDDmadrsM %>% mutate(Gender = "M")
MDDmadrsF = MDDmadrsF %>% mutate(Gender = "F")

MDDinf <- rbind(MDDinfM, MDDinfF)
MDDmadrs <- rbind(MDDmadrsM, MDDmadrsF)

MDDinf = data.frame(SubjectID = row.names(MDDinf), MDDinf, row.names = NULL)
MDDmadrs = data.frame(SubjectID = row.names(MDDmadrs), MDDmadrs, row.names = NULL)
ME_MDD_M = data.frame(SubjectID = row.names(ME_MDD_M), ME_MDD_M, row.names = NULL)
ME_MDD_F = data.frame(SubjectID = row.names(ME_MDD_F), ME_MDD_F, row.names = NULL)

MDDinf_long <- pivot_longer(MDDinf,
                            cols = -c(SubjectID, Gender),  
                            names_to = "Marker",          
                            values_to = "Marker Value")
MDDmadrs_long <- pivot_longer(MDDmadrs,
                            cols = -c(SubjectID, Gender),  
                            names_to = "Depression",          
                            values_to = "Depression Value")
ME_MDD_M_long <- pivot_longer(ME_MDD_M,
                              cols = -c(SubjectID, Gender),  
                              names_to = "Modules",          
                              values_to = "Modules Value")
ME_MDD_F_long <- pivot_longer(ME_MDD_F,
                              cols = -c(SubjectID, Gender),  
                              names_to = "Modules",          
                              values_to = "Modules Value")

ME_MDD_long <- rbind(ME_MDD_M_long, ME_MDD_F_long)

merged_data <- merge(MDDinf_long, MDDmadrs_long, by = "SubjectID", all = TRUE)
merged_data = merged_data[,-5]
merged_data <- merge(merged_data, ME_MDD_long, by = "SubjectID", all = TRUE)
merged_data = merged_data[,-7]

colnames(merged_data)=c("SubjectID","sex","marker","marker value","depression component","depression component value")

# save 
save(merged_data, file = "Longformat_MDD_Sex_Module_MADRS.RData")
