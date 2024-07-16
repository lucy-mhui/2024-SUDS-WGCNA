# the following code is for MDD - adjust as neede for Control

#############
# load data #
#############

# library
library(WGCNA)
options(stringsAsFactors=FALSE)

# load variables
lnames=load(file="/.../MDDinfMadrs.RData")
lnames
lnames=load(file="/.../02-MDD-networkConstruction.RData")
lnames




######################################
# quantify modeul trait associations #
######################################

# define numbers of genes and samples
nMarkers=ncol(MDDinfData5)
nSamples=nrow(MDDinfData5)

# recalculate MEs with color labels
MEs0=moduleEigengenes(MDDinfData5,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,MDDmadrs1,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

# display correlations and p values
textMatrix = paste(signif(moduleTraitCor,2), "\n(",
                   signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDDmadrs1),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))




###################################
# visualize network of eigengenes #
###################################

# Recalculate module eigengenes
MEs = moduleEigengenes(MDDinfData5, moduleColors)$eigengenes
# Isolate weight from the clinical traits
ef = as.data.frame(MDDmadrs1$`reported sadness`);
names(ef) = "reported sadness"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, ef))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
