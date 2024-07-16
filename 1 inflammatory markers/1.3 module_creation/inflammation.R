# the following code is for MDD - adjust as needed for Control

##################
# initial set up #
##################

# load library
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# load data frame 
lnames = load(file = "/.../MDDinfMadrs.RData")




#####################
# soft thresholding #
#####################

# choose a set of soft thresholding powers
powers=c(c(1:10), seq(from=12, to=40, by=2))

# call network topology analysis function
sft = pickSoftThreshold(MDDinfData5, powerVector = powers, verbose = 5)

# plot the result
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=0.9

# scale free topology fit index as a function of the soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")




#################################
# one step network construction #
#################################

net = blockwiseModules(MDDinfData5, power = 10,
                       TOMType = "signed", minModuleSize = 1,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)

# module sizes and number of modules
table(net$colors)




####################################################################
# plot hierarchial clustering dendrogram for module identification #
####################################################################

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
