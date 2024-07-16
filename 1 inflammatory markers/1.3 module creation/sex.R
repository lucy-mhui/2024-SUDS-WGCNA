# the following code is for MDD male - adjust as needed for MDD female and Control male and female

##############
# data input #
##############

# library
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# load data frame 
lnames = load(file = "/.../sex.RData")
lnames 




#####################
# soft thresholding #
#####################

# choose a set of soft thresholding powers
powers=c(c(1:10), seq(from=12, to=40, by=2))

# call network topology analysis function
sft = pickSoftThreshold(MDDinfM, powerVector = powers, verbose = 5)

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

net = blockwiseModules(MDDinfM, power = 6,
                       TOMType = "signed", minModuleSize = 1,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)

# module sizes and number of modules
table(net$colors)
