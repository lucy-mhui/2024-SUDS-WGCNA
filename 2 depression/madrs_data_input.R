# the following code is for MDD - adjust as needed for Control

# load library
library(dplyr)
library(WGCNA)

# import dataframes
MDDmadrs = read.csv("/.../CBN01_MADRS_DATA_Z3_01_V01_TRTMT.csv")

# select columns
MDDmadrs = MDDmadrs[,-c(2:4,16)]

# restructure
traitRows = match(rownames(MDDinfData5), MDDmadrs$SUBJLABEL)
MDDmadrs1 = MDDmadrs[traitRows, -1]
rownames(MDDmadrs1) = MDDmadrs[traitRows, 1]

names(MDDmadrs1) <- c("apparent sadness", "concentration difficulties", "inability to feel", 
                      "inner tension", "lassitude", "pessimistic thoughts", "reduced appetite", 
                      "reduced sleep", "reported sadness", "suicidal thoughts", "overall severity")

## visualize how the clinical traits relate to the sample dendrogram
sampleTree2 = hclust(dist(MDDinfData5), method = "average")
traitColors = numbers2colors(MDDmadrs1, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(MDDmadrs1), main = "Sample dendrogram and trait heatmap")

## save
save(MDDinfData5, MDDmadrs1, file = "MDDinfMadrs.RData")
