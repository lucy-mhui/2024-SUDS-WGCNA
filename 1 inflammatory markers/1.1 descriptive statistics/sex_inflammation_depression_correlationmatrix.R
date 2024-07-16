# following code is for MDD - adjust as needed for Control

#########################
# data frame processing #
#########################

# load libraries
library(WGCNA)
library(visdat)
library(mice)
library(tibble)

# load data frames
MDDinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
MDDsex = read.csv("/.../CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv")

# merge data frames 
merge = merge(MDDinfData, MDDsex, on='SUBJLABEL', how='inner')
MDDsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
MDDsexinf = MDDsexinf[1:30]
rownames(MDDsexinf) = MDDsex$SUBJLABEL

# remove missing values 
gsg = goodSamplesGenes(MDDinfM, verbose=3)
if (!gsg$allOK) {
  MDDinfM <- MDDinfM[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(MDDinfF, verbose=3)
if (!gsg$allOK) {
  MDDinfF <- MDDinfF[gsg$goodSamples, gsg$goodGenes]
}

# remove participant outlier 
sampleTree = hclust(dist(MDDinfM), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
MDDinfM = MDDinfM[keepSamples, ]
nGenes = ncol(MDDinfM)
nSamples = nrow(MDDinfM)

sampleTree = hclust(dist(MDDinfF), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
MDDinfF = MDDinfF[keepSamples, ]
nGenes = ncol(MDDinfF)
nSamples = nrow(MDDinfF)

# Multivariate Imputation by Chained Equations
impute = mice(MDDinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfM = complete(impute)

impute = mice(MDDinfF, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfF = complete(impute)




#########
# madrs #
#########

# load madrs data
lnames=load(file="/.../MDDinfMadrs.RData")

# subset by male and female
overlapping_ids <- intersect(row.names(MDDmadrs1), row.names(MDDinfF))
MDDmadrsF <- MDDmadrs1[overlapping_ids, ]

overlapping_ids <- intersect(row.names(MDDmadrs1), row.names(MDDinfM))
MDDmadrsM <- MDDmadrs1[overlapping_ids, ]

# correlation matrix for female
cor_matrices <- lapply(names(MDDinfF), function(col) {
  combined_data <- cbind(MDDmadrsF, MDDinfF[col])
  cor_matrix <- cor(combined_data)
  return(cor_matrix)
})
names(cor_matrices) <- names(MDDinfF)

for (name in names(cor_matrices)) {
  data <- cor_matrices[[name]]
  csv_file <- file.path("/Users/...", paste0(name, ".csv"))
  write.csv(data, file = csv_file, row.names = FALSE)
}

# correlation matrix for male
cor_matrices <- lapply(names(MDDinfM), function(col) {
  combined_data <- cbind(MDDmadrsM, MDDinfM[col])
  cor_matrix <- cor(combined_data)
  return(cor_matrix)
})
names(cor_matrices) <- names(MDDinfM)

for (name in names(cor_matrices)) {
  data <- cor_matrices[[name]]
  csv_file <- file.path("/Users/...", paste0(name, ".csv"))
  write.csv(data, file = csv_file, row.names = FALSE)
}




###########################
# load correlation matrix #
###########################

# load library
library(pheatmap)

# load csv
corrM = read.csv("/Users/.../male_correlation.csv")
corrF = read.csv("/Users/.../female_correlation.csv")

# Set row names to the first column values and remove the first column
rownames(corrM) <- corrM$domains
corrM <- corrM[, -1]

# Visualize the dataframe as a heatmap
pheatmap(corrM, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

# Set row names to the first column values and remove the first column
rownames(corrF) <- corrF$domains
corrF <- corrF[, -1]

# Visualize the dataframe as a heatmap
pheatmap(corrF, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = TRUE, 
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
