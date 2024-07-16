# MDD

## load libraries
library(WGCNA)
library(visdat)
library(mice)
library(tibble)

## load data frames
MDDinfData = read.csv("/Users/lucyhui/Downloads/SUDS/Individual/1 data input/Inflammatory markers/CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
MDDsex = read.csv("/Users/lucyhui/Downloads/SUDS/Individual/1 data input/Sex/CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv")

## combine data frames
merge = merge(MDDinfData, MDDsex, on='SUBJLABEL', how='inner')
MDDsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
MDDsexinf = MDDsexinf[1:30]
rownames(MDDsexinf) = MDDsex$SUBJLABEL

## subset data frames for male and female
MDDinfM <- MDDsexinf[MDDsexinf$SEX == 2, -30]
MDDinfF <- MDDsexinf[MDDsexinf$SEX == 1, -30]

## remove missing values for MDDinfM and MDDinfF 
gsg = goodSamplesGenes(MDDinfM, verbose=3)
if (!gsg$allOK) {
  MDDinfM <- MDDinfM[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(MDDinfF, verbose=3)
if (!gsg$allOK) {
  MDDinfF <- MDDinfF[gsg$goodSamples, gsg$goodGenes]
}

# remove participant outlier for MDDinfM and MDDinfF 
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

## multivariate imputation by chained equations for MDDinfM and MDDinfF 
impute = mice(MDDinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfM = complete(impute)

impute = mice(MDDinfF, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfF = complete(impute)

## calculate average for inflammatory markers for MDDinfM and MDDinfF 
MDDinfMmean = colMeans(MDDinfM)
MDDinfMmean = as.data.frame(MDDinfMmean)
MDDinfMmean <- tibble::rownames_to_column(MDDinfMmean, "marker")
MDDinfMmean = MDDinfMmean[,-1]

MDDinfFmean = colMeans(MDDinfF)
MDDinfFmean = as.data.frame(MDDinfFmean)
MDDinfFmean <- tibble::rownames_to_column(MDDinfFmean, "marker")
MDDinfFmean = MDDinfFmean[,-1]

## combine MDDinfM and MDDinfF 
marker = matrix(rbind(MDDinfMmean,MDDinfFmean),ncol=1)
marker = as.data.frame(marker)

## add gene names
genes = colnames(MDDinfM)
genes = as.data.frame(genes)
genes = genes[rep(seq_len(nrow(genes)), each = 2), ]
genes = as.data.frame(genes)

## add group labels
group = rep(c("MDD female", "MDD male"),times=29)
group = as.data.frame(group)

## create data frame
sex = data.frame(genes, group, marker)
colnames(sex) = c("marker", "sex", "measurement")

# Control
