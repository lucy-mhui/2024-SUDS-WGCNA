# the following code is for MDD male - adjust as needed for MDD female and Control male and female

#########################
# data frame processing #
#########################

# load library
library(WGCNA)
library(visdat)
library(mice)

# load data frame
MDDinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
MDDsex = read.csv("/.../CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv")

# combine data frame
merge = merge(MDDinfData, MDDsex, on='SUBJLABEL', how='inner')
MDDsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
MDDsexinf = MDDsexinf[1:30]
rownames(MDDsexinf) = MDDsex$SUBJLABEL

# subset for male
MDDinfM <- MDDsexinf[MDDsexinf$SEX == 2, -30]

# remove missing data 
gsg = goodSamplesGenes(MDDinfM, verbose=3)
vis_miss(MDDinfM)
if (!gsg$allOK) {
  MDDinfM <- MDDinfM[gsg$goodSamples, gsg$goodGenes]
}
vis_miss(MDDinfM)

# outlier participant removal
# perform hierarchial clustering
sampleTree = hclust(dist(MDDinfM), method = "average");

# plot sample tree
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
MDDinfM = MDDinfM[keepSamples, ]
nGenes = ncol(MDDinfM)
nSamples = nrow(MDDinfM)

# visualize
sampleTree_cleaned <- hclust(dist(MDDinfM), method = "average")
plot(sampleTree_cleaned, main = "Sample Clustering After Outlier Removal", 
     sub = "", xlab = "Samples", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# z normalize
MDDinfM=as.data.frame(scale(MDDinfM))

# boxplot
boxplot_data=list()
for(col in colnames(MDDinfM)){
  boxplot_data[[col]]=boxplot(MDDinfM[[col]],plot=FALSE)$out
}
boxplot(boxplot_data,las=2,main="Outliers in each marker", ylab="Values",col="lightblue")

# remove outliers
remove_outliers=function(x){
  mean_val=mean(x,na.rm=TRUE)
  sd_val=sd(x,na.rm=TRUE)
  lower_limit=mean_val-3*sd_val
  upper_limit=mean_val+3*sd_val
  x[which(x<lower_limit | x > upper_limit)]=NA
  return(x)
}
MDDinfM=apply(MDDinfM,2,remove_outliers)
MDDinfM=as.data.frame(MDDinfM)

# visualizing missing data
vis_miss(MDDinfM)

# MICE 
impute = mice(MDDinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfM = complete(impute)
vis_miss(MDDinfM)
