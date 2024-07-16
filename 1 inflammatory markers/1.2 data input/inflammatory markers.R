# the following code is for MDD - adjust for Control as needed

######################
# process dataframes #
######################

# load package
library(WGCNA)

# load data frame
MDDinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")

# select columns and rename rownames
MDDinfData1 = as.data.frame(MDDinfData[,-c(1:4,34)])
rownames(MDDinfData1) = MDDinfData$SUBJLABEL

# visualize missing data
vis_miss(MDDinfData1)

# remove patients with missing data
gsg = goodSamplesGenes(MDDinfData1, verbose=3)
if (!gsg$allOK) {
  MDDinfData2 <- MDDinfData1[gsg$goodSamples, gsg$goodGenes]
}
#vis_miss(MDDinfData2)

# apply multiple imputation by chained equation MICE method to impute the rest of the missing data
impute = mice(MDDinfData2, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfData3 = complete(impute)
#vis_miss(MDDinfData3)




################################
# sample hierarchal clustering #
################################

# perform hierarchal clustering
sampleTree = hclust(dist(MDDinfData3), method = "average");

# plot sample tree
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# choosing height cut to remove offending sample
keepSamples = (clust==1)
MDDinfData3 = MDDinfData3[keepSamples, ]
nGenes = ncol(MDDinfData3)
nSamples = nrow(MDDinfData3)

# visualize
sampleTree_cleaned <- hclust(dist(MDDinfData3), method = "average")
plot(sampleTree_cleaned, main = "Sample Clustering After Outlier Removal", 
     sub = "", xlab = "Samples", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# z normalize data frame
MDDinfData4=as.data.frame(scale(MDDinfData3))




######################################
# remove inflammatory marker outlier #
######################################

# boxplot
boxplot_data=list()
for(col in colnames(MDDinfData4)){
  boxplot_data[[col]]=boxplot(MDDinfData4[[col]],plot=FALSE)$out
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
MDDinfData5=apply(MDDinfData4,2,remove_outliers)
MDDinfData5=as.data.frame(MDDinfData5)

# visualizing missing data
vis_miss(MDDinfData5)

# mice
impute = mice(MDDinfData5, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfData5 = complete(impute)
vis_miss(MDDinfData5)
