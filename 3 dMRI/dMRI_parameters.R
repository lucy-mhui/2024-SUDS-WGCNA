# the following code is for MDD mdt - adjust as needed for fat, f, dtiMD, dtiFA, cdi and Control

#########
# setup #
#########

library(WGCNA)
options(stringsAsFactors=FALSE)
lnames=load(file="/.../MDDinfMadrs.RData")
lnames
lnames=load(file="/.../02-MDD-networkConstruction.RData")
lnames
lnames=load(file="/.../1 statistics/MDt.RData")
lnames

# subsetting dataframe
new_row_names <- gsub("^CBN01_", "", rownames(MDDinfData5))
rownames(MDDinfData5) <- new_row_names

rownames(MDD_MD_t_mean) <- MDD_MD_t_mean$SUBJ
MDD_MD_t_mean <- MDD_MD_t_mean[, -1]
new_row_names <- gsub("_01$", "", rownames(MDD_MD_t_mean), perl = TRUE)
rownames(MDD_MD_t_mean) <- new_row_names

# subset common subject IDs
common_subject_ids <- intersect(rownames(MDDinfData5), rownames(MDD_MD_t_mean))
MDDinfData5_common <- MDDinfData5[common_subject_ids, ]
MDD_MD_t_mean_common <- MDD_MD_t_mean[common_subject_ids, ]

# dimensions
dim(MDDinfData5_common)
dim(MDD_MD_t_mean_common)




#########################################
# quantifying module trait associations #
#########################################

# define numbers of genes and samples
nMarkers=ncol(MDDinfData5_common)
nSamples=nrow(MDDinfData5_common)

# recalculate MEs with color labels
MEs0=moduleEigengenes(MDDinfData5_common,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,MDD_MD_t_mean_common,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

# display correlations and p values
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")

# Create textMatrix with asterisks for significant correlations
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(12, 1, 1, 1))
par(las = 0)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_mean_common),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               xLabelsAngle = 80,
               plotLegend = FALSE)





################################################
# Blank: quantifying module trait associations #
################################################

# define numbers of genes and samples
nMarkers=ncol(MDDinfData5_common)
nSamples=nrow(MDDinfData5_common)

# recalculate MEs with color labels
MEs0=moduleEigengenes(MDDinfData5_common,moduleColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,MDD_MD_t_mean_common,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

# display correlations and p values
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")

# Create textMatrix with asterisks for significant correlations
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(0, 1, 0, 0))
par(las = 0)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_mean_common),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               xLabelsAngle = 80,
               plotLegend = FALSE)






##################################################
# top 10%: quantifying module trait associations #
##################################################

# Calculate the thresholds for top 10% and bottom 10%
top10_threshold <- quantile(moduleTraitCor, probs = 0.9)
bottom10_threshold <- quantile(moduleTraitCor, probs = 0.1)

# Display correlations and p-values
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")

# Create textMatrix for significant correlations in top 10% and bottom 10%
textMatrix <- ifelse(moduleTraitCor >= top10_threshold, paste0(signif(moduleTraitCor, 1), "\n", stars, "↑"), 
                     ifelse(moduleTraitCor <= bottom10_threshold, paste0(signif(moduleTraitCor, 1), "\n", stars, "↓"), 
                            ""))
# Ensure dimensions match
dim(textMatrix) <- dim(moduleTraitCor)

# Plot the labeled heatmap
par(mar = c(12, 1, 1, 1))
par(las = 0)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_mean_common),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               xLabelsAngle = 80,
               plotLegend = FALSE)




############################################################
# significant & 10%: quantifying module trait associations #
############################################################

# Calculate the threshold for top 10% and bottom 10%
top10_threshold <- quantile(moduleTraitCor, probs = 0.9)
bottom10_threshold <- quantile(moduleTraitCor, probs = 0.1)

# Create textMatrix for significant correlations in top 10% and bottom 10%, only if p-value is significant
textMatrix <- ifelse(moduleTraitPvalue < 0.05,
                     ifelse(moduleTraitCor >= top10_threshold, paste0(signif(moduleTraitCor, 1), "\n", "↑"), 
                            ifelse(moduleTraitCor <= bottom10_threshold, paste0(signif(moduleTraitCor, 1), "\n", "↓"), 
                                   "")),
                     "")

# Ensure dimensions match
dim(textMatrix) <- dim(moduleTraitCor)

# Plot the labeled heatmap
par(mar = c(12, 1, 1, 1))
par(las = 0)

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(MDD_MD_t_mean_common),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  xLabelsAngle = 80,
  plotLegend = FALSE
)





####################################################################################
# significant & 10% & no background colours: quantifying module trait associations #
####################################################################################

# Define numbers of genes and samples
nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_MD_t_mean_common, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Determine the top 10% and bottom 10% correlation thresholds
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)

# Create a mask for the top 10% and bottom 10% correlation values
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10

# Create a mask for significant p-values
significant_pvalue_mask <- moduleTraitPvalue < 0.05

# Combine the masks
combined_mask <- correlation_mask & significant_pvalue_mask

# Prepare the text matrix with asterisks for significant correlations
stars <- ifelse(moduleTraitPvalue < 0.05, "*", "")
textMatrix <- paste(signif(moduleTraitCor, 1), "\n", stars)
dim(textMatrix) <- dim(moduleTraitCor)

# Mask non-significant values
textMatrix[!combined_mask] <- ""

# Mask non-significant correlations for the heatmap
moduleTraitCor[!combined_mask] <- 0

# Plot the heatmap with the significant correlations only
par(mar = c(12, 1, 1, 1))
par(las = 0)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MDD_MD_t_mean_common),
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





###############################################
## table of module, ROI, correlation, p value #
###############################################

# Define numbers of genes and samples
nMarkers <- ncol(MDDinfData5_common)
nSamples <- nrow(MDDinfData5_common)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(MDDinfData5_common, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, MDD_MD_t_mean_common, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Determine the top 10% and bottom 10% correlation thresholds
threshold_top10 <- quantile(abs(moduleTraitCor), 0.9)
threshold_bottom10 <- quantile(abs(moduleTraitCor), 0.1)

# Create a mask for the top 10% and bottom 10% correlation values
correlation_mask <- abs(moduleTraitCor) >= threshold_top10 | abs(moduleTraitCor) <= threshold_bottom10

# Create a mask for significant p-values
significant_pvalue_mask <- moduleTraitPvalue < 0.05

# Combine the masks
combined_mask <- correlation_mask & significant_pvalue_mask

# Prepare the table with module, y-axis title, correlation, and p-value
module_names <- rownames(moduleTraitCor)
trait_names <- colnames(moduleTraitCor)
filtered_cor <- moduleTraitCor[combined_mask]
filtered_pval <- moduleTraitPvalue[combined_mask]

# Get the indices of the filtered values
filtered_indices <- which(combined_mask, arr.ind = TRUE)

# Create the table
significant_df <- data.frame(
  Module = module_names[filtered_indices[, 1]],
  Y_Axis_Title = trait_names[filtered_indices[, 2]],
  Correlation = filtered_cor,
  P_Value = filtered_pval
)

# Print the table
print(significant_df)





###################################################
# bar graph of positive and negative correlations #
###################################################

# Ensure the 'dplyr' package is installed and loaded
if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

# Create a dataframe categorizing "Module" by positive and negative "Correlation"
module_categorized <- significant_df %>%
  mutate(Correlation_Category = ifelse(Correlation > 0, "Positive", "Negative")) %>%
  group_by(Module, Correlation_Category) %>%
  summarise(Count = n()) %>%
  ungroup()

# Display the categorized dataframe
print(module_categorized)

# Assuming your data is stored in module_categorized dataframe
# First, ensure the column Correlation_Category is a factor with appropriate levels
module_categorized$Correlation_Category <- factor(module_categorized$Correlation_Category,
                                                  levels = c("Negative", "Positive"))

# Plotting a bar graph with colors based on Correlation_Category
library(ggplot2)

# Basic bar plot with custom colors
ggplot(module_categorized, aes(x = Module, y = Count, fill = Correlation_Category)) +
  geom_bar(stat = "identity", position = "stack", width=0.7) +
  scale_fill_manual(values = c("Negative" = "skyblue", "Positive" = "pink")) +
  labs(x = NULL, y = "Count") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  guides(fill = FALSE)




##########################
# count frequency of ROI #
###########################

# Ensure the 'dplyr' package is installed and loaded
if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

# Create a dataframe categorizing "Y_Axis_Title" by positive and negative "Correlation"
y_axis_categorized <- significant_df %>%
  mutate(Correlation_Category = ifelse(Correlation > 0, "Positive", "Negative")) %>%
  group_by(Y_Axis_Title, Correlation_Category) %>%
  summarise(Count = n()) %>%
  ungroup()

# Display the categorized dataframe
print(y_axis_categorized)




########################
# ROIs for each module #
########################

list(y_axis_categorized$Y_Axis_Title[y_axis_categorized$Correlation_Category=="Negative"])
list(y_axis_categorized$Y_Axis_Title[y_axis_categorized$Correlation_Category=="Positive"])
list(significant_df$Y_Axis_Title[significant_df$Module=="MEbrown"])
list(significant_df$Y_Axis_Title[significant_df$Module=="MEgreen"])
list(significant_df$Y_Axis_Title[significant_df$Module=="MEyellow"])
list(significant_df$Y_Axis_Title[significant_df$Module=="MEblue"])
list(significant_df$Y_Axis_Title[significant_df$Module=="MEturquoise"])




################################################
## brain plots: positive/negative correlations #
################################################

# correlation
library(ggplot2)
library(ggseg)

# cortical
library(dplyr)
somedata=tibble(
  region = c("amygdala","corpus callosum","hippocampus","fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital","lingual","parahippocampal","pars triangularis","supramarginal","temporal pole","transverse temporal"),
  p = sample(rep(1, length(region))))

somedata %>%
  ggseg(atlas=dk,
        position = "stacked",
        colour = "black",
        mapping = aes(fill=p),
        show.legend=FALSE)+
  scale_fill_continuous(low="salmon",high="red",na.value="transparent")

# subcortical
somedata %>%
  ggseg(atlas = aseg,
      mapping = aes(fill=p),
      colour="black",
      show.legend=FALSE)+
    scale_fill_continuous(low="salmon",high="red",na.value="transparent")




#########################
## brain plots: modules #
#########################

# change region labels appropriately 
library(ggplot2)
library(ggseg)

# cortical
library(dplyr)
somedata=tibble(
  region = c("hippocampus","amygdala","fusiform","inferior parietal","inferior temporal","isthmus cingulate","lingual","parahippocampal","pars triangularis","supramarginal","temporal pole","transverse temporal"),
  p = sample(rep(1, length(region))))

somedata %>%
  ggseg(atlas=dk,
        position = "stacked",
        colour = "black",
        mapping = aes(fill=p),
        show.legend=FALSE)+
  scale_fill_continuous(low="tan",high="brown",na.value="transparent")

somedata %>%
  ggseg(atlas = aseg,
        position = "stacked",
      mapping = aes(fill=p),
      colour="black",
      show.legend=FALSE)+
    scale_fill_continuous(low="tan",high="brown",na.value="transparent")
