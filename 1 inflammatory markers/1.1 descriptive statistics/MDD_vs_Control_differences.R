#########################
# data frame processing #
#########################

# load dataframes
MDD=read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
Control=read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_CNTRL.csv")

# remove missing value
gsg = goodSamplesGenes(MDD, verbose=3)
if (!gsg$allOK) {
  MDD <- MDD[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(Control, verbose=3)
if (!gsg$allOK) {
  Control <- Control[gsg$goodSamples, gsg$goodGenes]
}

# remove participant outlier
sampleTree = hclust(dist(MDD), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
MDD = MDD[keepSamples, ]
nGenes = ncol(MDD)
nSamples = nrow(MDD)

sampleTree = hclust(dist(Control), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
Control = Control[keepSamples, ]
nGenes = ncol(Control)
nSamples = nrow(Control)

# Multivariate Imputation by Chained Equations
impute = mice(MDD, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDD = complete(impute)

impute = mice(Control, m = 5, method = 'pmm', maxit = 50, seed = 500)
Control = complete(impute)




##########
# t-test #
##########
# load library
library(dplyr)

# initialize a dataframe to store t-test results
t_test_results <- data.frame(
  marker = character(),
  mean_mdd = numeric(),
  sd_mdd = numeric(),
  mean_control = numeric(),
  sd_control = numeric(),
  t_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# List of columns (markers) to test
markers <- colnames(MDD)

# Loop through each marker and perform a t-test
for (marker in markers) {
  # Extract the measurements for the current marker
  MDD_marker_data <- MDD[[marker]]
  Control_marker_data <- Control[[marker]]
  
  # Perform t-test
  t_test_result <- t.test(MDD_marker_data, Control_marker_data)
  
  # Calculate means and standard deviations
  mean_mdd <- mean(MDD_marker_data, na.rm = TRUE)
  sd_mdd <- sd(MDD_marker_data, na.rm = TRUE)
  mean_control <- mean(Control_marker_data, na.rm = TRUE)
  sd_control <- sd(Control_marker_data, na.rm = TRUE)
  
  # Append the result to the dataframe
  t_test_results <- rbind(t_test_results, data.frame(
    marker = marker,
    mean_mdd = mean_mdd,
    sd_mdd = sd_mdd,
    mean_control = mean_control,
    sd_control = sd_control,
    t_statistic = t_test_result$statistic,
    p_value = t_test_result$p.value
  ))
}



#############
# bar graph #
#############

# load library
library(ggplot2)
library(ggpubr)
library(reshape2)

# melt the dataframe to long format for easier plotting with ggplot2
t_test_long <- melt(t_test_results, id.vars = "marker", measure.vars = c("mean_mdd", "mean_control"))

# Create bar graphs for each marker
plot_list <- list()
for (i in 1:length(unique(t_test_long$marker))) {
  marker_name <- unique(t_test_long$marker)[i]
  
  # Filter data for the current marker
  marker_data <- t_test_long[t_test_long$marker == marker_name,]
  
  # Create bar graph
  p <- ggplot(marker_data, aes(x = variable, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = marker_name, x = "Sex", y = "Mean Level of Inflammatory Marker") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text(aes(label = round(value, 2)), vjust = -0.5)
  
  # Add significance annotation if p-value < 0.05
  p_value <- t_test_results[t_test_results$marker == marker_name, "p_value"]
  if (p_value < 0.05) {
    p <- p + 
      geom_signif(comparisons = list(c("mean_mdd", "mean_control")),
                  annotations = "*", 
                  y_position = max(marker_data$value) + 0.1 * max(marker_data$value),
                  tip_length = 0.01)
  }
  
  # Add plot to list
  plot_list[[i]] <- p
}

# Arrange all plots into one image
combined_plot <- ggarrange(plotlist = plot_list, ncol = 6, nrow = 5)
