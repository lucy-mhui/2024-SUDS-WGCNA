# Sample code for MDD - adjust as needed for Control




# data frame processing 

# load libraries
library(WGCNA)
library(visdat)
library(mice)
library(tibble)

# load data frames
MDDinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
MDDsex = read.csv("/.../CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv")

# combine data frames
merge = merge(MDDinfData, MDDsex, on='SUBJLABEL', how='inner')
MDDsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
MDDsexinf = MDDsexinf[1:30]
rownames(MDDsexinf) = MDDsex$SUBJLABEL

# subset data frames for male and female
MDDinfM <- MDDsexinf[MDDsexinf$SEX == 2, -30]
MDDinfF <- MDDsexinf[MDDsexinf$SEX == 1, -30]

# remove missing values for MDDinfM and MDDinfF 
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

# multivariate imputation by chained equations for MDDinfM and MDDinfF 
impute = mice(MDDinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfM = complete(impute)

impute = mice(MDDinfF, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfF = complete(impute)




# prepare data frame for t-test 

# calculate average for inflammatory markers for MDDinfM and MDDinfF 
MDDinfMmean = colMeans(MDDinfM)
MDDinfMmean = as.data.frame(MDDinfMmean)
MDDinfMmean <- tibble::rownames_to_column(MDDinfMmean, "marker")
MDDinfMmean = MDDinfMmean[,-1]

MDDinfFmean = colMeans(MDDinfF)
MDDinfFmean = as.data.frame(MDDinfFmean)
MDDinfFmean <- tibble::rownames_to_column(MDDinfFmean, "marker")
MDDinfFmean = MDDinfFmean[,-1]

# combine MDDinfM and MDDinfF 
marker = matrix(rbind(MDDinfMmean,MDDinfFmean),ncol=1)
marker = as.data.frame(marker)

# add gene names
genes = colnames(MDDinfM)
genes = as.data.frame(genes)
genes = genes[rep(seq_len(nrow(genes)), each = 2), ]
genes = as.data.frame(genes)

# add group labels
group = rep(c("MDD female", "MDD male"),times=29)
group = as.data.frame(group)

# create data frame
sex = data.frame(genes, group, marker)
colnames(sex) = c("marker", "sex", "measurement")




# t-test

# load libraries
library(dplyr)

# initialize data frame to store results
t_test_results <- data.frame(
  marker = character(),
  mean_male = numeric(),
  sd_male = numeric(),
  mean_female = numeric(),
  sd_female = numeric(),
  t_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# list of columns (markers) to test
markers <- colnames(MDDinfM)

# loop through each marker and perform a t-test
for (marker in markers) {
  
  # Extract the measurements for the current marker
  male_marker_data <- MDDinfM[[marker]]
  female_marker_data <- MDDinfF[[marker]]
  
  # Perform t-test
  t_test_result <- t.test(male_marker_data, female_marker_data)
  
  # Calculate means and standard deviations
  mean_male <- mean(male_marker_data, na.rm = TRUE)
  sd_male <- sd(male_marker_data, na.rm = TRUE)
  mean_female <- mean(female_marker_data, na.rm = TRUE)
  sd_female <- sd(female_marker_data, na.rm = TRUE)
  
  # Append the result to the dataframe
  t_test_results <- rbind(t_test_results, data.frame(
    marker = marker,
    mean_male = mean_male,
    sd_male = sd_male,
    mean_female = mean_female,
    sd_female = sd_female,
    t_statistic = t_test_result$statistic,
    p_value = t_test_result$p.value
  ))
}

# print the results table
print(t_test_results)




# bargraph

# load libraries
library(ggplot2)
library(ggpubr)
library(reshape2)

# melt data frame to long format
t_test_long <- melt(t_test_results, id.vars = "marker", measure.vars = c("mean_male", "mean_female"))

# create bar graphs for each marker
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
      geom_signif(comparisons = list(c("mean_male", "mean_female")),
                  annotations = "*", 
                  y_position = max(marker_data$value) + 0.1 * max(marker_data$value),
                  tip_length = 0.01)
  }
  
  # Add plot to list
  plot_list[[i]] <- p
}

# Arrange all plots into one image
combined_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 6)

# Display the combined plot
print(combined_plot)




# Bargraph for differences between average male marker levels and average female marker levels

# Load necessary library
library(dplyr)
library(ggplot2)

# Create a new dataframe with the differences between mean_male and mean_female
difference_df <- t_test_results %>%
  mutate(difference = mean_male - mean_female)

# Create a bar graph for the differences
p <- ggplot(difference_df, aes(x = marker, y = difference)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_text(aes(label = round(difference, 2)), vjust = -0.5) +
  labs(title = "Difference in Mean (MDD Male - MDDFemale) for Inflammatory Markers", x = "Inflammatory Marker", y = "Difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(p)




# Histogram for male and female levels of inflammatory markers 

# Load necessary libraries
library(ggplot2)
library(ggpubr)

# List of columns (markers) to plot
markers <- colnames(MDDinfM)

# Create histograms for each marker
plot_list <- list()
for (i in 1:length(markers)) {
  marker_name <- markers[i]
  
  # Extract the measurements for the current marker
  male_marker_data <- MDDinfM[[marker_name]]
  female_marker_data <- MDDinfF[[marker_name]]
  
  # Create dataframes for males and females
  male_df <- data.frame(measurement = male_marker_data, sex = "Male")
  female_df <- data.frame(measurement = female_marker_data, sex = "Female")
  
  # Combine data into a single dataframe for plotting
  combined_data <- rbind(male_df, female_df)
  
  # Create histogram
  p <- ggplot(combined_data, aes(x = measurement, fill = sex)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
    scale_fill_manual(values = c("Female" = "turquoise", "Male" = "salmon")) + # Flip colors
    labs(title = marker_name, x = "Measurement", y = "Count") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Add plot to list
  plot_list[[i]] <- p
}

# Arrange all plots into one image
combined_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 6)

# Display the combined plot
print(combined_plot)
