##############
# data input #
##############

# load library
library(readxl)

# load data frame MDD
MDD_cdi_mean = read_excel("/.../MDD_cdi_mean.xlsx", col_names=FALSE)
MDD_dti_FA_mean = read_excel("/.../MDD_dti_FA_mean.xlsx", col_names=FALSE)
MDD_dti_MD_mean = read_excel("/.../MDD_dti_MD_mean.xlsx", col_names=FALSE)
MDD_F_mean = read_excel("/.../MDD_F_mean.xlsx", col_names=FALSE)
MDD_FAt_mean = read_excel("/.../MDD_FAt_mean.xlsx", col_names=FALSE)
MDD_MD_t_mean = read_excel("/.../MDD_MD_t_mean.xlsx", col_names=FALSE)
MDD_subj = read_excel("/.../MDDsubj.xlsx")

# load data frame Control
CONTROL_MD_t_mean = read_excel("/.../CONTROL_MD_t_mean.xlsx", col_names=FALSE)
CONTROL_FAt_mean = read_excel("/.../CONTROL_FAt_mean.xlsx", col_names=FALSE)
CONTROL_F_mean = read_excel("/.../CONTROL_F_mean.xlsx", col_names=FALSE)
CONTROL_dti_MD_mean = read_excel("/.../CONTROL_dti_MD_mean.xlsx", col_names=FALSE)
CONTROL_dti_FA_mean = read_excel("/.../CONTROL_dti_FA_mean.xlsx", col_names=FALSE)
CONTROL_cdi_mean = read_excel("/Users/lucyhui/Downloads/SUDS/MRI/0 data/CONTROL_cdi_mean.xlsx", col_names=FALSE)
CONTROL_subj = read_excel(".../ControlSubj.xlsx")

# as.data.frame
dfs <- list(MDD_cdi_mean, MDD_dti_FA_mean, MDD_dti_MD_mean, MDD_F_mean, MDD_FAt_mean, MDD_MD_t_mean)
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(lapply(dfs[[i]], as.numeric))
}

dfs <- list(CONTROL_MD_t_mean, CONTROL_FAt_mean, CONTROL_F_mean, CONTROL_dti_MD_mean, CONTROL_dti_FA_mean, CONTROL_cdi_mean)
for (i in seq_along(dfs)) {
  dfs[[i]] <- as.data.frame(lapply(dfs[[i]], as.numeric))
}




###################
# file dimensions #
###################
# MDD
dim(MDD_cdi_mean)
dim(MDD_dti_FA_mean)
dim(MDD_dti_MD_mean)
dim(MDD_F_mean)
dim(MDD_FAt_mean)
dim(MDD_MD_t_mean)
dim(MDD_subj)

# Control
dim(CONTROL_MD_t_mean)
dim(CONTROL_FAt_mean)
dim(CONTROL_F_mean)
dim(CONTROL_dti_MD_mean)
dim(CONTROL_dti_FA_mean)
dim(CONTROL_cdi_mean)
dim(CONTROL_subj)




########################
# add ROI column names #
########################

# load library
library(readxl)

# read the Excel file
roi_data <- read_excel("/.../ROIs.xlsx", col_names = FALSE)

# Extract column names from the first column
new_colnames <- roi_data[[1]]

# Prepend "SUBJ" as the first column name
new_colnames <- c(new_colnames)

assign_column_names <- function(df, new_colnames) {
  colnames(df) <- new_colnames
  return(df)
}

# Example: Assign new column names to each dataframe individually
MDD_cdi_mean <- assign_column_names(MDD_cdi_mean, new_colnames)
MDD_dti_FA_mean <- assign_column_names(MDD_dti_FA_mean, new_colnames)
MDD_dti_MD_mean <- assign_column_names(MDD_dti_MD_mean, new_colnames)
MDD_F_mean <- assign_column_names(MDD_F_mean, new_colnames)
MDD_FAt_mean <- assign_column_names(MDD_FAt_mean, new_colnames)
MDD_MD_t_mean <- assign_column_names(MDD_MD_t_mean, new_colnames)
CONTROL_MD_t_mean <- assign_column_names(CONTROL_MD_t_mean, new_colnames)
CONTROL_FAt_mean <- assign_column_names(CONTROL_FAt_mean, new_colnames)
CONTROL_F_mean <- assign_column_names(CONTROL_F_mean, new_colnames)
CONTROL_dti_MD_mean <- assign_column_names(CONTROL_dti_MD_mean, new_colnames)
CONTROL_dti_FA_mean <- assign_column_names(CONTROL_dti_FA_mean, new_colnames)
CONTROL_cdi_mean <- assign_column_names(CONTROL_cdi_mean, new_colnames)




######################
# append subject IDs #
######################

# MDD
dataframes <- list(MDD_cdi_mean, MDD_dti_FA_mean, MDD_dti_MD_mean, MDD_F_mean, MDD_FAt_mean, MDD_MD_t_mean)
df_names <- c("MDD_cdi_mean", "MDD_dti_FA_mean", "MDD_dti_MD_mean", "MDD_F_mean", "MDD_FAt_mean", "MDD_MD_t_mean")

remove_last_row <- function(df) {
  return(df[-nrow(df), ])
}

dataframes <- lapply(dataframes, remove_last_row)

combined_dataframes <- lapply(dataframes, function(df) {
  cbind(MDD_subj, df)
})

for (i in 1:length(df_names)) {
  assign(df_names[i], combined_dataframes[[i]])
}

# Control
dataframes <- list(CONTROL_MD_t_mean, CONTROL_FAt_mean, CONTROL_F_mean, CONTROL_dti_MD_mean, CONTROL_dti_FA_mean, CONTROL_cdi_mean)
df_names <- c("CONTROL_MD_t_mean", "CONTROL_FAt_mean", "CONTROL_F_mean", "CONTROL_dti_MD_mean", "CONTROL_dti_FA_mean", "CONTROL_cdi_mean")

combined_dataframes <- lapply(dataframes, function(df) {
  cbind(CONTROL_subj, df)
})

for (i in 1:length(df_names)) {
  assign(df_names[i], combined_dataframes[[i]])
}




##################################
# remove missing values and -500 #
##################################

library(dplyr)

# dataframes
dfs <- list(MDD_cdi_mean, MDD_dti_FA_mean, MDD_dti_MD_mean, MDD_F_mean, MDD_FAt_mean, MDD_MD_t_mean, CONTROL_MD_t_mean, CONTROL_FAt_mean, CONTROL_F_mean, CONTROL_dti_MD_mean, CONTROL_dti_FA_mean, CONTROL_cdi_mean)

df_names <- c("MDD_cdi_mean", "MDD_dti_FA_mean", "MDD_dti_MD_mean", "MDD_F_mean", "MDD_FAt_mean", "MDD_MD_t_mean", "CONTROL_MD_t_mean", "CONTROL_FAt_mean", "CONTROL_F_mean", "CONTROL_dti_MD_mean", "CONTROL_dti_FA_mean", "CONTROL_cdi_mean")

remove_rows <- function(df) {
  df[df[, 2] != -500, ]
}

for (df_name in df_names) {
  df <- get(df_name)  # Retrieve the data frame by name
  df <- remove_rows(df)  # Remove rows where the first column is -500
  assign(df_name, df)  # Assign the modified data frame back to the original name
}




##################
# missing values #
##################

# List of dataframes
dfs <- list(MDD_cdi_mean, MDD_dti_FA_mean, MDD_dti_MD_mean, MDD_F_mean, MDD_FAt_mean, MDD_MD_t_mean, 
            CONTROL_MD_t_mean, CONTROL_FAt_mean, CONTROL_F_mean, CONTROL_dti_MD_mean, CONTROL_dti_FA_mean, CONTROL_cdi_mean)

# Function to find columns with missing data in a dataframe
find_columns_with_missing_data <- function(df) {
  colnames(df)[colSums(is.na(df)) > 0]
}

# Apply the function to each dataframe in the list
columns_with_missing_data_list <- lapply(dfs, find_columns_with_missing_data)

# Naming the list elements for easier reference
names(columns_with_missing_data_list) <- c("MDD_cdi_mean", "MDD_dti_FA_mean", "MDD_dti_MD_mean", "MDD_F_mean", "MDD_FAt_mean", "MDD_MD_t_mean", 
                                           "CONTROL_MD_t_mean", "CONTROL_FAt_mean", "CONTROL_F_mean", "CONTROL_dti_MD_mean", "CONTROL_dti_FA_mean", "CONTROL_cdi_mean")

# Display list of columns with missing data for each dataframe
columns_with_missing_data_list




###################
# file dimensions #
###################

# Function to remove specified columns from a dataframe
remove_columns <- function(df, columns_to_remove) {
  df_clean <- df[, !colnames(df) %in% columns_to_remove]
  return(df_clean)
}

# Columns to remove
columns_to_remove <- c("'cuneus';", "'wm_entorhinal';")

# Update each dataframe in the list dfs
MDD_cdi_mean <- remove_columns(MDD_cdi_mean, columns_to_remove)
MDD_dti_FA_mean <- remove_columns(MDD_dti_FA_mean, columns_to_remove)
MDD_dti_MD_mean <- remove_columns(MDD_dti_MD_mean, columns_to_remove)
MDD_F_mean <- remove_columns(MDD_F_mean, columns_to_remove)
MDD_FAt_mean <- remove_columns(MDD_FAt_mean, columns_to_remove)
MDD_MD_t_mean <- remove_columns(MDD_MD_t_mean, columns_to_remove)
CONTROL_MD_t_mean <- remove_columns(CONTROL_MD_t_mean, columns_to_remove)
CONTROL_FAt_mean <- remove_columns(CONTROL_FAt_mean, columns_to_remove)
CONTROL_F_mean <- remove_columns(CONTROL_F_mean, columns_to_remove)
CONTROL_dti_MD_mean <- remove_columns(CONTROL_dti_MD_mean, columns_to_remove)
CONTROL_dti_FA_mean <- remove_columns(CONTROL_dti_FA_mean, columns_to_remove)
CONTROL_cdi_mean <- remove_columns(CONTROL_cdi_mean, columns_to_remove)





############################
# table of mean, SD, range #
############################

# MDD
summary_stats = sapply(MDD_cdi_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_cdi_mean_summary = as.data.frame(t(summary_stats))

summary_stats <- sapply(MDD_dti_FA_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_dti_FA_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(MDD_dti_MD_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_dti_MD_mean_summary = as.data.frame(t(summary_stats))

summary_stats <- sapply(MDD_F_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_F_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(MDD_FAt_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_FAt_mean_summary = as.data.frame(t(summary_stats))

summary_stats <- sapply(MDD_MD_t_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
MDD_MD_t_mean_summary = as.data.frame(t(summary_stats))
  
# Control
summary_stats <- sapply(CONTROL_MD_t_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_MD_t_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(CONTROL_FAt_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_FAt_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(CONTROL_F_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_F_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(CONTROL_dti_MD_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_dti_MD_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(CONTROL_dti_FA_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_dti_FA_mean_summary = as.data.frame(t(summary_stats))
  
summary_stats <- sapply(CONTROL_cdi_mean[,-1], function(x) c(mean = mean(x), sd = sd(x), range = diff(range(x))))
CONTROL_cdi_mean_summary = as.data.frame(t(summary_stats))



                        
#############
# bargraphs #
#############

                        library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

# MDD
dfs <- list(MDD_cdi_mean_summary, MDD_dti_FA_mean_summary, MDD_dti_MD_mean_summary, MDD_F_mean_summary, MDD_FAt_mean_summary, MDD_MD_t_mean_summary)

df_names <- c("MDD_cdi_mean_summary", "MDD_dti_FA_mean_summary", "MDD_dti_MD_mean_summary", "MDD_F_mean_summary", "MDD_FAt_mean_summary", "MDD_MD_t_mean_summary")

extract_first_row <- function(df, name) {
  df_first_row <- df[1, , drop = FALSE] # Extract first row, keep it as a dataframe
  df_long <- pivot_longer(df_first_row, cols = everything(), names_to = "stat", values_to = "value")
  df_long$df <- name
  return(df_long)
}

combined_df <- bind_rows(lapply(seq_along(dfs), function(i) extract_first_row(dfs[[i]], df_names[i])))

plot_list <- lapply(df_names, function(df_name) {
  ggplot(combined_df[combined_df$df == df_name, ], aes(x = stat, y = value, fill = stat)) +
    geom_bar(stat = "identity") +
    ggtitle(df_name) +
    theme_minimal()
})

combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)

print(combined_plot)

# Control
dfs <- list(CONTROL_MD_t_mean_summary, CONTROL_FAt_mean_summary, CONTROL_F_mean_summary, CONTROL_dti_MD_mean_summary, CONTROL_dti_FA_mean_summary, CONTROL_cdi_mean_summary)

df_names <- c("CONTROL_MD_t_mean_summary", "CONTROL_FAt_mean_summary", "CONTROL_F_mean_summary", "CONTROL_dti_MD_mean_summary", "CONTROL_dti_FA_mean_summary", "CONTROL_cdi_mean_summary")

extract_first_row <- function(df, name) {
  df_first_row <- df[1, , drop = FALSE] # Extract first row, keep it as a dataframe
  df_long <- pivot_longer(df_first_row, cols = everything(), names_to = "stat", values_to = "value")
  df_long$df <- name
  return(df_long)
}

combined_df <- bind_rows(lapply(seq_along(dfs), function(i) extract_first_row(dfs[[i]], df_names[i])))

plot_list <- lapply(df_names, function(df_name) {
  ggplot(combined_df[combined_df$df == df_name, ], aes(x = stat, y = value, fill = stat)) +
    geom_bar(stat = "identity") +
    ggtitle(df_name) +
    theme_minimal()
})

combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)

print(combined_plot)



                                
#############
# histogram #
#############

library(ggplot2)
library(gridExtra)

# MDD
dfs <- list(MDD_cdi_mean[,2], MDD_dti_FA_mean[,2], MDD_dti_MD_mean[,2], MDD_F_mean[,2], MDD_FAt_mean[,2], MDD_MD_t_mean[,2])
df_names <- c("MDD_cdi_mean", "MDD_dti_FA_mean", "MDD_dti_MD_mean", "MDD_F_mean", "MDD_FAt_mean", "MDD_MD_t_mean")

create_histogram <- function(df, name) {
  ggplot(df, aes(x = df)) +
    geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(name) +
    theme_minimal() +
    labs(x = "Value", y = "Frequency")
}
plot_list <- lapply(seq_along(dfs), function(i) create_histogram(dfs[[i]], df_names[i]))
combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)
print(combined_plot)


# Control
dfs <- list(CONTROL_MD_t_mean[,-1], CONTROL_FAt_mean[,-1], CONTROL_F_mean[,-1], CONTROL_dti_MD_mean[,-1], CONTROL_dti_FA_mean[,-1], CONTROL_cdi_mean[,-1])
df_names <- c("CONTROL_MD_t_mean", "CONTROL_FAt_mean", "CONTROL_F_mean", "CONTROL_dti_MD_mean", "CONTROL_dti_FA_mean", "CONTROL_cdi_mean")

create_histogram <- function(df, name) {
  ggplot(df, aes(x = df[[2]])) +
    geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
    ggtitle(name) +
    theme_minimal() +
    labs(x = "Value", y = "Frequency")
}
plot_list <- lapply(seq_along(dfs), function(i) create_histogram(dfs[[i]], df_names[i]))
combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)
print(combined_plot)

# MDD
par(mfrow = c(2, 3))

hist(MDD_cdi_mean[,2], main = "MDD_cdi_mean", xlab = "Values", ylab = "Frequency")
hist(MDD_dti_FA_mean[,2], main = "MDD_dti_FA_mean", xlab = "Values", ylab = "Frequency")
hist(MDD_dti_MD_mean[,2], main = "MDD_dti_MD_mean", xlab = "Values", ylab = "Frequency")
hist(MDD_F_mean[,2], main = "MDD_F_mean", xlab = "Values", ylab = "Frequency")
hist(MDD_FAt_mean[,2], main = "MDD_FAt_mean", xlab = "Values", ylab = "Frequency")
hist(MDD_MD_t_mean[,2], main = "MDD_MD_t_mean", xlab = "Values", ylab = "Frequency")

# Control
par(mfrow = c(2, 3))

hist(CONTROL_cdi_mean[,2], main = "CONTROL_cdi_mean", xlab = "Values", ylab = "Frequency")
hist(CONTROL_dti_FA_mean[,2], main = "CONTROL_dti_FA_mean", xlab = "Values", ylab = "Frequency")
hist(CONTROL_dti_MD_mean[,2], main = "CONTROL_dti_MD_mean", xlab = "Values", ylab = "Frequency")
hist(CONTROL_F_mean[,2], main = "CONTROL_F_mean", xlab = "Values", ylab = "Frequency")
hist(CONTROL_FAt_mean[,2], main = "CONTROL_FAt_mean", xlab = "Values", ylab = "Frequency")
hist(CONTROL_MD_t_mean[,2], main = "CONTROL_MD_t_mean", xlab = "Values", ylab = "Frequency")



                    
########
# save #
########
save(MDD_dti_MD_mean, CONTROL_dti_MD_mean, file = "dtiMD.RData")
save(MDD_cdi_mean, CONTROL_cdi_mean, file = "cdi.RData")
save(MDD_dti_FA_mean, CONTROL_dti_FA_mean, file = "dtiFA.RData")
save(MDD_F_mean, CONTROL_F_mean, file = "F.RData")
save(MDD_FAt_mean, CONTROL_FAt_mean, file = "FAt.RData")
save(MDD_MD_t_mean, CONTROL_MD_t_mean, file = "MDt.RData")
