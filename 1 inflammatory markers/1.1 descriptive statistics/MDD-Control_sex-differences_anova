#################
# MDD load data #
#################

# libraries
library(WGCNA)
library(visdat)
library(mice)
library(tibble)

# dataframes
MDDinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_TRTMT.csv")
MDDsex = read.csv("/.../CBN01_DEMO_DATA_Z3_01_V01_TRTMT.csv")

# merge data frames
merge = merge(MDDinfData, MDDsex, on='SUBJLABEL', how='inner')
MDDsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
MDDsexinf = MDDsexinf[1:30]
rownames(MDDsexinf) = MDDsex$SUBJLABEL

# subset for male and female
MDDinfM <- MDDsexinf[MDDsexinf$SEX == 2, -30]
MDDinfF <- MDDsexinf[MDDsexinf$SEX == 1, -30]

# remove missing value
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

# mice
impute = mice(MDDinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfM = complete(impute)

impute = mice(MDDinfF, m = 5, method = 'pmm', maxit = 50, seed = 500)
MDDinfF = complete(impute)




#####################
# Control load data #
#####################

# dataframes
ControlinfData = read.csv("/.../CBN01_CYTORAWL_DATA_Z3_01_V01_CNTRL.csv")
Controlsex = read.csv("/.../CBN01_DEMO_DATA_Z3_01_V01_CNTRL.csv")

# merge data frames 
merge = merge(ControlinfData, Controlsex, on='SUBJLABEL', how='inner')
Controlsexinf = as.data.frame(merge[,-c(1:4,34:35,37:50)])
Controlsexinf = Controlsexinf[1:30]
rownames(Controlsexinf) = Controlsex$SUBJLABEL

# subset for male and female
ControlinfM <- Controlsexinf[Controlsexinf$SEX == 2, -30]
ControlinfF <- Controlsexinf[Controlsexinf$SEX == 1, -30]

# remove missing value
gsg = goodSamplesGenes(ControlinfM, verbose=3)
if (!gsg$allOK) {
  ControlinfM <- ControlinfM[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(ControlinfF, verbose=3)
if (!gsg$allOK) {
  ControlinfF <- ControlinfF[gsg$goodSamples, gsg$goodGenes]
}

# remove participant outlier
sampleTree = hclust(dist(ControlinfM), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
ControlinfM = ControlinfM[keepSamples, ]
nGenes = ncol(ControlinfM)
nSamples = nrow(ControlinfM)

sampleTree = hclust(dist(ControlinfF), method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 4000, minSize = 10)
keepSamples = (clust==1)
ControlinfF = ControlinfF[keepSamples, ]
nGenes = ncol(ControlinfF)
nSamples = nrow(ControlinfF)

# mice
impute = mice(ControlinfM, m = 5, method = 'pmm', maxit = 50, seed = 500)
ControlinfM = complete(impute)

impute = mice(ControlinfF, m = 5, method = 'pmm', maxit = 50, seed = 500)
ControlinfF = complete(impute)




###############################################
# create separate data frames for each marker #
###############################################

# ControlinfF
# Extract column names
column_names <- names(ControlinfF)
# Loop through each column to create individual data frames
for (col_name in column_names) {
  # Create new data frame with column and "sex" column
  new_df <- data.frame(
    value = ControlinfF[[col_name]],
    sex = "female",
    group = "control"
  )
  # Assign the new data frame to a variable with prefix "Control_"
  assign(paste0("ControlF_", col_name), new_df)
}

# ControlinfM
# Extract column names
column_names <- names(ControlinfM)
# Loop through each column to create individual data frames
for (col_name in column_names) {
  # Create new data frame with column and "sex" column
  new_df <- data.frame(
    value = ControlinfM[[col_name]],
    sex = "male",
    group = "control"
  ) 
  # Assign the new data frame to a variable with prefix "Control_"
  assign(paste0("ControlM_", col_name), new_df)
}

# MDDinfF
# Extract column names
column_names <- names(MDDinfF)
# Loop through each column to create individual data frames
for (col_name in column_names) {
  # Create new data frame with column and "sex" column
  new_df <- data.frame(
    value = MDDinfF[[col_name]],
    sex = "female",
    group = "MDD"
  )
  # Assign the new data frame to a variable with prefix "Control_"
  assign(paste0("MDDF_", col_name), new_df)
}

# MDDinfM
# Extract column names
column_names <- names(MDDinfM)
# Loop through each column to create individual data frames
for (col_name in column_names) {
  # Create new data frame with column and "sex" column
  new_df <- data.frame(
    value = MDDinfM[[col_name]],
    sex = "male",
    group = "MDD"
  )
  # Assign the new data frame to a variable with prefix "Control_"
  assign(paste0("MDDM_", col_name), new_df)
}




##########################################
# combine separate data frames by marker #
##########################################

# Initialize a list to store combined data frames
combined_dataframes <- vector("list", 29)

# List of unique suffixes to identify and combine
unique_suffixes <- unique(sub(".*_", "", ls(pattern = ".*_")))

# Loop through each unique suffix
for (i in seq_along(unique_suffixes)) {
  # Get all data frames with the current suffix
  relevant_dataframes <- mget(ls(pattern = paste0(".*_", unique_suffixes[i])))

  # Combine relevant data frames vertically into one
  combined_df <- do.call(rbind, relevant_dataframes)
  
  # Assign the combined data frame to a new variable name
  assign(unique_suffixes[i], combined_df)
}




#############################
# convert levels to factors #
#############################

column_names <- colnames(MDDinfM)

# Loop through each data frame name
for (df_name in column_names) {
  # Extract the current data frame
  current_df <- get(df_name)
  
  # Convert 'sex' and 'group' columns to factors if they exist
  if ("sex" %in% colnames(current_df)) {
    current_df$sex <- factor(current_df$sex)
  }
  if ("group" %in% colnames(current_df)) {
    current_df$group <- factor(current_df$group)
  }
  
  # Assign modified data frame back to its original name in the global environment
  assign(df_name, current_df, envir = .GlobalEnv)
}




#################
# two way anova #
#################

column_names <- colnames(MDDinfM)

# Function to conduct two-way ANOVA
conduct_anova <- function(df) {
  # Check if required columns exist
  if ("sex" %in% colnames(df) && "group" %in% colnames(df) && "value" %in% colnames(df)) {
    # Perform two-way ANOVA
    anova_result <- aov(value ~ sex * group, data = df)
    return(summary(anova_result))  # Return summary of ANOVA
  } else {
    # Print message if columns are missing
    print(paste("Columns 'sex', 'group', or 'value' missing in", deparse(substitute(df))))
    return(NULL)  # Return NULL if columns are missing
  }
}

# Loop through each data frame name
for (df_name in column_names) {
  # Extract the data frame
  df <- get(df_name)
  
  # Conduct ANOVA if columns exist
  anova_summary <- conduct_anova(df)
  
  # Print results if ANOVA was conducted
  if (!is.null(anova_summary)) {
    cat("ANOVA results for", df_name, ":\n")
    print(anova_summary)
    cat("\n")
  }
}
