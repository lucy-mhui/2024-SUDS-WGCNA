# load data frames
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/3 sex differences/dtiM.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/3 sex differences/dtiF.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/MDt.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/FAt.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/F.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/dtiMD.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/dtiFA.RData")
lnames=load(file="/Users/lucyhui/Downloads/SUDS/MRI/1 statistics/cdi.RData")

# DTI SUBJ labels
dataframes <- list(MDD_dti_MD_mean, CONTROL_dti_MD_mean,MDD_cdi_mean, CONTROL_cdi_mean,
                   MDD_dti_FA_mean, CONTROL_dti_FA_mean,MDD_F_mean, CONTROL_F_mean,
                   MDD_FAt_mean, CONTROL_FAt_mean, MDD_MD_t_mean, CONTROL_MD_t_mean)
df_names <- c("MDD_dti_MD_mean", "CONTROL_dti_MD_mean", "MDD_cdi_mean", "CONTROL_cdi_mean",
              "MDD_dti_FA_mean", "CONTROL_dti_FA_mean", "MDD_F_mean", "CONTROL_F_mean",
              "MDD_FAt_mean", "CONTROL_FAt_mean", "MDD_MD_t_mean", "CONTROL_MD_t_mean")

process_dataframe <- function(df) {
  new_row_names <- gsub("^CBN01_", "", rownames(df))
  rownames(df) <- new_row_names
  if ("SUBJ" %in% colnames(df)) {
    rownames(df) <- df$SUBJ
    df <- df[, -which(colnames(df) == "SUBJ")]
  }
  new_row_names <- gsub("_01$", "", rownames(df), perl = TRUE)
  rownames(df) <- new_row_names
  
  return(df)
}

modified_dataframes <- lapply(dataframes, process_dataframe)
names(modified_dataframes) <- df_names

MDD_dti_MD_mean <- modified_dataframes[["MDD_dti_MD_mean"]]
CONTROL_dti_MD_mean <- modified_dataframes[["CONTROL_dti_MD_mean"]]
MDD_cdi_mean <- modified_dataframes[["MDD_cdi_mean"]]
CONTROL_cdi_mean <- modified_dataframes[["CONTROL_cdi_mean"]]
MDD_dti_FA_mean <- modified_dataframes[["MDD_dti_FA_mean"]]
CONTROL_dti_FA_mean <- modified_dataframes[["CONTROL_dti_FA_mean"]]
MDD_F_mean <- modified_dataframes[["MDD_F_mean"]]
CONTROL_F_mean <- modified_dataframes[["CONTROL_F_mean"]]
MDD_FAt_mean <- modified_dataframes[["MDD_FAt_mean"]]
CONTROL_FAt_mean <- modified_dataframes[["CONTROL_FAt_mean"]]
MDD_MD_t_mean <- modified_dataframes[["MDD_MD_t_mean"]]
CONTROL_MD_t_mean <- modified_dataframes[["CONTROL_MD_t_mean"]]

# rename ROI
dataframes=list(MDD_dti_MD_mean, CONTROL_dti_MD_mean,MDD_cdi_mean, CONTROL_cdi_mean,
                MDD_dti_FA_mean, CONTROL_dti_FA_mean,MDD_F_mean, CONTROL_F_mean,
                MDD_FAt_mean, CONTROL_FAt_mean, MDD_MD_t_mean, CONTROL_MD_t_mean,
                MDD_cdi_F, MDD_dti_FA_F, MDD_dti_MD_F, MDD_F_F, MDD_FAt_F, MDD_MD_t_F,
                MDD_cdi_M, MDD_dti_FA_M, MDD_dti_MD_M, MDD_F_M, MDD_FAt_M, MDD_MD_t_M)

dataframe_names <- c("MDD_dti_MD_mean", "CONTROL_dti_MD_mean", "MDD_cdi_mean", "CONTROL_cdi_mean",
                     "MDD_dti_FA_mean", "CONTROL_dti_FA_mean", "MDD_F_mean", "CONTROL_F_mean",
                     "MDD_FAt_mean", "CONTROL_FAt_mean", "MDD_MD_t_mean", "CONTROL_MD_t_mean",
                     "MDD_cdi_F", "MDD_dti_FA_F", "MDD_dti_MD_F", "MDD_F_F", "MDD_FAt_F", "MDD_MD_t_F",
                     "MDD_cdi_M", "MDD_dti_FA_M", "MDD_dti_MD_M", "MDD_F_M", "MDD_FAt_M", "MDD_MD_t_M")

columnnames <- c(
  "thalamus", "caudate", "putamen", "pallidum", "brainstem", "hippocampus", "amygdala", "accumbens", 
  "corpus callosum posterior", "corpus callosom mid posterior", "corpus callosom central", 
  "corpus callosom mid anterior", "corpus callosom anterior", "caudal anterior cingulate", 
  "caudal middle frontal", "entorhinal", "fusiform", "inferior parietal", "inferior temporal", 
  "isthmus cingulate", "lateral occipital", "lateral orbitofrontal", "lingual", "medial orbitfrontal", 
  "middle temporal", "parahippocampal", "paracentral", "pars opercularis", "pars orbitalis", 
  "pars triangularis", "pericalcarine", "postcentral", "posterior cingulate", "precentral", "precuneus", 
  "rostral anterior cingulate", "rostral middle frontal", "superior frontal", "superior parietal", 
  "superior temporal", "supramarginal", "frontal pole", "temporal pole", "transverse temporal", "insula", 
  "caudal anterior cingulate wm", "caudal middle frontal wm", "cuneus wm", "fusiform wm", "inferior parietal wm", 
  "inferior temporal wm", "isthmus cingulate wm", "lateral occipital wm", "lateral orbitofrontal wm", 
  "lingual wm", "medial orbitfrontal wm", "middle temporal wm", "parahippocampal wm", "paracentral wm", 
  "pars opercularis wm", "pars orbitalis wm", "pars triangularis wm", "pericalcarine wm", "postcentral wm", 
  "posterior cingulate wm", "precentral wm", "precuneus wm", "rostral anterior cingulate wm", 
  "rostral middle frontal wm", "superior frontal wm", "superior parietal wm", "superior temporal wm", 
  "supramarginal wm", "frontal pole wm", "temporal pole wm", "transverse temporal wm", "insula wm"
)

for (i in seq_along(dataframes)) {
  colnames(dataframes[[i]]) <- columnnames
  assign(dataframe_names[i], dataframes[[i]])
}

# save 
save(MDD_dti_MD_mean, CONTROL_dti_MD_mean,MDD_cdi_mean, CONTROL_cdi_mean,
     MDD_dti_FA_mean, CONTROL_dti_FA_mean,MDD_F_mean, CONTROL_F_mean,
     MDD_FAt_mean, CONTROL_FAt_mean, MDD_MD_t_mean, CONTROL_MD_t_mean,
     MDD_cdi_F, MDD_dti_FA_F, MDD_dti_MD_F, MDD_F_F, MDD_FAt_F, MDD_MD_t_F,
     MDD_cdi_M, MDD_dti_FA_M, MDD_dti_MD_M, MDD_F_M, MDD_FAt_M, MDD_MD_t_M, file = "0_ALL_DTI.RData")
