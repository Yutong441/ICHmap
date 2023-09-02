library(magrittr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

inp_df_path <- args[1]
save_path <- args[2]

inp_df <- fread(inp_df_path)
ich_vol <- fread("data/stats/vol.csv")
colnames(ich_vol) <- c("ID", "ICH_vol")

inp_df <- inp_df %>% 
    dplyr::mutate(paths=paste("data/reg_lesion/", ID, ".nii.gz", sep=""))

if (!"ICH_vol" %in% colnames(inp_df)) {
    inp_df <- inp_df %>% merge(ich_vol, by="ID")
}

fwrite(inp_df, save_path)
