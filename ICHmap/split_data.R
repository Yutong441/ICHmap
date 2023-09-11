library(magrittr)
library(data.table)


all_df <- fread("data/labels/SCCAN.csv")
all_img <- list.files("data/reg_lesion", full.names=F)
img_df <- data.frame(paths=all_img) %>%
    dplyr::mutate(ID = gsub(".nii.gz", "", paths)) %>%
    dplyr::mutate(paths = paste("data/reg_lesion", paths, sep="/")) %>%
    merge(all_df, by="ID")

variables <- colnames(all_df)[colnames(all_df) != "ID"]
for (i in 1:length(variables)) {
    one_var <- variables[i]
    sel_df <- img_df %>%
        dplyr::filter(!is.na(!!as.symbol(one_var))) %>%
        dplyr::select(dplyr::all_of(c(one_var, "paths")))
    sel_df %>% dplyr::select(dplyr::all_of(one_var)) %>%
        write.table(paste("data/tmp/", one_var, "_behavior.txt", sep=""),
                    col.names=F, row.names=F, quote=F)
    sel_df %>% dplyr::select(dplyr::all_of("paths")) %>%
        write.table(paste("data/tmp/", one_var, "_img.txt", sep=""),
                    col.names=F, row.names=F, quote=F)
}
