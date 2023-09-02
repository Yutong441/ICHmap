# https://github.com/dorianps/LESYMAP/wiki/Testing-LESYMAP-installation
library(LESYMAP)
args <- commandArgs(trailingOnly = TRUE)

save_path <- args[1]
save_root <- args[2]
cogn <- args[3]

save_dir <- paste(save_root, cogn, sep="/")
if (!dir.exists(save_dir)) {
    dir.create(save_dir)
}

filenames <- data.table::fread(paste(save_path, "/", cogn, "_img.txt", sep=""),
                               header=FALSE) %>% dplyr::pull(V1)
behavior <- paste(save_path, "/", cogn, "_behavior.txt", sep="")

lsm = lesymap(filenames, behavior, method = 'sccan', 
              minSubjectPerVoxel = "2.5%", 
              optimizeSparseness=TRUE, saveDir=save_dir)
