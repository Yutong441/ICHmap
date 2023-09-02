# summarize connectome and demographics information
library(tidyverse)
library(data.table)

# prepare the data as: ID, metric, atlas, value, group
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
save_dir <- args[2]

patients <- list.files(root, full.names=T)
patients_brief <- list.files(root, full.names=F)

all_df <- as.list(patients) %>% lapply(function(x) {
    if (file.exists(paste(x, "/graph_metrics.csv", sep=""))) {
        one_df <- fread(paste(x, "/graph_metrics.csv", sep="")) %>%
            filter(measure == "count") %>%
            rename(Metric = V1) %>%
            rename(Value = unnormalized) %>%
            rename(Atlas = atlas) %>%
            select(Metric, Atlas, Value)
        one_df$ID <- basename(x)
        return(one_df)
    }
})

all_df <- do.call(rbind, all_df)
fwrite(all_df, paste(save_dir, "all_graph.csv", sep="/"))

# tract
all_df <- as.list(patients) %>% lapply(function(x) {
    if (file.exists(paste(x, "/metric/tract_num.csv", sep=""))) {
        one_df <- fread(paste(x, "/metric/tract_num.csv", sep="")) %>%
            rename(Tract = V1) %>%
            mutate(Tract = gsub(".trk", "", Tract)) %>%
            rename(Num = num)
        one_df$ID <- basename(x)
        return(one_df)
    }
})

all_df <- do.call(rbind, all_df)
fwrite(all_df, paste(save_dir, "all_tract.csv", sep="/"))
