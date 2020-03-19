pkgs2load <- c('dplyr', 'magrittr', 'tidyr', 'flowCore', 'data.table', 'readr')
sapply(pkgs2load, require, character.only = TRUE)
set.seed(20191007)

# directory and files
# data from Flow Repository https://flowrepository.org/, experiment-ID FR-FCM-Z24N
# data_dir <- "Z:/data/RA-HD data/Fcs_raw/"
data_dir <- "DATA/FlowRepository_FR-FCM-Z24N_files copy/UNSTIM/"
files <- list.files(data_dir, pattern = ".fcs")

# markers
lineage_markers <- c("147Sm_CD20", "170Er_CD3", "145Nd_CD4",
                     "146Nd_CD8a", "169Tm_CD45RA", "176Yb_CD56",
                     "148Nd_CD16", "160Gd_CD14", "209Bi_CD61",
                     "159Tb_CD11c", "151Eu_CD123", "174Yb_HLA-DR")

# rename column names
extract_exprs <- function(ff, subsample = 20000){
  exprs <- ff@exprs[sample(nrow(ff@exprs), min(subsample, nrow(ff@exprs))), ]
  name_match <- colnames(exprs) %in% ff@parameters@data$name
  ff@parameters@data$desc[is.na(ff@parameters@data$desc)] <-
    ff@parameters@data$name[is.na(ff@parameters@data$desc)]
  colnames(exprs)[name_match] <- ff@parameters@data$desc[name_match[!is.na(name_match)]]
  return(exprs)
}

# read and prepare data for analysis
dt <- rbindlist(lapply(files, function(file){
  dt_loc <- data.table(extract_exprs(read.FCS(paste0(data_dir, file))),
                       subsample=20000)[, .SD,.SDcols = lineage_markers]
  dt_loc[, id:=file]
  dt_loc
}))

# add ids and grouping information
dt2 <- dt %>%
  mutate(group = ifelse(grepl('KTR', id), 'control', 'diseased'),
         id = as.integer(as.factor(id))) %>%
  select(id, group, everything())

# write to file
write_csv(x = dt2, path = 'data/cell_data.csv')
