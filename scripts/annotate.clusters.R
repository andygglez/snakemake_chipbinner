suppressWarnings({
  suppressPackageStartupMessages({

library(ChIPbinner)
library(optparse)
library(yaml)
library(rtracklayer)
library(dplyr)

  })
})

option_list = list(
  make_option(c("--matrix"), type = "character", help = ""),
  make_option(c("--pooled_bed_file"), type = "character", help = ""),
  make_option(c("--hdbscan_output_file"), type = "character", help = ""),
  make_option(c("--out_dir"), type = "character", help = "")
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))

matrix <- strsplit(opt$matrix, "/")[[1]][4]

comparison <- paste0(strsplit(matrix, "_")[[1]][1], "_", strsplit(matrix, "_")[[1]][2])
mark <- strsplit(matrix, "_")[[1]][3]

print(paste0("Plotting clusters for comparison: ", comparison, " and mark: ", mark))

hdb_file <- strsplit(opt$hdbscan_output_file, "/")[[1]][5]
param1 <- strsplit(hdb_file, "\\.")[[1]][3]
param2 <- strsplit(hdb_file, "\\.")[[1]][4]


hdb_info <- read.csv(opt$hdbscan_output_file)
print(head(hdb_info))

print("eval:")
print(length(unique(hdb_info$X.1)))

if (length(unique(hdb_info$X.1))<=1){
  file.create(paste0(opt$out_dir, "/", comparison,"_",mark,".2clusters.", param1,".", param2, ".annotated_clusters.rda"))
  file.create(paste0(opt$out_dir, "/", comparison,"_",mark,".3clusters.", param1,".", param2, ".annotated_clusters.rda"))
} else {

  tryCatch({
    annotate_clust(
    number_of_clusters = 2,
    matrix_file = opt$matrix,
    pooled_bed_file = opt$pooled_bed_file,
    hdbscan_output_file = opt$hdbscan_output_file,
    output_filename = paste0(comparison, "_", mark, ".2clusters.",param1,".", param2),
    out_dir = opt$out_dir
  )

  annotate_clust(
    number_of_clusters = 3,
    matrix_file = opt$matrix,
    pooled_bed_file = opt$pooled_bed_file,
    hdbscan_output_file = opt$hdbscan_output_file,
    output_filename = paste0(comparison, "_", mark, ".3clusters.",param1,".", param2),
    out_dir = opt$out_dir
  )
  }, warning=function(w){
    print(paste("Problem detected for file", hdb_file))
  }, error=function(e){
    file.create(paste0(opt$out_dir, "/", comparison,"_",mark,".2clusters.", param1,".", param2, ".annotated_clusters.rda"))
    file.create(paste0(opt$out_dir, "/", comparison,"_",mark,".3clusters.", param1,".", param2, ".annotated_clusters.rda"))
  })
}
