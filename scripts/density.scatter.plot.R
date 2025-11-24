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
  make_option(c("--out_name"), type = "character", help = ""),
  make_option(c("--out_dir"), type = "character", help = ""),
  make_option(c("--organism"), type = "character", help = ""),
  make_option(c("--bigWigWT"), type = "character", help = ""),
  make_option(c("--bigWigTest"), type = "character", help = ""),
  make_option(c("--mark"), type = "character", help = ""),
  make_option(c("--annotated_clust"), type = "character", help = ""),
  make_option(c("--clust_num"), type = "numeric", help = ""),
  make_option(c("--wt_condition"), type = "character", help = ""),
  make_option(c("--test_condition"), type = "character", help = "")
)


opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))


anno_clus_file <- strsplit(opt$annotated_clust, '\\.')[[1]]
param1 <- strsplit(opt$annotated_clust, '\\.')[[1]][3]
param2 <- strsplit(opt$annotated_clust, '\\.')[[1]][4]

out_name <- paste0(opt$out_name, ".", param1, ".", param2, ".plots")


############################## extract {comparison} and {mark}

ctrl <- basename(opt$bigWigWT)
ctrl_condition <- strsplit(ctrl, "\\.")[[1]][1]

treatment <- basename(opt$bigWigTest)
treatment_condition <- strsplit(treatment, "\\.")[[1]][1]

comparison <- paste0(treatment_condition,"_",ctrl_condition)

mark <- strsplit(treatment, "\\.")[[1]][2]

yaml_entry <- paste0(comparison,"_",mark)

config_yaml <- yaml.load_file("config/config.yaml")
  
values <- config_yaml[["comparison_parameters"]][[yaml_entry]]
print(yaml_entry)
print(values)

values <- strsplit(values, "_")[[1]]

pos_min <- as.numeric(values[1])
pos_max <- as.numeric(values[2])
pow <- as.numeric(values[3])


if ((file.size(opt$annotated_clust)==0) | (! file.exists(opt$annotated_clust))){

  file.name <- paste0(opt$out_dir, out_name, "_density_scatterplot.pdf")
  print(file.name)
  file.create(file.name)

} else {

  density_based_scatterplot(
    out_dir = opt$out_dir,
    genome_assembly = opt$organism,
    are_R_objects = FALSE,
    output_filename = out_name,
    wildtype_samp_norm_bw = opt$bigWigWT,
    treated_samp_norm_bw = opt$bigWigTest,
    cell_line = "",
    histone_mark = opt$mark,
    annotated_clusters = opt$annotated_clust,
    number_of_clusters = opt$clust_num,
    title_of_plot = opt$mark,
    pow = pow,
    min_x = pos_min,
    min_y = pos_min,
    max_x = pos_max,
    max_y = pos_max,
    hexbins = 150,
    show_scales = FALSE,
    xaxis_label = opt$wt_condition,
    yaxis_label = opt$test_condition,
    height_of_figure = 6,
    width_of_figure = 15,
    include_additional_density_plot = TRUE,
    filter_extreme_bins = TRUE
  )

}