suppressWarnings({
  suppressPackageStartupMessages({

  libraries <- c(
    "ChIPbinner",
    "optparse",
    "yaml",
    "rtracklayer",
    "dplyr",
    "reticulate"
  )

  lapply(libraries, require, character.only = TRUE)

  })
})

option_list = list(
  make_option(c("--bigWigWT"), type="character", help="File path to bigWig1"),
  make_option(c("--bigWigTest"), type="character", help="File path to bigWig2"),
  make_option(c("--mark"), type="character", help="Name of the histone mark"),
  make_option(c("--pos_max"), type="numeric", help="Max limit for x and y axes"),
  make_option(c("--pos_min"), type="numeric", help="Max limit for x and y axes"),
  make_option(c("--organism"), type="character", help="Max limit for x and y axes"),
  make_option(c("--wt_condition"), type="character", help="Max limit for x and y axes"),
  make_option(c("--test_condition"), type="character", help="Max limit for x and y axes"),
  make_option(c("--pow"), type="numeric", help="pow argument"),
  make_option(c("--output_folder"), type="character", help="Max limit for x and y axes")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

############################## Find arguments in yaml

config_yaml <- yaml.load_file("config/config.yaml")
if (length(config_yaml$comparisons)==1){
  config_yaml$comparisons <- list(config_yaml$comparisons)
}

############################## extract {comparison} and {mark}

ctrl <- basename(opt$bigWigWT)
ctrl_condition <- strsplit(ctrl, "\\.")[[1]][1]

treatment <- basename(opt$bigWigTest)
treatment_condition <- strsplit(treatment, "\\.")[[1]][1]

comparison <- paste0(treatment_condition,"_",ctrl_condition)

mark <- strsplit(treatment, "\\.")[[1]][2]

yaml_entry <- paste0(comparison,"_",mark)

# Busca si esta la comparacion
if (!("comparison_parameters" %in% names(config_yaml)) || !(yaml_entry %in% names(config_yaml[["comparison_parameters"]]))) {
# Si la comparacion no esta, crea los valores y escribelos
  
  ctrl <- rtracklayer::import.bw(opt$bigWigWT)
  min_signal_ctrl <- min(ctrl$score)
  max_signal_ctrl <- max(ctrl$score)
  
  treatment <- rtracklayer::import.bw(opt$bigWigTest)
  min_signal_treatment <- min(treatment$score)
  max_signal_treatment <- max(treatment$score)
  
  min_signal <- min(c(min_signal_ctrl, min_signal_treatment))
  max_signal <- max(c(max_signal_ctrl, max_signal_treatment))
  
  pos_min <- min_signal
  pos_max <- max_signal
  pow <- 0.75
  
  parameters <- paste0(round(pos_min, digits=2),"_",round(pos_max, digits=2),"_",pow)
  
  config_yaml[['comparison_parameters']] <- c(config_yaml[['comparison_parameters']],   
                                              setNames(list(parameters), yaml_entry))
  
  write_yaml(config_yaml, "config/config.yaml")
  
} else if (yaml_entry %in% names(config_yaml[["comparison_parameters"]])){
  
  config_yaml <- yaml.load_file("config/config.yaml")
  
  values <- config_yaml[["comparison_parameters"]][[yaml_entry]]
  values <- strsplit(values, "_")[[1]]
  
  pos_min <- as.numeric(values[1])
  pos_max <- as.numeric(values[2])
  pow <- as.numeric(values[3])
  
}
print(paste0("Pos_min taken: ", pos_min))
print(paste0("Pos_max taken: ", pos_max))
print(paste0("Pow taken: ", pow))

out_dir = paste0(opt$output_folder, "/", opt$test_condition,"_",opt$wt_condition)
print("Create")
print(out_dir)


if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


genic_intergenic_scatterplot(out_dir = out_dir, genome_assembly = opt$organism,
    cell_line = " ", wildtype_samp_norm_bw = opt$bigWigWT,
    treated_samp_norm_bw = opt$bigWigTest,
    are_R_objects = FALSE, histone_mark = opt$mark, output_filename = paste0(opt$test_condition,"_",opt$wt_condition,"_",opt$mark),
    title_of_plot = paste0("Change in ",opt$mark), xaxis_label = opt$wt_condition,
    yaxis_label = opt$test_condition, max_x = pos_max, max_y = pos_max, min_x = pos_min,
    min_y = pos_min, pow = pow, show_scales = TRUE, show_legend = FALSE,
    legend_pos = "right")

############################################################## Clustering ###############################################################

pre_clust(
  out_dir = out_dir,
  treated_samp_norm_bw = opt$bigWigTest,
  wildtype_samp_norm_bw = opt$bigWigWT,
  output_filename = paste0(opt$test_condition,"_",opt$wt_condition,"_",opt$mark, "_preclust"),
  are_R_objects = FALSE
)
