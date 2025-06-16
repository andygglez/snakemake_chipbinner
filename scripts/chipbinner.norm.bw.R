suppressWarnings({
  suppressPackageStartupMessages({

  library(ChIPbinner)
  library(optparse)

  })
})
option_list = list(
  make_option(c("--binned_bed"), type = "character", help = "Path to the count matrix"),
  make_option(c("--outdir"), type = "character", help = "Path to the count matrix"),
  make_option(c("--organism"), type = "character", help = "Path to the column Info file"),
  make_option(c("--MS_factor"), type = "numeric", help = "Path to the count matrix"),
  make_option(c("--use_input"), type = "logical", help = "Path to the count matrix"),
  make_option(c("--input_bed"), type = "character", help = "Path to the count matrix")
  # make_option(c("--blacklist"), type = "character", help = "Path to the column Info file"),
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))

# calculate chip_samp_library_size
# sample <- read.csv(opt$binned_bed, sep="\t", header=FALSE)
# sample_signal <- sum(sample$V4)

# blacklist <- read.csv(opt$blacklist, sep="\t", header=FALSE)
# blacklist_signal <- sum(blacklist$V4)

# library_size <- sample_signal - blacklist_signal

if(!opt$use_input){

  print("Normalizing bw using only MS factor")
  norm_bw(out_dir = opt$outdir,
              genome_assembly = opt$organism,
              immunoprecipitated_binned_file = opt$binned_bed,
              use_input = FALSE,
              depth_norm = TRUE,
              raw_count_cutoff = 100,
              pseudocount = 1e-15,
              scaling_factor = opt$MS_factor)

}

if (opt$use_input){


  print("Normalizing bw using MS factor and input file")
  norm_bw(out_dir = opt$outdir,
              genome_assembly = opt$organism,
              immunoprecipitated_binned_file = opt$binned_bed,
              use_input = TRUE,
              input_binned_file = opt$input_bed,
              depth_norm = TRUE,
              raw_count_cutoff = 100,
              pseudocount = 1e-15,
              scaling_factor = opt$MS_factor)

}
