suppressWarnings({
  suppressPackageStartupMessages({

  library(ChIPbinner)
  library(optparse)
  library(rtracklayer)

  })
})
option_list = list(
  make_option(c("--rep1"), type = "character", help = ""),
  make_option(c("--rep2"), type = "character", help = ""),
  make_option(c("--outdir"), type = "character", help = ""),
  make_option(c("--outname"), type = "character", help = "")
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))

merge_norm_bw(
  rep1 = opt$rep1,
  rep2 = opt$rep2,
  merged_samples_label = opt$outname,
  out_dir = opt$outdir
)

