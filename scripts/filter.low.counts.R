suppressWarnings({
  suppressPackageStartupMessages({

  library(ChIPbinner)
  library(optparse)

  })
})
option_list = list(
  make_option(c("--binned_bed"), type = "character", help = ""),
  make_option(c("--outdir"), type = "character", help = ""),
  make_option(c("--organism"), type = "character", help = ""),
  make_option(c("--cutoff"), type = "numeric", help = "")
)

opt <- parse_args(OptionParser(usage = "Usage: %prog [options]", option_list = option_list))

filter_low_counts(opt$outdir, opt$organism, sample_bedfiles = opt$binned_bed, cutoff = opt$cutoff)