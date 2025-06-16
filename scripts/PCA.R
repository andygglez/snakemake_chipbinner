

suppressWarnings({
  suppressPackageStartupMessages({

  library(ChIPbinner)

  })
})




plot_PCA(
  out_dir = "outdir",
  output_filename = "PCA_plot",
  colors = c("blue", "blue3", "red", "red4", "forestgreen", "green1"),
  sample_labels = c("samp_1", "samp_2", "samp_3", "samp_4", "samp_5", "samp_6"),
  plot_title = "",
  plot_height = 8,
  plot_width = 12,
  "samp1_10kb.norm.bw",
  "samp2_10kbb.norm.bw",
  "samp3_10kbb.norm.bw",
  "samp4_10kbb.norm.bw",
  "samp5_10kbb.norm.bw",
  "samp6_10kbb.norm.bw",
)