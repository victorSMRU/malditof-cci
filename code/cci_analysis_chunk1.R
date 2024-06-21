# Load packages
library(data.table)
library(MALDIquant)

# Set working directory
my_wd = 'path_to_root_directory' # this needs editing
setwd(my_wd)

# Load custom function
source("code/cci_algorithm.R")

# Load the processed mass spectra library
msl = readRDS("Rdata/msl_processed.RDS")

# Make an output data frame consisting of all possible pairwise comparisons
output = data.frame(matrix(combn(names(msl), 2), ncol = 2, byrow = T))
colnames(output) = c("spectrum.in", "spectrum.out")

# Subset the first 1 million lines
output = output[1:1000000, ]
output$tag = paste(output$spectrum.in, output$spectrum.out)

# Check if some CCI values were already computed
out_file = 'cci_output/out1.csv'
if (file.exists(out_file)) {
  df_done = read.csv(file = out_file, col.names = c("spectrum.in","spectrum.out","ccf1","ccf2","ccf3","ccf4","ccf5","ccf6","ccf7","ccf8","ccf9","ccf10","ccf11","ccf12","ccf13","ccf14","ccf15","ccf16","ccf17","ccf18","ccf19","ccf20","ccf21","ccf22","ccf23","ccf24","ccf25","ccf26","ccf27","ccf28","ccf29","ccf30","ccf31","ccf32","ccf33","ccf34","ccf35","ccf36"))
  df_done$tag = paste(df_done$spectrum.in, df_done$spectrum.out)
} else {
  df_done = data.frame(tag = NA)
}

# Remove the lines already done
output = output[!output$tag %in% df_done$tag, ]
output$tag = NULL
rm(df_done)

# 
n = nrow(output)
for(i in 1:n) {
  cci.algorithm(msl[output$spectrum.in[i]], msl[output$spectrum.out[i]], outfile = "cci_output/out1.csv", min.mass = 2000, max.mass = 20000, interval = 500)
}
