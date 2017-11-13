### Set working directory
#setwd("C:/My work/")

### Viral genome-pc CSV file
# The CSV file MUST contain the headers: contig_id and pc_id
# Please use blanks for unhitted PCs and they will be counted as singletons
pc_csv_file  <- "pc2.csv"
#pc_csv_file  <- "pc.csv"

### Threshold for Sig. value
thresh       <- 2.0

### Number of randomization in permutation test
num_random   <- 0

source("permutation.R")
permutation(pc_csv_file,thresh,num_random)
