### Set working directory
#setwd("C:/My work/")

### Viral genome-pc CSV file
# The CSV file MUST contain the headers: contig_id and pc_id
# Please use blanks for unhitted PCs and they will be counted as singletons
pc_csv_file  <- "pc2.csv"

### Threshold for Sig. value
thresh       <- 1.0

### 
sig_csv_file <- "sig2_test.csv"

source("pcprofiles.R")
pcprofiles(pc_csv_file,thresh,sig_csv_file)
