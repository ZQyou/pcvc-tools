### Set working directory
#setwd("C:/My work/")

### Viral genome-pc CSV file
# The CSV file MUST contain the headers: contig_id,VC,pc_id
# Please use blanks for unhitted PCs and they will be counted as singletons

#pc_csv_file  <- "remark_vc.csv"
#pc_csv_file  <- "Modified_remark_vc.csv"
pc_csv_file  <- "Modified_profile.csv"

### Uncomment if you want to use previous format of CSV file
#source("VCscript.R")
#VCscript(pc_csv_file)

### For new input
source("VCscript_mod.R")
VCscript_mod(pc_csv_file)
