# Purpose: Function to merge tables by concatenation.
# Updated: 2020-12-03
library(optparse)

# Options.
opt_list <- list()

opt <- make_option(c("--dir"), type = "character", help = "directory", default = NULL)
opt_list <- c(opt_list, opt)

# Option parsing
parsed_opts <- optparse::OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)
input_dir <- params$dir

# -----------------------------------------------------------------------------
# Merger.
# -----------------------------------------------------------------------------

if (dir.exists(input_dir)) {
  setwd(input_dir)
  
  # Directory contents.
  files <- dir()
  
  # File stems.
  stems <- unique(gsub(pattern = "(.*)_[0-9]+\\.rds", replacement = "\\1", x = files))
  stems <- stems[!grepl(x = stems, pattern = "Master")]
  
  # Loop over different file stems.
  for (stem in stems) {
    
    # Current pattern.
    px <- paste0("^", stem, "_[0-9]+.rds")
    
    # Subset files to merge.
    files_to_merge <- sort(files[grepl(pattern = px, x = files)])
    
    if (length(files_to_merge) > 0) {
      
      # Master file.
      master_file <- paste0(stem, "_Master.rds")
      if (file.exists(master_file)) {
        master <- readRDS(file = master_file)
      } else {
        
        # Set first file as master.
        in_file <- files_to_merge[1]
        master <- readRDS(file = in_file)
        saveRDS(object = master, file = master_file)
        
        # Delete input.
        file.remove(in_file)
        files_to_merge <- files_to_merge[-1]
      }
      
      # Merge.
      for (file in files_to_merge) {
        
        # Import.
        next_file <- readRDS(file = file)
        
        # Merge into master.
        master <- rbind(master, next_file)
        
        # Delete input
        file.remove(file)
      }
      
      ## Export
      saveRDS(object = master, file = master_file)
      }
    } 
  
} else {
  
  warning("Directory does not exist.")
  
}