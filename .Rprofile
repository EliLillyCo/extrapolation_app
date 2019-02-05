# ensure no symlinks are present in lib-R subdirectory
sfiles <- list.files(paste0("packrat/lib-R/x86_64-pc-linux-gnu/", getRversion()), full.names = TRUE)
if (length(sfiles) > 0L) {
  tmp <- file.remove(sfiles)
    rm(tmp)
    }
    rm(sfiles)

#### -- Packrat Autoloader (version 0.5.0) -- ####
source("packrat/init.R")
#### -- End Packrat Autoloader -- ####
