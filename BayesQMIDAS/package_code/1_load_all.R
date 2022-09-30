##### A script that loads all necessary functions and libraries
##### Daniel Dempsey

### Load in all libraries
# mlVAR is for simulations
# statmod is for Inverse Gaussian
# pROC is for ROC calculations
# tictoc is for benchmarking
# matrixStats is for extra matrix functionality (though I can't remember where exactly it's used)
# mvtnorm is for Multivariate Normal
libs_to_load <- c('mlVAR', 'statmod', 'pROC', 'matrixStats', 'mvtnorm', 'GIGrvg')
lapply(libs_to_load, library, character.only = TRUE)

### Load MIDAS Scripts
MIDAS_files <- list.files('package_code/ALD_MIDAS', pattern = ".R$", full.names = TRUE)
invisible( sapply(MIDAS_files, source) )

### Remove files with no more use
rm(libs_to_load, MIDAS_files)
