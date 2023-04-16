source( 'package_code/1_load_all.R' )

setwd( 'package_code/Simulations/IRTS_MIDAS_Sims/Output/' )

stat_fun <- function( x ) {
  
  load( paste0(x, '/res.Rdata') )
  tab <- matrix( '0', 3, 5 )
  for( i in 1:3 ) {
    
    xx <- res[[i]]
    mean_res <- apply( xx[, 1:3], 2, mean )
    se_res <- apply( xx[, 1:3], 2, sd ) / sqrt(1000)
    
    DLF1_mean <- mean( xx[,4] / (-2 * xx[,5]) )
    DLF2_mean <- mean( xx[,6] / (-2 * xx[,7]) )
    
    DLF1_se <- sd( xx[,4] / (-2 * xx[,5]) ) / sqrt(1000)
    DLF2_se <- sd( xx[,6] / (-2 * xx[,7]) ) / sqrt(1000)
    
    res1 <- round(c( mean_res, DLF1 = DLF1_mean, DLF2 = DLF2_mean ), 2)
    res2 <- round(c( se_res, DLF1 = DLF1_se, DLF2 = DLF2_se ), 2)
    
    tab[i, ] <- paste0( '& ', res1, ' \newline (', res2, ')' )
    
  }
  
  rownames( tab ) <- paste0( 'n', c(100, 500, 1000) )
  colnames( tab ) <- names( res1 )
  
  t( tab )
  
}

nms <- substr(list.dirs(recursive = FALSE), 3, 100)

all_res <- lapply( nms, stat_fun )
names( all_res ) <- nms

all_res$Irreg7
all_res$Irreg30
all_res$TW7
all_res$TW25
all_res$Noise

