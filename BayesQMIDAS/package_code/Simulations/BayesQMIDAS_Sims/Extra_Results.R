##### Writing extra results from the simulation study
##### Daniel Dempsey
library( readr )
library( corrplot )
library( xtable )
setwd( 'package_code/Simulations/BayesQMIDAS_Sims/Output/' )

### Colour pallette for plots, taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",  
           "#0072B2", "#D55E00", "#CC79A7", "#999999", '#FFFFFF')
scales::show_col(mycol)

filenames <- list.files()

nms <- c( 'Beta Simulation', 'Hyperprior (1, 1)', 'Hyperprior(3, 3)', 
          'Standard Deviation 3', 'Time Window 25', 'Time Window 7' )
allnms <- rep( nms, 4 )

all_ppd <- list()
for ( nm in filenames ) {
  
  load( paste0( nm, '/res.Rdata' ) )
  dats <- data.frame( training = all_res$ppd_train, test = all_res$ppd_test )
  all_ppd[[nm]] <- dats
  write_csv( dats, file = paste0(nm, '/ppd_res.csv') )
  
}

mean_ppd <- lapply( all_ppd, function(x){ sapply( x, mean ) } )
ppd_tab <- data.frame( Experiment = allnms, do.call( 'rbind', mean_ppd )*100 )
xtable( ppd_tab, digits = 0 )

### Simulation Correlation matrices
colfunc <- colorRampPalette(mycol[c(1,9,2)])
load( 'Large2_2_multilevel_33/res.Rdata' )
res1 <- all_res
load( 'Large_2_with_multilevel_33/res.Rdata' )
res2 <- all_res

pdf( '../Extra/sim1_cor.pdf' )
corrplot( res1$data_correlation, col = colfunc(200), tl.col = 'black',
          method = 'number' )
dev.off()
pdf( '../Extra/sim2_cor.pdf' )
corrplot( res2$data_correlation, col = colfunc(200), tl.col = 'black',
          method = 'number' )
dev.off()

#### Variable Selection Boxplots
boxplot_make <- function( x, graph_title ) {
  
  load( paste0(x, '/res.Rdata') )
  MCMC_len <- 80000 
  burn_in <- 40000
  thin <- 10
  keep <- seq( burn_in, MCMC_len, thin )
  
  varsels_list <- lapply( all_res$MCMC, function( x ) { x$vsres[keep, ] } )
  varsels_mean_list <- lapply( varsels_list, function(x) { apply(x, 2, mean) } )
  varsels_mean <- do.call( 'rbind', varsels_mean_list )
  quantlab <- all_res$MCMC[[1]]$quantile
  
  use_col <- c( rep(mycol[2], 4), rep(mycol[1], 14) )
  
  pdf( paste0(x, '_boxplot.pdf') )
  boxplot( varsels_mean, ylim = c(0, 1), col = use_col, outcol = use_col, 
           pch = 20, main = paste0( graph_title ),
           ylab = 'Selection Probability', las = 2 )
  abline( h = seq(0.2, 1, 0.2), lty = 3, col = 'lightgray')
  dev.off()
  
  NULL
  
}

Map( boxplot_make, x = filenames, graph_title = allnms )


### Beta distribution plot
xx <- seq( 0.001, 0.999, 0.001 )
yy1 <- dbeta( xx, 1, 19 )
yy2 <- dbeta( xx, 3, 17 )
yy3 <- dbeta( xx, 1, 1 )
yy4 <- dbeta( xx, 3, 3 )

pdf( 'Beta_Hyperprior.pdf' )
plot( xx, yy1 / sum(yy1), type = 'l', lwd = 2, col = mycol[1], ylab = 'Beta Distribution Weight', xlab = '',  )
lines( xx, yy3 / sum(yy3), lwd = 2, col = mycol[1], lty = 3 )
lines( xx, yy2 / sum(yy2), lwd = 2, col = mycol[2] )
lines( xx, yy4 / sum(yy4), lwd = 2, col = mycol[2], lty = 3 )
grid()
text( 0.1, 0.015, 'Beta(1, 19)', col = mycol[1] )
text( 0.9, 0.0015, 'Beta(1, 1)', col = mycol[1] )
text( 0.25, 0.005, 'Beta(3, 17)', col = mycol[2] )
text( 0.5, 0.0023, 'Beta(3, 3)', col = mycol[2] )
dev.off()
