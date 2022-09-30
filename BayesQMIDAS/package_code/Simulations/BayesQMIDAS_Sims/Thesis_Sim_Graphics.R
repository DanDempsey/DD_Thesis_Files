##### Code for creating Bayes-Q-MIDAS simulation graphics
##### Daniel Dempsey

library( ggplot2 )
library( dplyr )
library( xtable )

mycol <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",  
           "#0072B2", "#D55E00", "#CC79A7", "#999999", '#FFFFFF')
scales::show_col(mycol)

### Load all data
setwd( 'package_code/Simulations/BayesQMIDAS_Sims/' )
all_files <- list.files( pattern = '*Rdata' )

retained <- seq( 50001, 100000, 10 )

### Variable Selection Boxplots
varsel_plot <- function( ind, keep = retained, modlab = 1, colour = 'dodgerblue', fn = NULL, ... ) {
  
  load( all_files[ ind ] )
  
  varsels_list <- lapply( res$MCMC, function( x ) { x$vsres[keep, ] } )
  varsels_mean_list <- lapply( varsels_list, function(x) { apply(x, 2, mean) } )
  varsels_mean <- do.call( 'rbind', varsels_mean_list )
  quantlab <- res$MCMC[[1]]$quantile
  
  if (!is.null(fn)) { pdf( paste0( 'Graphs/Variable_Selection_Boxplots/', fn ) ) }
  boxplot( varsels_mean, ylim = c(0, 1), col = colour, outcol = colour, 
           main = paste0( 'Model ', modlab, ', Quantile = ', quantlab ),
           ylab = 'Proportion of Variable Selected', xaxt = 'n', ... )
  if (!is.null(fn)) { dev.off() }
  
}

col_scheme2 <- c( rep( 'dodgerblue', 5 ), rep( 'orange', 15 ) )

varsel_plot( 1, fn = 'Box1_5.pdf' )
varsel_plot( 2, fn = 'Box1_2.pdf' )
varsel_plot( 3, colour = col_scheme2, modlab = 2, fn = 'Box2_5.pdf' )
varsel_plot( 4, colour = col_scheme2, modlab = 2, fn = 'Box2_2.pdf' )

### PPD AUC Histograms
ppd_plot <- function( ind, modlab = 1, fn = NULL, ... ) {
  
  load( all_files[ ind ] )
  
  ppd_test <- sapply( res$ppd_test, function(x) { x$ROC$auc[1] } )
  ppd_train <- sapply( res$ppd_train, function(x) { x$ROC$auc[1] } )
  quantlab <- res$MCMC[[1]]$quantile
  
  ppd_dat <- data.frame( Set = c( rep( "Train", length( ppd_train ) ), rep( "Test", length( ppd_test ) ) ),
                         Data = c( ppd_train, ppd_test ) )
  
  if (!is.null(fn)) { pdf(paste0( 'Graphs/PPD_AUC_Histograms/Train_', fn)) }
  hist( ppd_train, col = 'yellow', xlab = '', xlim = c( 0, 1 ),
        main = paste0( 'Model ', modlab, ', Quantile = ', quantlab ), ... )
  if (!is.null(fn)) { dev.off() }
  
  if (!is.null(fn)) { pdf(paste0( 'Graphs/PPD_AUC_Histograms/Test_', fn)) }
  hist( ppd_test, col = 'yellow', xlab = '', xlim = c( 0, 1 ),
        main = paste0( 'Model ', modlab, ', Quantile = ', quantlab ), ... )
  if (!is.null(fn)) { dev.off() }
  
  ppd_dat
  
}

dat1 <- ppd_plot( 1, fn = 'Hist1_5.pdf' )
dat2 <- ppd_plot( 2, fn = 'Hist1_2.pdf' )
dat3 <- ppd_plot( 3, modlab = 2, fn = 'Hist2_5.pdf' )
dat4 <- ppd_plot( 4, modlab = 2, fn = 'Hist2_2.pdf' )

all_dat <- list( dat1, dat2, dat3, dat4 )

ppd_agg <- function( dat, type, fun ) {
  
  x <- filter( dat, Set == type )
  fun( x$Data )
  
}

modlabs <- rep( 1:2, each = 2 )
quants_dat <- rep( quants, 2 )

train_dat <- data.frame( modlabs, quants_dat, do.call( 'rbind', lapply( all_dat, ppd_agg, type = 'Train', fun = range ) ) )
colnames( train_dat ) <- c( 'Model', 'Quantile', 'Lower', 'Upper' )

test_dat <- data.frame( modlabs, quants_dat, do.call( 'rbind', lapply( all_dat, ppd_agg, type = 'Test', fun = range ) ) )
colnames( test_dat ) <- colnames( train_dat )

cap <- paste0( 'Range of AUC values for the ', c( 'training', 'test' ), ' simulation data.' )
tablabs <- paste0( 'tab:', c( 'train', 'test' ), '_AUC' )

xtable( train_dat, caption = cap[ 1 ], label = tablabs[ 1 ] )
xtable( test_dat, caption = cap[ 2 ], label = tablabs[ 2 ] )

### Beta Traceplots
varsel_change <- function( x, y ) {
  x[which( !y )] <- NA
  x
}

hline <- function( x, col = 'dodgerblue', ... ) {
  abline( h = x, lty = 2, col = col, ... )
}

beta_plot <- function( res, i, keep = retained, fn = NULL, ... ) {
  
  true <- res$true_beta[ i + 1 ]
  modlab <- ifelse( length( res$true_beta ) < 4, 1, 2 )
  quantlab <- res$q * 10
  
  betas_list <- lapply( res$MCMC, function(x) { x$betares[keep, i+1] } )
  if ( i == 0 ) { vsres_list <- lapply( seq_along( res$MCMC ), function( x ) { rep( TRUE, length( keep ) ) } ) }
  else {
    vsres_list <- lapply( res$MCMC, function(x) { x$vsres[keep, i] } )
  }
  
  beta_list_NA <- Map( varsel_change, x = betas_list, y = vsres_list )
  betas <- do.call( 'cbind', beta_list_NA )
  
  if ( !is.null(fn) ) { pdf( paste0( 'Graphs/Beta_Traceplots/Fit', modlab, '_', quantlab, '/', fn, '_c', i, '.pdf' ) ) }
  plot( 0, type = 'n', ylim = range( betas, na.rm = TRUE ), xlim = c(0, length( keep ) ),
        ylab = '', main = paste0( 'Covariate ', i, ' beta' ), ... )
  grid()
  for ( j in 1:ncol(betas) ) {
    lines( betas[, j], col = 'lightgray' )
  }
  
  beta_mean <- rowSums( betas, na.rm = TRUE ) / rowSums( do.call('cbind', vsres_list) )
  lines( beta_mean, col = 'dodgerblue' )
  hline( true, col = 'orange', lwd = 2 )
  if ( !is.null(fn) ) { dev.off() }
  
}

beta_plot_wrapper <- function( ind, maxval = 6, keep = retained, fn = NULL, ... ) {
  
  load( all_files[ ind ] )
  nvars <- 1:min( ncol( res$MCMC[[1]]$betares ), maxval ) - 1
  
  for ( i in nvars ) {
    beta_plot( res = res, i = i, keep = keep, fn = fn, ... )
  }
  
}

lapply( 1:4, beta_plot_wrapper, keep = retained, fn = 'BetaTrace' )

### Theta Traceplots
DLF_plot <- function( res, i, keep = retained, fn = NULL, ... ) {
  
  DLF_list <- lapply( res$MCMC, function(x) { x$DLFres[[i]] } )
  DLF1_all <- as.data.frame( do.call( 'cbind', lapply(DLF_list, '[', i = keep, j = 1) ) )
  DLF2_all <- as.data.frame( do.call( 'cbind', lapply(DLF_list, '[', i = keep, j = 2) ) )
  
  vsres_list <- lapply( res$MCMC, function(x) { x$vsres[keep, i] } )
  
  DLF1 <- as.data.frame( Map( varsel_change, x = DLF1_all, y = vsres_list ) )
  DLF2 <- as.data.frame( Map( varsel_change, x = DLF2_all, y = vsres_list ) )
  
  all_DLF <- list( DLF1, DLF2 )
  true <- res$true_DLF[[ i ]]
  modlab <- ifelse( length( res$true_beta ) < 4, 1, 2 )
  quantlab <- res$q * 10
  
  if (!is.null(fn)) { pdf(paste0('Graphs/Theta_Traceplots/Fit', modlab, '_', quantlab, '/', fn, '_c', i, '.pdf')) }
  par( mfrow = c(2, 1) )
  for ( k in 1:2 ) {
    DLF <- all_DLF[[k]]
    plot( DLF[, 1], type = 'n', ylim = range(DLF, na.rm = TRUE), main = '', xlab = '',
          ylab = paste0('Covariate ', i, ' Theta ', k), ... )
    grid()
    for ( j in 1:ncol(DLF) ) {
      lines( DLF[, j], col = 'lightgray' )
    }
    DLF_mean <- rowSums( DLF, na.rm = TRUE ) / rowSums( do.call('cbind', vsres_list) )
    
    lines( DLF_mean, col = 'dodgerblue' )
    hline( true[k], col = 'orange', lwd = 2 )
  }
  if (!is.null(fn)) { dev.off() }
  
}

DLF_plot_wrapper <- function( ind, maxval = 6, keep = retained, fn = NULL, ... ) {
  
  load( all_files[ ind ] )
  nvars <- 2:min( ncol( res$MCMC[[1]]$betares ), maxval ) - 1
  
  for ( i in nvars ) {
    DLF_plot( res = res, i = i, keep = keep, fn = fn, ... )
  }
  
}

lapply( 1:4, DLF_plot_wrapper, keep = retained, fn = 'DLFTrace' )

