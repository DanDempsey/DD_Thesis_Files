##### Easy Sim Code
##### Daniel Dempsey

source('package_code/1_load_all.R')
library( readr )
library( dplyr )
#setwd( 'package_code/Simulations/Big_Runs/' )

### Colour pallette for plots, taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",  
           "#0072B2", "#D55E00", "#CC79A7", "#999999", '#FFFFFF')
scales::show_col(mycol)

### Read in SKT data, scale it and add white noise
skt_dat <- read_csv( 'Flare_Modelling_Application/Data/Processed_Data/Final_Environmet_Flare_Data.csv' ) %>%
  select( DATE, RKD_ID, SKT ) %>% filter( RKD_ID == 103 ) %>% select( -RKD_ID )
skt_dat$SKT_scale <- as.numeric( scale( skt_dat$SKT ) )
plot( SKT_scale ~ DATE, data = skt_dat, type = 'l' )
date_len <- length( skt_dat$DATE )

effect_size <- 1.5
true_DLF <- c(7, -1)

set.seed( 100 )
skt_fit <- vector( mode = "list", length = 15 )
for ( i in 1:15 ) {
  skt_dat$noise <- rnorm( nrow(skt_dat) )
  #plot( noise ~ DATE, data = skt_dat, type = 'l' )

  ### Simulate flare data based on SKT
  date_inds <- date_len + 1
  while ( any(date_inds > date_len) ) {
    date_inds <- round( cumsum( rexp(150, 1/5) ) ) + 30
  }
  flare_dates <- skt_dat$DATE[ date_inds ]
  
  skt_regime1 <- midas_regimes( value = skt_dat[3:4], time = skt_dat$DATE, 
                                time_window = 30, DLF_parameters = true_DLF )
  
  # Simulate the response index and find the time windows
  skt_regime2 <- lapply( skt_regime1, "[[<-", i = "response_index", value = flare_dates )
  skt_regime <- lapply( skt_regime2, time_delta_list )
  
  # Compute the midas design matrix and thus retrieve the linear predictor
  WX <- cbind( 1, midas_design_matrices(skt_regime) )
  lin_pred <- WX %*% c(0, effect_size, 0)
  
  # Now simulate the y values based on lin_pred
  y <- ifelse( rALD(n = length( lin_pred ), mu = lin_pred, p = 0.5) >= 0, 1, 0 )
  
  #plot( SKT_scale ~ DATE, data = skt_dat, type = 'l' )
  #points( flare_dates, y, pch = 20, col = 'blue' )
  
  ### Fit the model
  mrtest <- midas_regimes( value = skt_dat[, 3:4], time = skt_dat$DATE, 
                           DLF_parameters = list(c(0, -1)), time_window = 30 )
  rtest <- reponse_regime( value = y, time = flare_dates )
  form <- rtest ~ mrtest
  
  skt_fit[[i]] <- IRTS_MIDAS_AuxVar( formula = form, data = skt_dat, quantile = 0.5,
                                     response_dist = "ald", varsel = TRUE, MCMC_length = 10000 )
}
 
save( skt_fit, file = 'package_code/Simulations/Single_Dataset_Runs/skt_res.Rdata' )
#load( file = 'package_code/Simulations/Single_Dataset_Runs/skt_res.Rdata' ) 

retained <- seq( 5000, 10000, 10 )

### Graphing functions
varsel_plot <- function( res, keep = retained, colour = 'dodgerblue', fn = NULL, ... ) {
  
  varsels_list <- lapply( res, function(x) { x$vsres[keep, ] } )
  varsels_mean_list <- lapply( varsels_list, function(x) { apply(x, 2, mean) } )
  varsels_mean <- do.call( 'rbind', varsels_mean_list )
  
  if (!is.null(fn)) { pdf(fn) }
  boxplot( varsels_mean, ylim = c(0, 1), col = colour, outcol = colour, ... )
  if (!is.null(fn)) { dev.off() }
  
}

varsel_change <- function( x, y ) {
  x[which( !y )] <- NA
  x
}

hline <- function( x, col = 'red', ... ) {
  abline( h = x, lty = 2, col = col, ... )
}

beta_plot <- function( res, i, true, keep = retained, fn = NULL, ... ) {
  
  betas_list <- lapply( res, function(x) { x$betares[keep, i+1] } )
  if ( i == 0 ) { vsres_list <- TRUE }
  else {
    vsres_list <- lapply( res, function(x) { x$vsres[keep, i] } )
  }
  
  beta_list_NA <- Map( varsel_change, x = betas_list, y = vsres_list )
  betas <- do.call( 'cbind', beta_list_NA )
  
  if ( !is.null(fn) ) { pdf(fn) }
  plot( betas[, 1], type = 'n', ylim = range(betas, na.rm = TRUE), ylab = '',
        main = paste0('Covariate ', i, ' beta'), ... )
  grid()
  for ( j in 1:ncol(betas) ) {
    lines( betas[, j], col = 'lightgray' )
  }
  
  if ( i == 0 ) {
    beta_mean <- rowSums( betas, na.rm = TRUE ) / 15
  }
  else {
    beta_mean <- rowSums( betas, na.rm = TRUE ) / rowSums( do.call('cbind', vsres_list) )
  }
  
  lines( beta_mean, col = 'dodgerblue' )
  hline( true, col = 'orange', lwd = 2 )
  if ( !is.null(fn) ) { dev.off() }
  
}

DLF_plot <- function( res, i, true, keep = retained, fn = NULL, ... ) {
  
  DLF_list <- lapply( res, function(x) { x$DLFres[[i]] } )
  DLF1_all <- as.data.frame( do.call( 'cbind', lapply(DLF_list, '[', i = keep, j = 1) ) )
  DLF2_all <- as.data.frame( do.call( 'cbind', lapply(DLF_list, '[', i = keep, j = 2) ) )
  
  vsres_list <- lapply( res, function(x) { x$vsres[keep, i] } )
  
  DLF1 <- as.data.frame( Map( varsel_change, x = DLF1_all, y = vsres_list ) )
  DLF2 <- as.data.frame( Map( varsel_change, x = DLF2_all, y = vsres_list ) )
  
  all_DLF <- list( DLF1, DLF2 )
  
  for ( k in 1:2 ) {
    if (!is.null(fn)) { pdf(paste0('DLF', k, '/', fn)) }
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
    if (!is.null(fn)) { dev.off() }
  }
  
}

### Make plots
setwd('package_code/Simulations/Graphs/SKT_Sim_Graphs/')

pdf( 'Variable_Selection_Plot.pdf', width = 10 )
varsel_plot( skt_fit, col = mycol[c(2, 1)] )
dev.off()

pdf( 'Beta_int.pdf' )
beta_plot( skt_fit, 0, 0 )
dev.off()

pdf( 'Beta_SKT.pdf' )
beta_plot( skt_fit, 1, effect_size )
dev.off()

pdf( 'Beta_junk.pdf' )
beta_plot( skt_fit, 2, 0 )
dev.off()

pdf( 'DLF_pars.pdf' )
par( mfrow = c(2, 1) )
DLF_plot( skt_fit, 1, true_DLF )
par( mfrow = c(1, 1) )
dev.off()

