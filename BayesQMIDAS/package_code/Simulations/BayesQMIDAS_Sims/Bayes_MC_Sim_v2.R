##### Big Sim Code
##### Daniel Dempsey

source( 'package_code/1_load_all.R' )
setwd( 'package_code/Simulations/BayesQMIDAS_Sims/Output' )
library( parallel )
library( corrplot )
library( ggplot2 ) # For plotting ridgeplots
library( ggridges ) # Ditto ^
library( tictoc )

### Colour pallette for plots, taken from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
mycol <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",  
           "#0072B2", "#D55E00", "#CC79A7", "#999999", '#FFFFFF')
scales::show_col(mycol)

### Utility functions
ALD_innov <- function(n, x, q = 0.5) {
  ifelse( rALD(n, mu = x, p = q) >= 0, 1, 0 )
}

varsel_change <- function( x, y ) {
  x[which( !y )] <- NA
  x
}

beta_extract <- function( x, k, i ) {
  x$betares[k, i]
}

theta_extract <- function( x, k, i, j ) {
  x$DLFres[[j]][k, i]
}

gamma_extract <- function( x, k, i ) {
  x$vsres[k, i]
}

run_sim <- function(b = c(0, -1.5, 1.5), DLF_fun = 'irts_nealmon', DLF = list(c(14, -1), c(0.5, -0.1)), q = 0.5, ny = 1000,
                    s = c(7, 2.5), p = 30, MCMC_len = 50000, burn_in = 25000, thin = 10, tw_sim = 15, 
                    tw_mod = 15, use_prior = NULL, rr = 1/15, file_path = 'Test', sb = 1, DLF_mod_pars = list(c(0, -1)),
                    reps = 1, n_cores = 1) {

  # Simulate data
  n_vars <- length(b) - 1
  #if ( missing(which_adjust) ) { which_adjust <- 1:n_vars }
  sim_dat <- irts_midas_sim(n_y = ny, innov_fun = ALD_innov, shift = s, periodicity = p,
                            n_vars = n_vars, beta = b, DLF_parameters = DLF, DLF = DLF_fun,
                            seasonal_adjust = TRUE, innov_pars = list(q = q),
                            time_window = tw_sim, response_rate = rr,
                            sigma_base = diag(sb, n_vars))$Data
  
  # Split into training/test sets
  y_full <- na.omit( sim_dat[, c('INDEX', 'y')] )
  
  INDEX_train <- y_full[ 1:(0.8*ny), 'INDEX' ]
  sim_dat_train <- sim_dat[ which( sim_dat$INDEX <= max( INDEX_train ) ), ]
  
  INDEX_test <- setdiff( y_full$INDEX, INDEX_train )
  sim_dat_test <- sim_dat[ which( sim_dat$INDEX > ( min( INDEX_test ) - tw_mod ) ), ]
  drop_inds <- which( !is.na( sim_dat_test$y ) & ( sim_dat_test$INDEX < min( INDEX_test ) )  )
  if( length(drop_inds) > 0 ) { sim_dat_test <- sim_dat_test[ -drop_inds, ] }
  
  imba <- sum( y_full$y ) / length( y_full$y )
  rm( y_full, INDEX_train, INDEX_test, drop_inds )
  
  # Fit the model
  mrtest <- midas_regimes( value = sim_dat_train[, seq(1, n_vars)], time = sim_dat_train$INDEX, 
                           DLF_parameters = DLF_mod_pars, time_window = tw_mod )
  rtest <- response_regime( value = sim_dat_train$y, time = sim_dat_train$INDEX )
  form <- rtest ~ mrtest
  
  cor_dat = cor( na.omit( sim_dat_train[, 1:n_vars] ) )
  
  corrplot( cor_dat, method = 'number' )
  
  run_fun <- function( xx ) {
    IRTS_MIDAS_AuxVar( formula = form, data = sim_dat_train, quantile = q,
                       prior = use_prior, response_dist = "ald", varsel = TRUE,
                       MCMC_length = MCMC_len )[c(1:5, 8)]
  }
  res <- mclapply( 1:reps, run_fun, mc.cores = n_cores )
  
  #browser()
  
  # Training ppd
  ppd_res_train <- mclapply( res, ppd_ALD_vs, new_form = form, 
                             burn_in = burn_in, thin = thin, 
                             mc.cores = n_cores )
  
  ppd_res_train <- sapply( lapply(ppd_res_train, '[[', i = 'ROC'), '[[', i = 'auc' )
  
  # Test ppd
  mrtest_test <- midas_regimes( value = sim_dat_test[, seq(1, n_vars)], 
                                time = sim_dat_test$INDEX, 
                                DLF_parameters = list(c(0, -1)), 
                                time_window = tw_mod )
  rtest_test <- response_regime( value = sim_dat_test$y, time = sim_dat_test$INDEX )
  form_test <- rtest_test ~ mrtest_test
  
  ppd_res_test <- mclapply( res, ppd_ALD_vs, new_form = form_test, 
                            burn_in = burn_in, thin = thin, 
                            mc.cores = n_cores )
  
  ppd_res_test <- sapply( lapply(ppd_res_test, '[[', i = 'ROC'), '[[', i = 'auc' )
  
  keep <- seq( burn_in, MCMC_len, thin )
  if( !dir.exists(file_path) ) { dir.create(file_path) } 
  setwd( file_path )
  
  ### Selection plots
  varsels_list <- lapply( res, function( x ) { x$vsres[keep, ] } )
  varsels_mean_list <- lapply( varsels_list, function(x) { apply(x, 2, mean) } )
  varsels_mean <- do.call( 'rbind', varsels_mean_list )
  quantlab <- res[[1]]$quantile
  
  use_col <- ifelse( b[-1] == 0, mycol[1], mycol[2] )
  
  pdf( 'Variable_Selection_Boxplot.pdf' )
  boxplot( varsels_mean, ylim = c(0, 1), col = use_col, outcol = use_col, 
           pch = 20, main = paste0( 'Quantile = ', quantlab ),
           ylab = 'Selection Probability', las = 2 )
  abline( h = seq(0.2, 1, 0.2), lty = 3, col = 'lightgray')
  dev.off()
  
  ### Ridgeplots
  len <- length( keep )
  for ( j in 1:(n_vars+1) ) {
    
    beta_vec <- do.call( 'c', lapply( res, beta_extract, k = keep, i = j ) )
    if ( j == 1 ) {
      gamma_vec <- rep( FALSE, len )
    }
    else {
      gamma_vec <- !(do.call( 'c', lapply( res, gamma_extract, k = keep, i = j-1 ) ))
    }
    
    beta_vec[gamma_vec] <- NA
    if(all(is.na(beta_vec))) { next }
    graph_dat <- data.frame( beta = beta_vec, ind = factor(rep(1:reps, each = len)) )
    
    rplot <- ggplot(graph_dat, aes(x = beta, y = ind)) +
      geom_density_ridges( fill = mycol[2] ) +
      theme_ridges() + xlab('') + ylab('') +
      geom_vline(xintercept = b[j], colour = mycol[1]) +
      ggtitle( paste0('Variable ', j - 1) ) +
      theme(legend.position = "none")
    
    pdf( paste0('Beta_Ridgeplot_', j-1, '.pdf') )
    print( rplot )
    dev.off()
    
  }
  
  ### Theta Traceplots
  for ( j in 1:n_vars ) {
    
    for ( k in 1:2 ) {
      
      theta_vec <- do.call( 'c', lapply( res, theta_extract, k = keep, i = k, j = j ) )
      gamma_vec <- !(do.call( 'c', lapply( res, gamma_extract, k = keep, i = j ) ))
      
      theta_vec[gamma_vec] <- NA
      if(all(is.na(theta_vec))) { next }
      graph_dat <- data.frame( theta = theta_vec, ind = factor(rep(1:reps, each = len)) )
      
      rplot <- ggplot(graph_dat, aes(x = theta, y = ind)) +
        geom_density_ridges( fill = mycol[2] ) +
        theme_ridges() + xlab('') + ylab('') +
        geom_vline(xintercept = DLF[[j]][k], colour = mycol[1]) +
        ggtitle( paste0('Variable ', j, ', Theta ', k) ) +
        theme(legend.position = "none")
      
      pdf( paste0('Theta_Ridgeplot_', j, '_', k, '.pdf') )
      print( rplot )
      dev.off()
      
    }
    
  }
  
  keep_beta <- lapply( res, function(x) { x$betares[, -1] } )
  DLF_attempts <- lapply( keep_beta, function(x) { colSums(x != 0) } )
  DLF_accepts <- lapply( res, function(x){ x$DLF_accept } )
  DLF_accept_prop <- Map( '/', x = DLF_accepts, y = DLF_attempts )
  
  # Compile results
  all_res <- list( MCMC = res, ppd_train = ppd_res_train, ppd_test = ppd_res_test, 
                   DLF_accept_prop = DLF_accept_prop, data_correlation = cor_dat,
                   #train_dat = sim_dat_train, test_dat = sim_dat_test, 
                   imbalance = imba, true_beta = b, true_DLF = DLF, q = q )
  save( all_res, file = "res.Rdata" )
  
  setwd( '..' )
  
}

run_sim_loop_wrapper <- function( x, ... ) {
  run_sim( ... )
}

### Set simulation parameters
b_test1 <- c(0, 1.5, -1.5)
b_test2 <- c(0, 1.5, -1.5, 1, -1, rep(0, 14))
b_test3 <- c(0, 1.5, -1.5, rep(0, 16))

blen1 <- length( b_test1 )
blen2 <- length( b_test2 )

myprior1 <- list( beta0 = rep(0, blen1), V0 = 10 * diag(blen1), 
                  DLF_pars = list( DLF1 = c(15, 5), DLF2 = c(1, 1) ),
                  model_selection = c(3, 3) )
myprior2 <- list( beta0 = rep(0, blen2), V0 = 10 * diag(blen2), 
                  DLF_pars = list( DLF1 = c(15, 5), DLF2 = c(1, 1) ),
                  model_selection = c(3, 3) )
myprior2_11 <- list( beta0 = rep(0, blen2), V0 = 10 * diag(blen2), 
                  DLF_pars = list( DLF1 = c(15, 5), DLF2 = c(1, 1) ),
                  model_selection = c(1, 1) )


DLF1 <- list( c(14, -1), c(0.5, -0.1) )
DLF2 <- rep( list( c(14, -1), c(0.5, -0.1) ), (blen2 - 1)/2 )

#reps <- 6
reps <- 30
n_cores <- 10

### Run simulations
#set.seed( 111 )
#tic()
#run_sim( b = b_test1, use_prior = myprior1, DLF = DLF1, 
#         q = 0.5, n_cores = n_cores, reps = reps,
#         MCMC_len = 20000, burn_in = 10000,
#         #MCMC_len = 1000, burn_in = 500,
#         file_path = 'Small_5' )
#toc() # 4518.627 sec elapsed
#
#set.seed( 222 )
#tic()
#run_sim( b = b_test1, use_prior = myprior1, DLF = DLF1, 
#         q = 0.2, n_cores = n_cores, reps = reps,
#         MCMC_len = 20000, burn_in = 10000, 
#         file_path = 'Small_2' )
#toc() # 4484.367 sec elapsed

set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2_11, DLF = DLF2, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large2_5_multilevel_11' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2_11, DLF = DLF2,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_multilevel_11' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2_11, DLF = DLF2, 
         #p = c( rep(30, 2), rep(15, 2), rep(5, 14) ),
         q = 0.5, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large_5_multilevel_11' )
toc() # 40572.806 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2_11, DLF = DLF2, 
         q = 0.2, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_multilevel_11' )
toc() # 57157.414 sec elapsed



set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large2_5_multilevel_33' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_multilevel_33' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         #p = c( rep(30, 2), rep(15, 2), rep(5, 14) ),
         q = 0.5, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large_5_multilevel_33' )
toc() # 40572.806 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.2, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_multilevel_33' )
toc() # 57157.414 sec elapsed



### Try different sigmas
set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps, sb = 3,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large2_5_sb3' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps, sb = 3,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_sb3' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         #p = c( rep(30, 2), rep(15, 2), rep(5, 14) ),
         q = 0.5, n_cores = n_cores, reps = reps, sb = 3,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 10000, burn_in = 5000,
         file_path = 'Large_5_sb3' )
toc() # 52599.612 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.2, n_cores = n_cores, reps = reps, sb = 3,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_sb3' )
toc() # 60383.523 sec elapsed




### Try different simulation DLF
a1 <- 30
b1 <- 12

a2 <- 2
b2 <- 4

DLF1_beta <- list( c(a1, b1), c(a2, b2) )
DLF2_beta <- rep( list( c(a1, b1), c(a2, b2) ), (blen2 - 1)/2 )

set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2_beta, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps, DLF_fun = 'irts_beta',
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_5_beta' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2_beta,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps, DLF_fun = 'irts_beta',
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_beta' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2_beta, 
         q = 0.5, n_cores = n_cores, reps = reps, DLF_fun = 'irts_beta',
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_5_beta' )
toc() # 68412.977 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2_beta, 
         q = 0.2, n_cores = n_cores, reps = reps, DLF_fun = 'irts_beta',
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_beta' )
toc() # 54391.596 sec elapsed


### Try different time window
set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps, tw_mod = 25,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_5_tw25' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps, tw_mod = 25,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_tw25' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.5, n_cores = n_cores, reps = reps, tw_mod = 25,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_5_tw25' )
toc() # 68412.977 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.2, n_cores = n_cores, reps = reps, tw_mod = 25,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_tw25' )
toc() # 54391.596 sec elapsed

set.seed( 111 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.5, n_cores = n_cores, reps = reps, tw_mod = 7,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_5_tw7' )
toc() # 40572.806 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2,
         s = sample(1:20, 18, replace = TRUE),
         p = sample(15:60, 18, replace = TRUE),
         q = 0.2, n_cores = n_cores, reps = reps, tw_mod = 7,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large2_2_tw7' )
toc() # 57157.414 sec elapsed

set.seed( 333 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.5, n_cores = n_cores, reps = reps, tw_mod = 7,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_5_tw7' )
toc() # 68412.977 sec elapsed

set.seed( 444 )
tic()
run_sim( b = b_test2, use_prior = myprior2, DLF = DLF2, 
         q = 0.2, n_cores = n_cores, reps = reps, tw_mod = 7,
         MCMC_len = 80000, burn_in = 40000, 
         file_path = 'Large_2_tw7' )
toc() # 54391.596 sec elapsed







### Extra test (easy)
set.seed( 111 )
tic()
run_sim( b = b_test3, use_prior = myprior2, DLF = DLF2, 
         s = c(7, 2.5, sample(1:10, 16, replace = TRUE)),
         p = c(rep(30, 2), sample(10:60, 16, replace = TRUE)), 
         q = 0.5, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000,
         #MCMC_len = 1000, burn_in = 500,
         file_path = 'Large_5_LowMulitcol' )
toc() # 68412.977 sec elapsed

set.seed( 222 )
tic()
run_sim( b = b_test3, use_prior = myprior2, DLF = DLF2, 
         s = c(7, 2.5, sample(1:10, 16, replace = TRUE)),
         p = c(rep(30, 2), sample(10:60, 16, replace = TRUE)),
         q = 0.2, n_cores = n_cores, reps = reps,
         MCMC_len = 80000, burn_in = 40000, 
         #MCMC_len = 1000, burn_in = 500,
         file_path = 'Large_2_LowMulitcol' )
toc() # 54391.596 sec elapsed
