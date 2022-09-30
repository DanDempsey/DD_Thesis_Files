### Run Simulation
library(SparseM)
library(corrplot)
library(tictoc)
library(matrixcalc)

source('package_code/1_load_all.R')

ALD_innov <- function(n, x, q = 0.5) {
  ifelse( rALD(n, mu = x, p = q) >= 0, 1, 0 )
}

vs_mean <- function(x, k) {
  mean( x[k] )
}

calculate_time <- function(x) {
  round( x$toc - x$tic )
}

hline <- function( x, col = 'red', ... ) {
  abline( h = x, lty = 2, col = col, ... )
}

vline <- function( x, col = 'red', ... ) {
  abline( v = x, lty = 2, col = col, ... )
}

run_sim <- function(b = c(0, -2, 2), DLF = list(c(14, -1)), q = 0.5, ny = 1000, nx = 5000, seed,
                    MCMC_len = 5000L, burn_in = 1000, thin = 10, time_window = 50, 
                    DLF_fixed = FALSE, which_adjust, seasonal_adjust = FALSE, ...) {
  
  #browser()
  if (!missing(seed)) { set.seed(seed) }
  
  # Simulate data
  n_vars <- length(b) - 1
  if ( missing(which_adjust) ) { which_adjust <- 1:n_vars }
  sim <- irts_midas_sim(n_y = ny, n_x_full = nx, innov_fun = ALD_innov, n_num = n_vars, beta = b, DLF_parameters = DLF,
                        innov_pars = list(q = q), seasonal_adjust = seasonal_adjust, which_adjust = which_adjust)
  sim_dat <- sim$Data
  
  # Split into training/test sets
  y_full <- na.omit( sim_dat[, c('INDEX', 'y')] )
  N <- nrow( y_full )
  
  INDEX_train <- y_full[ 1:(0.8*N), 'INDEX' ]
  sim_dat_train <- sim_dat[ 1:which( sim_dat$INDEX == max(INDEX_train) ), ]
  
  INDEX_test <- y_full[ ((0.8*N) + 1):N, ]
  sim_dat_test <- sim_dat[ which(sim_dat$INDEX >= min(INDEX_test$INDEX) - time_window), ]

  imba <- sum(y_full$y) / length(y_full$y)
  rm(y_full, N, INDEX_train, INDEX_test)
  
  # Fit the model
  mrtest <- midas_regimes(value = sim_dat_train[, seq(1, n_vars)], time = sim_dat_train$INDEX, 
                          DLF_parameters = list(c(0, -1)), time_window = time_window)
  rtest <- reponse_regime(value = sim_dat_train$y, time = sim_dat_train$INDEX)
  form <- rtest ~ mrtest
  
  cat( 'Data simulated. Now running MCMC inference.\n' )
  
  corrplot( cor(na.omit( sim_dat_train[, 1:n_vars])), method = 'number' )
  
  tic()
  res <- IRTS_MIDAS_AuxVar(form, data = sim_dat_train, quantile = q,
                           response_dist = "ald", varsel = TRUE,
                           nburn = burn_in, DLF_fixed = DLF_fixed,
                           MCMC_length = MCMC_len, ...)
  timer_MCMC <- calculate_time( toc(quiet = TRUE) )
  
  if (DLF_fixed) {
    
    all_res <- list()
    all_res$MCMC <- res
    
    all_res$train_dat <- sim_dat_train
    all_res$test_dat <- sim_dat_test
    all_res$imbalance <- imba
    
    return( all_res )
    
  }
  
  # Training ppd
  cat( paste0('Complete. Time taken: ', timer_MCMC, ' seconds.\nNow computing PPD for training data.\n') )
  tic()
  ppd_res_train <- ppd_ALD_vs( MCMC_res = res, new_form = form, 
                               burn_in = burn_in, thin = thin )
  timer_ppd_train <- calculate_time( toc(quiet = TRUE) )
  cat( paste0('Complete. Time taken: ', timer_ppd_train, ' seconds.\nNow computing PPD for test data.\n') )
  
  # Test ppd
  mrtest_test <- midas_regimes(value = sim_dat_test[, seq(1, n_vars)], 
                               time = sim_dat_test$INDEX, 
                               DLF_parameters = list(c(0, -1)), 
                               time_window = 10)
  rtest_test <- reponse_regime(value = sim_dat_test$y, time = sim_dat_test$INDEX)
  form_test <- rtest_test ~ mrtest_test
  
  tic()
  ppd_res_test <- ppd_ALD_vs( MCMC_res = res, new_form = form_test, 
                              burn_in = burn_in, thin = thin )
  timer_ppd_test <- calculate_time( toc(quiet = TRUE) )
  cat( paste0('Complete. Time taken: ', timer_ppd_test, ' seconds.\n') )
  
  # Compile results
  all_res <- list()
  all_res$MCMC <- res
  
  all_res$ppd_train <- ppd_res_train
  all_res$ppd_test <- ppd_res_test
  
  ttt <- sum( timer_MCMC, timer_ppd_train, timer_ppd_test )
  all_res$timer_MCMC <- timer_MCMC
  all_res$timer_ppd_train <- timer_ppd_train
  all_res$timer_ppd_test <- timer_ppd_test
  all_res$total_time_taken <- ttt
  
  all_res$train_dat <- sim_dat_train
  all_res$test_dat <- sim_dat_test
  all_res$imbalance <- imba
  
  cat( paste0('Simulation study completed. Total time taken: ', round( ttt / 3600, 2 ), ' hours.\n') ) 
  all_res
  
}

if (FALSE) {
  
myprior <- list( beta0 = rep(0, 3), V0 = 100 * diag(3), 
                 DLF_pars = list( DLF1 = c(14, 7), DLF2 = c(1, 1) ),
                 model_selection = c(1, 1) )

#res <- run_sim(b = c(0, 0.8, -0.8, rep(0, 7)), q = 0.5, MCMC_len = 100000, 
#               seed = 15, burn_in = 50000, thin = 10, prior = myprior)
res <- run_sim(b = c(0, 0.8, -0.8), DLF = c(0, -0.1), q = 0.5, MCMC_len = 10000, 
               seed = 10, burn_in = 5000, thin = 10, prior = myprior)
#save(res, file = 'package_code/Simulations/seed10_q2_myprior.Rdata')
#load(file = 'package_code/Simulations/seed10_q5.Rdata')

### Accepted draws
#keep <- seq( 50000, 100000, 10 )
keep <- seq( 5000, 10000, 10 )
nvars <- ncol( res$MCMC$betares ) - 1

### AUC
res$ppd_train
res$ppd_test

### Variable selection 
barplot( apply( res$MCMC$vsres, 2, vs_mean, k = keep ), las = 2, col = 'yellow' )

round( res$MCMC$birth_acceptances, 2 )
round( res$MCMC$death_acceptances, 2 )
round( res$MCMC$all_acceptances, 2 )

res$MCMC$birth_oppurtunity
res$MCMC$birth
res$MCMC$death_oppurtunity
res$MCMC$death

### Beta traceplots
plot( res$MCMC$betares[keep, 1], type = 'l' )
hline( 0 )
plot( density(res$MCMC$betares[keep, 1]) )
vline( 0 )
plot( res$MCMC$betares[keep, 2], type = 'l' )
hline( 0.8 )
plot( density(res$MCMC$betares[keep, 2]) )
vline( 0.8 )
plot( res$MCMC$betares[keep, 3], type = 'l' )
hline( -0.8 )
plot( density(res$MCMC$betares[keep, 3]) )
vline( -0.8 )

### DLF parameter traceplots
plot( res$MCMC$DLFres$num_1[keep, 1], type = 'l' )
hline( 0 )
plot( density(res$MCMC$DLFres$num_1[keep, 1]) )
vline( 0 )
plot( res$MCMC$DLFres$num_1[keep, 2], type = 'l' )
hline( -0.1 )
plot( density(res$MCMC$DLFres$num_1[keep, 2]) )
vline( -0.1 )

plot( res$MCMC$DLFres$num_2[keep, 1], type = 'l' )
hline( 0 )
plot( density(res$MCMC$DLFres$num_2[keep, 1]) )
vline( 0 )
plot( res$MCMC$DLFres$num_2[keep, 2], type = 'l' )
hline( -0.1 )
plot( density(res$MCMC$DLFres$num_2[keep, 2]) )
vline( -0.1 )


### Model selection probability (gamma prior parameter)
plot( res$MCMC$model_sel_prob[keep], ylim = c(0, 1), type = 'l' )
hline( 1 )
plot( density(res$MCMC$model_sel_prob[keep]), xlim = c(0, 1) )
vline( 1 )

### DLF parameter acceptance rates
round( res$MCMC$DLF_accept / keep[length(keep)], 2 )[1:2]

##### Correlation Matrices Image plots
corrplot( cor(na.omit( res$train_dat[, 1:nvars])), method = 'number' )

image( t( as.matrix.csr(res$MCMC$vsres[keep, ]) ) )

im_compare <- function( ind = 2, actual = 2, keep = seq( 50000, 100000, 10 ) ) {
  par(mfrow = c(2, 1))
  image( as.matrix.csr(t(res$MCMC$vsres[keep,])) )
  plot( res$MCMC$betares[keep, ind], type = 'l', ylab = '', xaxs = 'i', yaxs = 'i' )
  hline( actual )
  par(mfrow = c(1, 1))
}

im_compare(2, 2)
im_compare(3, -2)

### Theta covariance
pdf('output/DLF1_posterior_covariance.pdf')
plot( res$MCMC$DLFres$num_1[keep, 1], log( -res$MCMC$DLFres$num_1[keep, 2] ) )
hline( -1 )
vline( 10 )
dev.off()
cor( res$MCMC$DLFres$num_1[keep, 1], log( -res$MCMC$DLFres$num_1[keep, 2] ) )
cov( res$MCMC$DLFres$num_1[keep, 1], log( -res$MCMC$DLFres$num_1[keep, 2] ) )
var( res$MCMC$DLFres$num_1[keep, 1] )
var( log( -res$MCMC$DLFres$num_1[keep, 2] ) )

plot( res$MCMC$DLFres$num_2[keep, 1], res$MCMC$DLFres$num_2[keep, 2] )
hline( -1 )
vline( 10 )
cor( res$MCMC$DLFres$num_2[keep, 1], res$MCMC$DLFres$num_2[keep, 2] )
cov( res$MCMC$DLFres$num_2[keep, 1], -exp(res$MCMC$DLFres$num_2[keep, 2]) )
var( res$MCMC$DLFres$num_2[keep, 1] )
var( res$MCMC$DLFres$num_2[keep, 2] )

#### What does nealmon look like?
xx <- seq(0.1, 10, 0.1)
yy1 <- irts_nealmon(xx, c(6, -0.6))
yy2 <- irts_nealmon(xx, c(18, -1.8))
yyreal <- irts_nealmon(xx, c(10, -1))

pdf('output/DLF1_posterior_covariance_nealmon.pdf')
plot(xx, yy2, type = 'l', col = 'dodgerblue')
lines(xx, yy1, type = 'l', col = 'orange')
lines(xx, yyreal, type = 'l', col = 'green')
dev.off()
}
