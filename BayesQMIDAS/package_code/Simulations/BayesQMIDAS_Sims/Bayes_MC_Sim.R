##### Big Sim Code
##### Daniel Dempsey

source( 'package_code/1_load_all.R' )
setwd( 'package_code/Simulations/BayesQMIDAS_Sims/' )
library( parallel )
library( corrplot )
library( tictoc )

ALD_innov <- function(n, x, q = 0.5) {
  ifelse( rALD(n, mu = x, p = q) >= 0, 1, 0 )
}

run_sim <- function(b = c(0, -2, 2), DLF = list(c(14, -1)), q = 0.5, ny = 1000, nx = 5000,
                    MCMC_len = 50000, burn_in = 25000, thin = 10, time_window = 30, 
                    reps = 1, n_cores, ...) {
  
  # Simulate data
  n_vars <- length(b) - 1
  #if ( missing(which_adjust) ) { which_adjust <- 1:n_vars }
  sim <- irts_midas_sim(n_y = ny, n_x_full = nx, innov_fun = ALD_innov, 
                        n_num = n_vars, beta = b, DLF_parameters = DLF,
                        seasonal_adjust = TRUE, innov_pars = list(q = q),
                        time_window = time_window )
  sim_dat <- sim$Data
  
  # Split into training/test sets
  y_full <- na.omit( sim_dat[, c('INDEX', 'y')] )
  
  INDEX_train <- y_full[ 1:(0.8*ny), 'INDEX' ]
  sim_dat_train <- sim_dat[ which( sim_dat$INDEX <= max( INDEX_train ) ), ]
  
  INDEX_test <- setdiff( y_full$INDEX, INDEX_train )
  sim_dat_test <- sim_dat[ which( sim_dat$INDEX >= ( min( INDEX_test ) - time_window ) ), ]
  drop_inds <- which( !is.na( sim_dat_test$y ) & ( sim_dat_test$INDEX < min( INDEX_test ) )  )
  sim_dat_test <- sim_dat_test[ -drop_inds, ]
  
  imba <- sum( y_full$y ) / length( y_full$y )
  rm( y_full, INDEX_train, INDEX_test, drop_inds )
  
  # Fit the model
  mrtest <- midas_regimes( value = sim_dat_train[, seq(1, n_vars)], time = sim_dat_train$INDEX, 
                           DLF_parameters = list(c(0, -1)), time_window = time_window )
  rtest <- response_regime( value = sim_dat_train$y, time = sim_dat_train$INDEX )
  form <- rtest ~ mrtest
  
  cor_dat = cor( na.omit( sim_dat_train[, 1:n_vars] ) )
  
  corrplot( cor_dat, method = 'number' )
  
  run_fun <- function( xx ) {
    IRTS_MIDAS_AuxVar( formula = form, data = sim_dat_train, quantile = q,
                       response_dist = "ald", varsel = TRUE,
                       MCMC_length = MCMC_len )
  }
  res <- mclapply( 1:reps, run_fun, mc.cores = n_cores )
  
  # Training ppd
  ppd_res_train <- mclapply( res, ppd_ALD_vs, new_form = form, 
                             burn_in = burn_in, thin = thin, 
                             mc.cores = n_cores )
  
  # Test ppd
  mrtest_test <- midas_regimes( value = sim_dat_test[, seq(1, n_vars)], 
                                time = sim_dat_test$INDEX, 
                                DLF_parameters = list(c(0, -1)), 
                                time_window = time_window )
  rtest_test <- response_regime( value = sim_dat_test$y, time = sim_dat_test$INDEX )
  form_test <- rtest_test ~ mrtest_test
  
  ppd_res_test <- mclapply( res, ppd_ALD_vs, new_form = form_test, 
                            burn_in = burn_in, thin = thin, 
                            mc.cores = n_cores )
  
  # Compile results
  list( MCMC = res, ppd_train = ppd_res_train, ppd_test = ppd_res_test, 
        train_dat = sim_dat_train, test_dat = sim_dat_test, imbalance = imba,
        data_correlation = cor_dat, true_beta = b, true_DLF = DLF, q = q )
  
}

run_sim_loop_wrapper <- function( x, ... ) {
  run_sim( ... )
}

### Set simulation parameters
b_test1 <- c(0, 0.8, -0.8)
b_test2 <- c(0, 0.2, -0.2, 0.5, -0.5, 0.8, rep(0, 15))

blen1 <- length( b_test1 )
blen2 <- length( b_test2 )

myprior1 <- list( beta0 = rep(0, blen1), V0 = 100 * diag(blen1), 
                  DLF_pars = list( DLF1 = c(15, 10), DLF2 = c(1, 1) ),
                  model_selection = c(1, 1) )
myprior2 <- list( beta0 = rep(0, blen2), V0 = 100 * diag(blen2), 
                  DLF_pars = list( DLF1 = c(15, 10), DLF2 = c(1, 1) ),
                  model_selection = c(1, 1) )

DLF1 <- list( c(14, 1), c(0.5, -0.1) )
DLF2 <- rep( list( c(14, 1), c(0.5, -0.1) ), 10 )

reps <- 10
ncore <- 5
quants <- c( 0.5, 0.2 )

### Run simulations
set.seed( 222 )
tic()
for ( i in 1:2 ) {
  
  for ( j in 1:2 ) {
    
    b_test <- get( paste0('b_test', i) )
    prior_test <- get( paste0('myprior', i) )
    DLF_test <- get( paste0('DLF', i) )
    q_test <- quants[ j ]
    
    res <- run_sim( b = b_test, prior = prior_test, DLF = DLF_test, 
                    q = q_test, n_cores = ncore, reps = reps,
                    MCMC_len = 100000, burn_in = 50000 )
    
    save( res, file = paste0( 'res_', i, j, '.Rdata' ) )
    rm( res )
    
  }
  
}
toc() # 37339.404 sec elapsed / roughly 10 hours

