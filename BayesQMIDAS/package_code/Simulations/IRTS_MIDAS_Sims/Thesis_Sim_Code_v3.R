##### Sim Example for Thesis
##### Daniel Dempsey

source( 'package_code/1_load_all.R' )
setwd( 'package_code/Simulations/IRTS_MIDAS_Sims/Output' )

library( parallel )
library( tictoc )

binom_innov <- function(n, x) {
  rbinom(n, prob = binomial()$linkinv(x), size = 1)
}

sim_test <- function( r, b, t, n_s, tw_sim = 15, tw_mod = 15, 
                      rr = 1/15, sb = 1, DLFm = "irts_nealmon", DLFs = "irts_nealmon",
                      DLFp = list(c(0, -0.1)) ) {
  
  nvar <- length(b) - 1
  test <- irts_midas_sim( n_y = n_s, n_vars = nvar, seasonal_adjust = TRUE,
                          beta = b, innov_fun = binom_innov, time_window = tw_sim, shift = c(7, 2.5),
                          DLF = DLFs, DLF_parameters = t, response_rate = rr, sigma_base = diag(sb, nvar) )$Data
  
  mrtest <<- midas_regimes( value = test[1:nvar], time = test$INDEX, 
                            DLF_parameters = DLFp, time_window = tw_mod,
                            DLF = DLFm )
  rtest <<- response_regime(value = test$y, time = test$INDEX)
  form <- rtest ~ mrtest
  
  if( DLFm == 'irts_nealmon' ) {
    gr_mod <- midas_gradient
    method_mod <- 'BFGS'
    lower_mod <- -Inf
  }
  else {
    gr_mod <- NULL
    method_mod <- 'L-BFGS-B'
    lower_mod <- 0 # <- bug
  }
  
  irts_midas( form, data = test, family = "binomial", method = method_mod, 
              gr = gr_mod, lower = lower_mod )$par
  
  #fit$regime_object <- NULL # to save space
  
  #list( fit = fit, b = b, t = t, ns = n_s, tw_sim = tw_sim, tw_mod = tw_mod, 
        #rr = rr, sb = sb, DLF_mod = DLF_mod )
  
}

density_output <- function( x, trup, file_path ) {
  
  mn <- apply( x, 2, mean )
  nm <- colnames( x )
  for ( i in 1:ncol(x) ) {
    
    pdf( paste0(file_path, '/density_', nm[i], '.pdf') )
    plot( density(x[, i]), main = nm[i] )
    abline( v = trup[i], lty = 2, col = 'orange' )
    abline( v = mn[i], lty = 2, col = 'dodgerblue' )
    dev.off()
    
  }
  
  theta_var1 <- -x[, 4] / (2 * x[, 5]) 
  pdf( paste0(file_path, '/density_theta_var1.pdf') )
  plot( density(theta_var1), main = 'Distributed Lag Variable 1' )
  abline( v = 7, lty = 2, col = 'orange' )
  abline( v = mean(theta_var1), lty = 2, col = 'dodgerblue' )
  dev.off()
  
  theta_var2 <- -x[, 6] / (2 * x[, 7]) 
  pdf( paste0(file_path, '/density_theta_var2.pdf') )
  plot( density(theta_var2), main = 'Distributed Lag Variable 2' )
  abline( v = 2.5, lty = 2, col = 'orange' )
  abline( v = mean(theta_var2), lty = 2, col = 'dodgerblue' )
  dev.off()
  
}

sim_output <- function( reps = 1000, n_samps = c(100, 500, 1000), file_path,
                        betap = c( 0, 1.5, -1.5 ), 
                        thetap = list( c( 14, -1 ), c( 0.5, -0.1 ) ), 
                        time_wind_sim = 15, time_wind_mod = 15, 
                        response_rate = 1/15, sigma_base = 1, DLF_mod = "irts_nealmon", 
                        DLF_sim = "irts_nealmon", DLF_pars = list(c(0, -0.1)), n_cores = 6 ) {
  
  len <- length( n_samps )
  res <- vector( 'list', len )
  names( res ) <- n_samps
  
  const_args <- list( b = betap, t = thetap, tw_sim = time_wind_sim,
                      tw_mod = time_wind_mod, rr = response_rate, sb = sigma_base,
                      DLFm = DLF_mod, DLFs = DLF_sim, DLFp = DLF_pars )
  
  for ( i in 1:len ) {
    const_args$n_s <- n_samps[i]
    res[[i]] <- do.call( 'rbind', mcMap( sim_test, r = 1:reps,
                                         MoreArgs = const_args, 
                                         mc.cores = n_cores ) )
    colnames( res[[i]] ) <- c( 'Intercept', 'beta1', 'beta2', 'var1theta1',
                               'var1theta2', 'var2theta1', 'var2theta2' )
  }
  
  true_vals <- c( betap, unlist(thetap) )
  if( !dir.exists(file_path) ) { dir.create(file_path) }
  setwd( file_path )
  dir_path <- paste0('res_', names(res))
  for ( i in 1:length(res) ) {
    if( !dir.exists(dir_path[i]) ) { dir.create(dir_path[i]) }
  }
  
  Map( density_output, x = res, file_path = dir_path, MoreArgs = list(trup = true_vals) )
  
  res$input <- const_args
  save( res, file = 'res.Rdata' )
  setwd( '..' )
  
}

### Test
set.seed(11)
tic()
sim_output( reps = 6, n_cores = 1, file_path = 'Test' )
toc() # 43.988 sec elapsed

gc()

set.seed(44)
tic()
sim_output( reps = 1000, n_cores = 6, response_rate = 1/7, file_path = 'Irreg7' )
toc() # 5814.186 sec elapsed

gc()

### Change Noise
set.seed(55)
tic()
sim_output( reps = 1000, n_cores = 6, sigma_base = 3, file_path = 'Noise' )
toc() # 7130.912 sec elapsed

gc()

### Change Time Window
set.seed( 66 )
tic()
sim_output( reps = 1000, n_cores = 6, time_wind_mod = 7, file_path = 'TW7' )
toc() # 3826.518 sec elapsed

gc()

### Extra
set.seed( 88 )
tic()
sim_output( reps = 1000, n_cores = 6, response_rate = 1/30, file_path = 'Irreg30' )
toc() # 9723.695 sec elapsed

gc()

set.seed( 99 )
tic()
sim_output( reps = 1000, n_cores = 6, time_wind_mod = 25, file_path = 'TW25' )
toc() # 8178.584 sec elapsed

gc()

