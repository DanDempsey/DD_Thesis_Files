##### Sim Example for Thesis
##### Daniel Dempsey

source('package_code/1_load_all.R')

library( parallel )
library( tictoc )

binom_innov <- function(n, x) {
  fam <- binomial()
  rbinom(n, prob = fam$linkinv(x), size = 1)
}

sim_test <- function( b, t ) {
  
  nvar <- length(b) - 1
  test <- irts_midas_sim( n_y = 1000, n_num = nvar, seasonal_adjust = TRUE,
                          beta = b, innov_fun = binom_innov,
                          DLF_parameters = t )$Data
  
  mrtest <<- midas_regimes(value = test[1:nvar], time = test$INDEX, DLF_parameters = list(c(0, -0.1)))
  rtest <<- response_regime(value = test$y, time = test$INDEX)
  form <- rtest ~ mrtest
  
  fit <- irts_midas( form, data = test, family = "binomial", 
                     method = "BFGS", gr = midas_gradient )
  
  list( fit = fit, sim_dat = test, b = b, t = t )
    
}

b_list <- list( c( 0, 0.8, -0.8 ),
                c( -2, 0.5, -0.5, 0.8, 0, 0 ) )

b_list_all <- rep( b_list, each = 30 )

t_list <- list( list( c( 14, -1 ), c( 0.5, -0.1 ) ), 
                list( c( 14, -1 ), c( 0.5, -0.1 ), c( 14, -1 ), 
                      c( 0.5, -0.1 ), c( 14, -1 ) ) ) 

t_list_all <- rep( t_list, each = 30 )

set.seed( 1608 )
tic()
all_irts_midas_sims <- mcMap( sim_test, b = b_list_all, t = t_list_all,
                              mc.cores = 6 )
toc() # 365.594 sec elapsed

save( all_irts_midas_sims, file = 'package_code/Simulations/IRTS_MIDAS_Sims/fits.Rdata' )
load( 'package_code/Simulations/IRTS_MIDAS_Sims/fits.Rdata' )

### Results
# beta
beta_extract <- function( x, ind ) {
  
  x$fit$par[ ind ]
  
}

true1 <- c( 0, 0.8, -0.8 )

pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/Beta_sim1.pdf', height = 5, width = 15 )
par( mfrow = c( 1, 3 ) )
for ( i in 1:3 ) {
  
  
  bet <- sapply( all_irts_midas_sims[1:30], beta_extract, ind = i )
  plot( density(bet), xlab = paste0( 'beta ', i - 1 ), main = '' )
  abline( v = true1[i], col = 'dodgerblue' )
  
}
dev.off()

true2 <- c( -2, 0.5, -0.5, 0.8, 0, 0 )

pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/Beta_sim2.pdf', height = 10, width = 15 )
par( mfrow = c( 2, 3 ) )
for ( i in 1:6 ) {
  
  bet <- sapply( all_irts_midas_sims[31:60], beta_extract, ind = i )
  plot( density(bet), xlab = paste0( 'beta ', i - 1 ), main = '' )
  abline( v = true2[i], col = 'dodgerblue' )
  
}
dev.off()

# theta
pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/Theta_sim11.pdf', width = 10 )
t11 <- sapply( all_irts_midas_sims[1:30], function( x ) { x$fit$par[4] } )
plot( density(t11), xlab = 'Theta 1', main = '' )
abline( v = 14, col = 'dodgerblue' )
dev.off()

pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/Theta_sim1mode.pdf', width = 10 )
t12 <- sapply( all_irts_midas_sims[1:30], function( x ) { x$fit$par[5] } )
plot( density( t11 / (-2 * t12) ), xlab = 'Nealmon Mode', main = '' )
abline( v = 7, col = 'dodgerblue' )
dev.off()

### Example covariate graph
pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/sim_dat_ex.pdf', width = 10 )
plot( na.omit( all_irts_midas_sims[[1]]$sim_dat$num_1 ), type = 'l', 
      ylab = '', main = 'Example Simulated Covariate' )
dev.off()

### Example Linear Predictor
ex_dat <- midas_design_matrices( all_irts_midas_sims[[1]]$fit$regime_object )[ , 1 ]
pdf( 'package_code/Simulations/IRTS_MIDAS_Sims/Graphs/lin_pred.pdf', width = 10 )
plot( ex_dat, type = 'l' )
dev.off()
