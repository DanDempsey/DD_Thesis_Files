##### Simulate Mixed Regime Data and IRTS-MIDAS response variable
##### Daniel Dempsey

library( dplyr )

### Function to make step adjustments to covariates
step_adjust_fun <- function( x, step_seq, step_len, step_size ) {
  
  step_add <- step_seq + step_len
  steps <- rep( seq( step_seq + 1, length(x)-step_add, step_add ), each = step_len ) + 
    seq(0, step_len-1, 1)
  x <- x - step_size
  x[steps] <- x[steps] + 2*step_size
  x
  
}

### Function to make seasonal adjustments to covariates
seasonal_adjust_fun <- function( x, p, s ) {
  
  x + (2*cos( (2*pi/p)*(seq_along(x)+s) ))
  
}


### Simulation Code
irts_midas_sim <- function(n_y = 1000, n_vars = 1, 
                           beta = NULL, time_window = 30, DLF = "irts_nealmon",
                           DLF_parameters = list( c(14, -1) ), pars = NULL, sigma_base = NULL,
                           response_index = NULL, covar_sample_index_prop = 0.5, 
                           response_rate = 0.25, innov_fun = rnorm, innov_pars = NULL, 
                           step_adjust = FALSE, seasonal_adjust = FALSE, periodicity = 30, shift = 7,
                           which_adjust, ...) {
  
  ### Prepare simulation parameters
  if ( missing(which_adjust) ) {
    which_adjust <- 1:n_vars
  }
  var_seq <- 1:n_vars
  if (is.null(pars)) { pars <- diag(0.5, n_vars) }
  if (!is.matrix(pars)) { pars <- diag(pars, n_vars) }
  #power_elements <- abs( rep(var_seq, n_vars) - rep(var_seq, each = n_vars) )
  #power_mat <- matrix(power_elements, nrow = n_vars)
  if(is.null(sigma_base)) { sigma_base <- diag(1, n_vars) }
  if( all(dim(sigma_base) == 1) ) { sigma_base <- as.vector(sigma_base) } 
  
  ### Simulate response indices
  if (is.null(response_index)) {
    response_index <- cumsum( rexp(n_y, response_rate) ) + max(time_window)
  }
  n_x_full <- max(response_index)
  
  ### Simulate data and assign appropriate names
  X_full <- simulateVAR(pars = pars, Nt = n_x_full, residuals = sigma_base, ...)
  if( step_adjust ) {
    X_full[which_adjust] <- Map( step_adjust_fun, x = X_full[which_adjust],
                                 step_seq = time_window, step_len = time_window,
                                 step_size = 5 )
    X_full <- as.data.frame( X_full )
  }
  if( seasonal_adjust ) {
    #X_full[which_adjust] <- lapply( X_full[which_adjust], seasonal_adjust_fun, p = periodicity )
    X_full <- as.data.frame( Map( seasonal_adjust_fun, x = X_full, p = periodicity, s = shift ) )
  }
  names(X_full) <- paste0("num_", var_seq)
  
  ### Retrieve covariate indices
  INDEX <- 1:n_x_full
  X <- as.data.frame( apply(X_full, 2, scale) )
  X$INDEX <- X_full$INDEX <- INDEX
  
  ### Now compute the linear predictor
  regime0 <- midas_regimes(value = X[var_seq], time = X$INDEX, time_window = time_window,
                           DLF = DLF, DLF_parameters = DLF_parameters)
  
  # Create regime object
  regime1 <- lapply( regime0, "[[<-", i = "response_index", value = response_index )
  regime <- lapply( regime1, time_delta_list )

  # Compute the midas design matrix and thus retrieve the linear predictor
  WX <- cbind( 1, midas_design_matrices(regime) )
  WX_p <- ncol(WX)
  lin_pred <- WX %*% beta
  
  ### Now simulate the y values based on lin_pred
  innov_pars_complete <- c(n = n_y, x = list(lin_pred), innov_pars) 
  y <- do.call("innov_fun", innov_pars_complete)
  
  ### Compile all the simulated data and return the result
  response_dat <- data.frame(INDEX = response_index, y = y)
  final_dat <- full_join(X, response_dat, by = "INDEX")
  
  list( Data = final_dat[ order(final_dat$INDEX), ], 
        fit_parameters = list(beta = beta, DLF_parameters = DLF_parameters) )
  
}

if (FALSE) {
  
  ALD_innov <- function( n, x, q = 0.5 ) {
    ifelse( rALD( n, mu = x, p = q ) >= 0, 1, 0 )
  }
  
  binom_innov <- function( n, x ) {
    rbinom( n, prob = binomial()$linkinv( x ), size = 1 )
  }
  
  set.seed( 101 )
  nvar <- 2
  tw <- 10
  qq <- 0.5
  mc_len <- 5000
  bi <- 2500
  thn <- 10
  keep <- seq( bi, mc_len, thn )
  ri <- seq( tw, 5000, tw/2 )
  xx <- irts_midas_sim( n_y = length(ri), time_window = tw, n_vars = nvar, pars = diag(0.05, nvar),
                        beta = c( 0, 2, -2 ), innov_fun = ALD_innov, innov_pars = list(q = qq),
                        DLF_parameters = list( c(14, -1) ), seasonal_adjust = TRUE,
                        step_adjust = FALSE, response_index = ri, which_adjust = 1:nvar )
  
  cor( xx$Data[1:nvar] )
  
  plot( na.omit( xx$Data$num_1[1:500] ), type = 'l' )
  points( xx$Data$y )
  
  mrtest <- midas_regimes( xx$Data[1:nvar], xx$Data$INDEX, DLF_parameters = c(0.5, -0.1),
                           time_window = 30 )
  rtest <- response_regime( xx$Data$y, xx$Data$INDEX )
  form <- rtest ~ mrtest
  
  test_fit <- IRTS_MIDAS_AuxVar( formula = form, data = xx$Data, quantile = qq,
                                 response_dist = "ald", varsel = TRUE,
                                 MCMC_length = mc_len )
  
  #XX <- midas_design_matrices( test_fit$regime_object )
  #plot( XX[, 1], type = 'l', ylim = range(na.omit( xx$Data$num_1 )) )
  
  hist( test_fit$betares[ keep, 1 ] )
  ppd_ALD_vs( test_fit, form, burn_in = bi, thin = thn )
  
}

##### Function explanation:
## Simulate a multivariate time series using the simulateVAR function from the mlVAR package.
## If factor variables were requested, they are created by bucketing numerical values by their quantiles.
## The response index is simulated as the cumulative sum of exponential random variables.
## The value of the response is created based on the MIDAS model.
## Output: A list: The first element is a data frame containing an index column, covariate columns, and response columns.
##         The second element says whether the response was generated from the data

##### Function arguement explanations:
## n_y: Number of response observations
## n_x_full: Length of multivariate time series sampled by mlVAR::simulateVAR
## n_num: Number of numerical covariates
## n_fac: Number of factor covariates
## lev_threshold_quantile: Factor levels are simulated by first simulating a numeric variable, then bucketing
##                         on the quantiles. This parameter tells the function where the thresholds should be.
## beta: beta parameters. The first value is assumed to be the intercept. The next values are assumed to correspond to
##       the numerical covariates. The ones after those should correspond to the factor levels.
## time_window: Size of time window to be used when creating response
## DLF: Distributed Lag Functions to be used when creating response
## DLF_parameters: Parameters for the DLF
## pars: The parameter arguement for the mlVAR::simulateVAR function. See ?mlVAR::simulateVAR for more details
## sigma_base: The base of the covariance matrix for the simulated covariate distribution
## covar_index: User supplied indices for the covariates
## response_index: User supplied indices for the response
## covar_sample_index_prop: The proportion of the the VAR that is sampled if the user did not supply covar_index
## The rate parameter of the exponential distribution if the user did not supply response_index
## innov_fun: How the response is randomly generated after the linear predictor is computed. For this to work,
##            the following format must be adhered to: the first arguement should be the length of the data,
##            and the second argument is the mean
## innov_pars: Extra arguements for the innovation function if required
## post_sample_response: If TRUE, response values are based on variables AFTER sampling the covariate indices.
##                       Otherwise, use the entire VAR simulation values.
## ...: Extra arguements for the mlVAR::simulateVAR function.

