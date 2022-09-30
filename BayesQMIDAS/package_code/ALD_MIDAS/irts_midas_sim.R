##### Simulate Mixed Regime Data and IRTS-MIDAS response variable
##### Daniel Dempsey

### Function to make step adjustments to covariates
step_adjust_fun <- function( x, step_seq, step_len, step_size ) {
  
  steps <- rep( seq( step_seq, length(x)-step_len, step_seq ), each = step_len ) + 
    seq(0, step_len-1, 1)
  x[steps] <- step_size + x[steps]
  x
  
}

### Function to make seasonal adjustments to covariates
seasonal_adjust_fun <- function( x ) {
  
  len <- runif( 1, 0.01, 0.1 )
  amp <- runif( 1, 1, 3 )
  x + ( amp * cos( len * seq_along(x) ) )
  
}


### Simulation Code
irts_midas_sim <- function(n_y = 1000, n_x_full = 5000, n_num = 1, n_fac = 0, lev_threshold_quantile = list(0.5), 
                           beta = NULL, time_window = 30, DLF = "irts_nealmon",
                           DLF_parameters = list( c(14, -1) ), pars = NULL, sigma_base = 0,
                           response_index = NULL, covar_sample_index_prop = 0.5, 
                           response_rate = 0.25, innov_fun = rnorm, innov_pars = NULL, post_sample_response = TRUE, 
                           #step_adjust = FALSE, 
                           seasonal_adjust = FALSE, which_adjust, 
                           ...) {
  
  ### Prepare simulation parameters
  n_vars <- sum(n_num, n_fac)
  if ( missing(which_adjust) ) {
    which_adjust <- 1:n_vars
  }
  var_seq <- 1:n_vars
  if (is.null(pars)) { pars <- diag(0.5, n_vars) }
  if (!is.matrix(pars)) { pars <- diag(pars, n_vars) }
  power_elements <- abs( rep(var_seq, n_vars) - rep(var_seq, each = n_vars) )
  power_mat <- matrix(power_elements, nrow = n_vars)
  covar_mat <- ( sigma_base ^ power_mat )
  if( all(dim(covar_mat) == 1) ) { covar_mat <- as.vector(covar_mat) } 
  
  ### Simulate data and assign appropriate names
  X_full <- simulateVAR(pars = pars, Nt = n_x_full, residuals = covar_mat, ...)
  #if( step_adjust ) {
  #  step_starts <- sample(seq(50, 200, 10), replace = TRUE, n_vars)
  #  step_lens <- sample(seq(10, 20, 1), replace = TRUE, n_vars)
  #  step_sizes <- sample(c(-3, 3), replace = TRUE, n_vars)
  #  X_full[which_adjust] <- Map( step_adjust_fun, x = X_full[which_adjust], 
  #                               step_seq = step_starts, step_len = step_lens, 
  #                               step_size = step_sizes )
  #  X_full <- as.data.frame( X_full )[which_adjust]
  #}
  if( seasonal_adjust ) {
    X_full[which_adjust] <- lapply( X_full[which_adjust], seasonal_adjust_fun )
    X_full <- as.data.frame( X_full )
  }
  if( n_num > 0 ) { num_names <- paste0("num_", 1:n_num) }
  else { num_names <- character(0) }
  if( n_fac > 0 ) { fac_names <- paste0("fac_", 1:n_fac) }
  else { fac_names <- character(0) }
  names(X_full) <- c(num_names, fac_names)
  
  ### Convert to factor if user requested
  if (n_fac > 0) {
    
    fac_var_inds <- (n_vars - n_fac + 1):n_vars
    fac_vars <- X_full[fac_var_inds]
    lev_threshold_quantile <- as.list(lev_threshold_quantile)
    if (length(lev_threshold_quantile) > n_fac) { 
      warning("Too many lev_threshold_quantile elements. Truncating list to match the number of factor variables.")
      lev_threshold_quantile <- lev_threshold_quantile[1:n_fac]
    }
    
    # Bucket data based on their quantile
    quantile_breaks <- Map( quantile, x = fac_vars, probs = lev_threshold_quantile )
    quantile_breaks <- Map( function(x, y) { c( min(y), x, max(y) ) }, x = quantile_breaks, y = fac_vars)
    labs <- lapply(lev_threshold_quantile, function(x) { 0:length(x) })
    fac_vars_bucketed <- Map( cut, x = fac_vars, breaks = quantile_breaks, labels = labs,
                                 MoreArgs = list(include.lowest = TRUE))
    
    # Append this to original data
    X_full[fac_var_inds] <- fac_vars_bucketed
    
  }
   
  ### Retrieve covariate indices
  INDEX <- 1:n_x_full
  #if (is.null(covar_index)) {
  #  size = floor(n_x_full * covar_sample_index_prop)
  #  covar_index <- replicate( n_vars, sort( sample(x = INDEX, size = size) ), simplify = FALSE )
  #}
  
  covar_index <- rep( list(INDEX), n_vars )
  if (length(covar_index) > n_vars) { 
    warning("Too many covar_index elements. Truncating list to match the number of variables.")
    covar_index <- covar_index[var_seq]
  }
  
  # Convert observations at the non-chosen indices to NA
  anti_covar_ind <- Map(setdiff, y = covar_index, MoreArgs = list(x = INDEX))
  X_list <- Map("[<-", x = X_full, i = anti_covar_ind, MoreArgs = list(value = NA) )
  X <- do.call("data.frame", X_list)
  X <- as.data.frame( lapply(X, scale) )
  X$INDEX <- X_full$INDEX <- INDEX
  res <- list( Data = X, Full = X_full, model_data = ifelse(post_sample_response, "partial", "full") )
  
  ### Now compute the linear predictor
  mod_dat <- res[[ ifelse(post_sample_response, 1, 2) ]]
  regime0 <- midas_regimes(value = mod_dat[var_seq], time = mod_dat$INDEX, time_window = time_window,
                           DLF = DLF, DLF_parameters = DLF_parameters)
  
  # Simulate the repsonse index and find the time windows
  if (is.null(response_index)) {
    response_index <- list( All = n_x_full + 1 )
    while ( max(response_index$All) > n_x_full ) { # To avoid having a larger index than is available
      response_index <- list( All = cumsum( rexp(n_y, response_rate) ) + max(time_window) )
    }
  }
  regime1 <- lapply( regime0, "[[<-", i = "response_index", value = response_index )
  regime <- lapply( regime1, time_delta_list )

  # Compute the midas design matrix and thus retrieve the linear predictor
  WX <- cbind( 1, midas_design_matrices(regime) )
  WX_p <- ncol(WX)
  if (is.null(beta)) { beta <- rnorm( WX_p, mean = 0, sd = 3 ) }
  if (length(beta) > WX_p) {
    warning("Too many beta elements. Truncating variables to match the number of design matrix columns.")
    beta <- beta[1:WX_p]
  }
  if (length(beta) < WX_p) { stop(paste0("Too few beta elements. Missing ", WX_p - length(beta), " values.")) }
  lin_pred <- WX %*% beta
  
  ### Now simulate the y values based on lin_pred
  innov_pars_complete <- c(n = n_y, x = list(lin_pred), innov_pars) 
  y <- try(do.call("innov_fun", innov_pars_complete), silent = TRUE)
  if (inherits(y, "try-error"))
    stop( paste0("The innovation function failed with the following message:\n", y[1]) )
  
  ### Compile all the simulated data and return the result
  response_dat <- data.frame(INDEX = unlist(response_index), y = y)
  final_dat <- dplyr::full_join(X, response_dat, by = "INDEX")
  res[[1]] <- final_dat[ order(final_dat$INDEX), ]
  res$fit_parameters <- list(beta = beta, DLF_parameters = DLF_parameters)
  
  res
  
}

if (FALSE) {
  
  ALD_innov <- function( n, x, q = 0.5 ) {
    ifelse( rALD( n, mu = x, p = q ) >= 0, 1, 0 )
  }
  
  binom_innov <- function( n, x ) {
    rbinom( n, prob = binomial()$linkinv( x ), size = 1 )
  }
  
  set.seed( 100 )
  nvar <- 2
  xx <- irts_midas_sim(n_y = 1000, n_num = nvar, pars = diag(0.95, nvar),
                 beta = c( 0, 0.8, -0.8 ), innov_fun = ALD_innov,
                 DLF_parameters = list( c(14, -1) ), seasonal_adjust = TRUE )
  
  #plot( na.omit( xx$Data$num_1 ), type = 'l' )
  
  mrtest <- midas_regimes( xx$Data[1:nvar], xx$Data$INDEX, DLF_parameters = c(0.5, -0.1) )
  
  rtest <- response_regime( xx$Data$y, xx$Data$INDEX )
  
  form <- rtest ~ mrtest
  
  test_fit <- IRTS_MIDAS_AuxVar(formula = form, data = xx$Data, quantile = 0.5,
                    response_dist = "ald", varsel = TRUE,
                    MCMC_length = 1000 )
  
  #XX <- midas_design_matrices( test_fit$regime_object )
  #plot( XX[, 1], type = 'l', ylim = range(na.omit( xx$Data$num_1 )) )
  
  ppd_ALD_vs( test_fit, form, burn_in = 500 )
  
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

