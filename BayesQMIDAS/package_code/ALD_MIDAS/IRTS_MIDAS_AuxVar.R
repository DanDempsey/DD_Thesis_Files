###### Auxiliary Variable Binary Inference
###### Daniel Dempsey

IRTS_MIDAS_AuxVar <- function(formula, data, response_dist = "ald", quantile = 0.5, prior, 
                              beta_start, gamma_start, incl_DLF_pars = TRUE, varsel = FALSE, 
                              MCMC_length = 10000L, ...) {
  
  ### Initialise environment to be passed into model.frame
  if ( missing(data) ) { Zenv <- environment() }
  else { Zenv <- as.environment(data) }
  parent.env(Zenv) <- environment(formula)
  
  ### Check response_dist is valid
  valid_settings <- c("normal", "logistic", "ald", "error")
  chosen_setting <- pmatch(response_dist, valid_settings, 4)
  response_dist <- valid_settings[chosen_setting]
  if ( response_dist == "error" ) {
    stop("Supported values for response_dist are 'normal', 'logistic' and 'ald'.")
  }
    
  ### Check if DLF parameter inference is desired
  if ( !incl_DLF_pars ) {
    response_dist <- paste0( response_dist, '_fixedDLF' )
  }
    
  ### Check if variable selection is desired
  if ( varsel ) { response_dist <- paste0( response_dist, '_vs' ) }
  
  ### Unpack the formula
  form <- midas_formula_unpack(formula, Zenv)
  
  ### Check prior
  if ( missing(prior) ) {
    prior <- get(paste0("defaultPrior_", response_dist))(form, Zenv)
  }
  
  if( is.null(prior) ) {
    prior <- get(paste0("defaultPrior_", response_dist))(form, Zenv)
  }
    
  ### Check starting values
  # Beta values
  n_var <- length( prior$beta0 )
  if ( missing(beta_start) ) { beta_start <- prior$beta0 }
  if ( length(beta_start) < n_var ) {
    warning("Too few values in beta_start. Using mean of prior distribution instead.")
    beta_start <- prior$beta0
  }
  if ( length(beta_start) > n_var ) {
    warning("Too many values in beta_start. Truncating to suitable length.")
    beta_start <- beta_start[1:n_var]
  }
  
  if ( missing(gamma_start) ) {
    gamma_start <- sample( c(TRUE, FALSE), n_var - 1, replace = TRUE, prob = c(0.5, 0.5) )
  }
  
  ### Run the MCMC routine
  routine_fun <- get(paste0("MIDAS_MCMC_", response_dist))
  routine_fun( formula = form, start = list( beta = beta_start, gamma = gamma_start ), 
               prior = prior, quantile = quantile, MCMC_length = MCMC_length, ... )
  
}

