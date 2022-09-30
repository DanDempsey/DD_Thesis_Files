###### Default Priors for IRTS MIDAS Parameters
###### Daniel Dempsey

defaultPrior <- function(formula, data, incl_DLF_pars, varsel) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm )
  
  # Beta prior
  prior_beta <- list( beta0 = rep(0, n_beta), 
                      V0 = 100 * diag(n_beta) )
  
  prior <- list(beta = prior_beta)
  
  # Variable inclusion prior
  if ( varsel )
    prior$vars <- c( 1, rep(1 / (n_beta + 1), n_beta - 1) )
  
  # DLF prior
  if ( incl_DLF_pars )
    prior$DLF_pars <- list( DLF1 = c(0, 1), DLF2 = c(0, 1) )
  
  prior
  
}

defaultPrior_ald <- function(formula, data) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm ) 
  list( beta0 = rep(0, n_beta), 
        V0 = 100 * diag(n_beta) )
  
}

defaultPrior_ald_vs <- function(formula, data) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm )
  
  # Beta prior
  prior <- list( beta0 = rep(0, n_beta), V0 = 100 * diag(n_beta) )
  
  # DLF prior
  prior$DLF_pars <- list( DLF1 = c(15, 10), DLF2 = c(1, 1) )
  
  # Variable inclusion prior
  prior$model_selection <- c( 1, 1 )
  
  prior
  
}
