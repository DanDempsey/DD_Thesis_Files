##### Main function that performs inference via MCMC
MIDAS_MCMC_ald_vs_DLF_fixed <- function(formula, start, prior, quantile, MCMC_length) {
  
  ### Prepare data environment
  Zenv <- environment( formula )
  regime_object <- Zenv$regime_object
  nvar <- length(regime_object)
  
  ### MCMC setup (first iteration)
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  delta <- 2 + ( ( psi^2 ) / omega )
  V0i_full <- chol2inv( chol(prior$V0) )
  V0i <- V0i_full[1, 1]
  V0ib0_full <- V0i_full %*% prior$beta0
  V0ib0 <- V0ib0_full[1, 1]
  
  int <- X_mat <- matrix(1, nrow = n_y, ncol = 1)
  varsel <- rep( FALSE, nvar )
  nu <- rep( 1, n_y )
  Ni <- 1 / ( omega * nu )
  z <- rep( 0, n_y )
  z_pn <- z - ( psi * nu ) 
  
  # Beta initialization
  betares <- matrix( 0, ncol = nvar+1, nrow = MCMC_length )
  vsres <- matrix( FALSE, ncol = nvar, nrow = MCMC_length )
  mod_sel_prior <- prior$model_selection
  DLF_prior_mean1 <- prior$DLF_pars[[1]][1]
  DLF_prior_mean2 <- log( prior$DLF_pars[[2]][1] / prior$DLF_pars[[2]][2] )
  DLF_prior_mean <- c( DLF_prior_mean1, DLF_prior_mean2 )
  model_sel_prob <- rep( mod_sel_prior[1] / sum( mod_sel_prior ), length = MCMC_length )
  varnames <- names( Zenv$regime_object )
  colnames(betares) <- c( "Int", varnames )
  colnames(vsres) <- varnames
  betares[1, ] <- start
  varsel <- vsres[1, ]
  varsel_int <- c( TRUE, varsel )
  
  # DLF parameter initialization
  var_extract <- paste0('regime_object$', varnames, '$DLF_parameters')
  DLFres_cmd <- Map( parse, text = var_extract )
  DLF_1 <- lapply( DLFres_cmd, eval, envir = Zenv )
  names(DLF_1) <- varnames
  DLFres <- lapply( DLF_1, function(x) { matrix(rep(x, each = MCMC_length), nrow = MCMC_length) } )
  
  var_inds <- 1:nvar
  birth_oppurtunity <- birth <- death_oppurtunity <- death <- accept <- rep( 0, nvar )
  jump_proposal_DLF_sigma <- diag( c(5, 2) )
  covar_sd <- ( 2.38^2 )/2
  small_pos <- diag( 0.000001, 2 )
  DLF_last_covar <- replicate( nvar, matrix(c(6, 0.5, 0.5, 0.05), ncol = 2), simplify = FALSE )
  DLF_last_mean <- replicate( nvar, DLF_prior_mean, simplify = FALSE )
  names( accept ) <- names( birth_oppurtunity ) <- 
    names( birth ) <- names( death_oppurtunity ) <- 
    names( death ) <- varnames 
  rm( var_extract, DLFres_cmd, DLF_1, varnames, DLF_prior_mean1, DLF_prior_mean2 )
  
  ### Fix DLF values (for bugfixing)
  regime_object <- lapply( regime_object, function(x) { x$DLF_parameters <- c(14, -1); x } )
  regime_object <- lapply( regime_object, function(x) { x$WX <- midas_design_matrix( x ); x } )
  
  ### Main Loop
  for ( i in 2:MCMC_length ) {
    
    ### Update covariate indicator
    V0i <- V0i_full[varsel_int, varsel_int]
    V0ib0 <- V0ib0_full[varsel_int]
    XtNi <- t( X_mat * Ni )
    
    V_posti <- V0i + XtNi%*%X_mat
    V_post <- chol2inv( chol(V_posti) )
    B_post <- V_post%*%( V0ib0 + XtNi%*%z_pn )
    
    if ( i%%10 != 0 ) { 
      
      vsres[i, ] <- vsres[i-1, ]
      model_sel_prob[i] <- model_sel_prob[i-1]
      
    } else {
      
      ### Propose dimension change
      change_ind <- sample( var_inds, 1 )
      varsel_star <- varsel
      varsel_star[change_ind] <- !varsel_star[change_ind]
      varsel_star_int <- c( TRUE, varsel_star )
      change_name <- names( regime_object )[change_ind]
      
      if( varsel[change_ind] ) {
        death_oppurtunity[change_name] <- death_oppurtunity[change_name] + 1
      } else {
        birth_oppurtunity[change_name] <- birth_oppurtunity[change_name] + 1
      }
      
      # Holmes and Held (2006) ratio
      X_mat_star <- make_design_mat( regime_object, int, varsel_star )
      V0i_star <- V0i_full[varsel_star_int, varsel_star_int]
      V0ib0_star <- V0ib0_full[varsel_star_int]
      XtNi_star <- t( X_mat_star * Ni )
      
      V_posti_star <- V0i_star + XtNi_star%*%X_mat_star
      V_post_star <- chol2inv( chol( V_posti_star ) ) 
      B_post_star <- V_post_star%*%( V0ib0_star + XtNi_star%*%z_pn )
      
      ldet_V_post <- sum( log(diag(chol(V_post))) ) 
      ldet_V_post_star <- sum( log(diag(chol(V_post_star))) )
      ldet_V0i <- sum( log( diag(chol(V0i)) ) )
      ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )
      
      varsel_lprior <- sum( dbinom( varsel, 1, model_sel_prob[i-1], log = TRUE ) )
      varsel_lprior_star <- sum( dbinom( varsel_star, 1, model_sel_prob[i-1], log = TRUE ) )
      
      lkernel <- crossprod( B_post, V_posti )%*%B_post/2
      lkernel_star <- crossprod( B_post_star, V_posti_star )%*%B_post_star/2
      
      ldenom <- sum( ldet_V_post, ldet_V0i, lkernel, varsel_lprior )
      lnum <- sum( ldet_V_post_star, ldet_V0i_star, lkernel_star, varsel_lprior_star ) 
      
      # Acceptance probability
      if ( (lnum - ldenom) > log(runif(1)) ) {
        if( !varsel[change_ind] ) {
          birth[change_name] <- birth[change_name] + 1
        } else {
          death[change_name] <- death[change_name] + 1
        }
        vsres[i, ] <- varsel <- varsel_star
        varsel_int <- varsel_star_int
        X_mat <- X_mat_star
        V_post <- V_post_star
        B_post <- B_post_star
      } else {
        vsres[i, ] <- vsres[i-1, ]
      }
      
      ### Update the prior success probability p
      model_sel_prob[i] <- rbeta( 1, mod_sel_prior[1] + sum(varsel), 
                                  mod_sel_prior[2] + nvar - sum(varsel) )
      
    }
    
    ### Update ALD parameters
    # Update beta
    bet <- betares[i, varsel_int] <- t( rmvnorm( 1, mean = B_post, sigma = V_post, checkSymmetry = FALSE ) )
    
    # Update z
    Xb <- X_mat %*% bet
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
  
    # Update nu
    chi <- (z - Xb)^2 / omega
    nu <- 1/rinvgauss( n_y, mean = sqrt(delta/chi), shape = delta )
    #for (k in 1:n_y) {
    #  nu[k] <- rgig( n = 1, lambda = 0.5, chi = chi[k], psi = delta )
    #}
      
    ### Update DLF parameters
    z_pn <- z - (psi * nu)
    Ni <- 1 / ( omega * nu )
    ll <- ( z_pn * Ni ) %*% Xb - crossprod(Xb*Ni, Xb)/2
    
    # Print progress
    if ( !i%%500 )
      cat( paste0("Current iteration: ", i, "\n") )
    
  }
  
  birth_acceptances <- birth / birth_oppurtunity
  death_acceptances <- death / death_oppurtunity
  all_acceptances <- ( birth + death ) / ( birth_oppurtunity + death_oppurtunity )
  
  list(betares = betares, DLFres = DLFres, vsres = vsres,
       model_sel_prob = model_sel_prob,
       quantile = quantile, response_vec = Zenv$response_vector, 
       prior = prior, accept = accept,
       birth_oppurtunity = birth_oppurtunity,
       death_oppurtunity = death_oppurtunity,
       birth = birth, death = death, 
       birth_acceptances = birth_acceptances,
       death_acceptances = death_acceptances,
       all_acceptances = all_acceptances,
       regime_object = regime_object)
  
}
