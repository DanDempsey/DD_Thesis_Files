##### Compute posterior predictive distribution for ALD_vs
##### Daniel Dempsey

##### Function that builds the MIDAS design matrix
run_design_mat <- function(regime_object, offset_only, vars_ind) {
  
  if ( !any(vars_ind) ) 
    return( offset_only )
  
  cbind( offset_only, midas_design_matrices(regime_object[vars_ind] ) )
  
}
  
### Function that creates design matrices for test set
create_test_X <- function(reg_obj, var_sel, DLF, n_y) {
  
  int <- matrix(1, nrow = n_y, ncol = 1)
  N <- ncol(var_sel)
  res <- vector('list', N) 
  for ( i in 1:N ) {
    reg_obj_update <- Map( function(x, y, i) { x$DLF_parameters <- y[i, ] ; 
                           return(x) }, 
                           x = reg_obj, y = DLF, MoreArgs = list(i = i) )
    res[[i]] <- run_design_mat(reg_obj_update, int, var_sel[, i])
  }
  
  res
  
}

### Posterior Predictive Distribution for ALD Variable Selection 
ppd_ALD_vs <- function(MCMC_res, new_form, burn_in = 1000, thin = 10, ...) {
 
  ### Extract list of parameters from model fit
  MCMC_length <- nrow( MCMC_res$betares )
  MCMC_extract <- seq( burn_in + 1, MCMC_length, thin )
  gam <- as.data.frame( t( MCMC_res$vsres[MCMC_extract, ] ) )
  gam_int <- rbind( TRUE, gam )
  bet_full <- as.data.frame( t( MCMC_res$betares[MCMC_extract, ] ) )
  bet <- Map( '[', x = bet_full, y = gam_int )
  DLF <- lapply(MCMC_res$DLFres, function(x, y) { x[y, ] }, y = MCMC_extract )
  
  ### Design matrices
  Zenv <- environment()
  parent.env( Zenv ) <- environment( new_form )
  form <- midas_formula_unpack(new_form, Zenv)
  regime_object <- Zenv$regime_object
  ynew <- Zenv$response_vector
  X_mat_test <- create_test_X(regime_object, gam, DLF, length(ynew)) 
  Xb_test <- do.call('cbind', Map( '%*%', x = X_mat_test, y = bet )) 
  rm( X_mat_test, regime_object )
  
  ### Compute posterior predictive probabilities and ROC
  test_p <- pALD( Xb_test, p = MCMC_res$quantile )
  ppp <- apply( test_p, 1, mean )
  ROC <- pROC::roc( ynew, ppp, ... )
  
  ### Return results
  list( predicted = ppp, actual = ynew, ROC = ROC )
  
}

