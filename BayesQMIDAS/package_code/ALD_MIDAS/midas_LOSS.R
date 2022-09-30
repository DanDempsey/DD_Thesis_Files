##### Computes the negative log posterior at given co-ordinates
##### Daniel Dempsey

midas_LOSS <- function( pars, pinds, formula, family ) {
  
  ### Extract formula environment
  Zenv <- environment(formula)
  
  ### Set the parameters parameters
  par_list <- split(pars, pinds)
  Zenv$regime_object <- Map( "[[<-", x = Zenv$regime_object, 
                             i = "DLF_parameters", value = par_list[-1] )
  
  ### Compute the linear predictor
  WX <- model.matrix( formula, Zenv )
  WXb <- WX %*% par_list[[1]]
  eta <- family$linkinv(WXb)
  
  ### Return log posterior density
  sum( family$dev.resids(Zenv$response_vector, eta, 1) )
  
}

