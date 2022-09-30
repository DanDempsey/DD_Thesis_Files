##### Computes the gradient of the log-likelihood
##### Daniel Dempsey

midas_gradient <- function( pars, pinds, formula, family ) {
  
  ### Extract formula environment
  Zenv <- environment(formula)
  
  ### Set the parameters
  par_list <- split(pars, pinds)
  
  ### Compute necessary parameters
  WX <- model.matrix( formula, Zenv )
  eta <- WX %*% par_list[[1]]
  p <- family$linkinv( eta )
  yv <- Zenv$response_vector
  
  ### Compute common components
  dldp <- ( yv / p ) - ( ( 1 - yv ) / ( 1 - p ))
  dpde <- ( ( 1 + exp( -eta ) )^(-2) ) * exp( -eta )
  dlde <- dldp * dpde
  
  ### Beta components
  dldb <- crossprod( dlde, WX )
  
  ### Theta components
  WX_grad <- Map( midas_design_matrices, M = 1:2, 
                  MoreArgs = list(regime_object = Zenv$regime_object, gr = TRUE) )
  dlde_WX_grad <- do.call( 'rbind', lapply( WX_grad, crossprod, x = dlde ) )
  dldt_mat <- sweep( dlde_WX_grad, 2, par_list[[ 1 ]][-1], '*' )
  dldt_list <- split( dldt_mat, rep( 1:ncol( dldt_mat ), each = 2 ) )
  
  ### Return result
  # Minus two since objective is -2 * log likelihood
  -2 * c( dldb, unlist( dldt_list ) )
  
}
