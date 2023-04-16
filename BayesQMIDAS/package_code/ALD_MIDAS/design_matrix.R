##### A function that computes the IRTS-MIDAS design matrix
##### Daniel Dempsey

### Apply DLF
weight_matrix <- function( DLF, DLF_parameters, time_window_list, ... ) {
  
  DLF_weights <- Map( DLF, delta_t = time_window_list$x, 
                      MoreArgs = list(theta = DLF_parameters, ...) )
  Matrix::spMatrix( nrow = time_window_list$nrow, ncol = time_window_list$ncol, 
                    i = time_window_list$i, j = time_window_list$j, 
                    unsplit(DLF_weights, time_window_list$i) )
  
}

### Design matrix for individual regime objects
midas_design_matrix <- function( regime_object, gr = FALSE, ... ) {
  
  # Compute the weight matrix
  weight_mat_list <- Map(weight_matrix, time_window_list = regime_object$time_delta,
                         MoreArgs = list(DLF = ifelse( gr, regime_object$DLF_gradient, 
                                                       regime_object$DLF ), 
                                         DLF_parameters = regime_object$DLF_parameters, ...))
  
  # Compute the WX matrix
  WX_list <- Map( "%*%", x = weight_mat_list, y = regime_object$model_matrix )
  do.call( "rbind", WX_list )
  
}


### Wrapper function to apply across multiple regime objects
midas_design_matrices <- function( regime_object, gr = FALSE, ... ) {

  design_list <- Map( midas_design_matrix, regime_object = regime_object, 
                      MoreArgs = list( gr = gr, ... ) )
  as.matrix( do.call("cbind", design_list) )
  
}

