##### Creates a matrix of time discrepencies from the time window
##### Daniel Dempsey

time_delta_matrix <- function(series, response_index, time_window) {
  
  # Calculate differences between the response indices and the covariate indices
  xt <- series$time
  rownum <- length(response_index)
  colnum <- length(xt)
  yt <- matrix( rep(response_index, each = colnum), nrow = rownum, byrow = TRUE )
  diff_mat <- sweep(yt, 2, xt)
  
  # Determine which values are inside the time window, store results as a list
  within <- ( 0 <= diff_mat ) & ( diff_mat < time_window )
  inds <- which(within, arr.ind = TRUE)
  list( nrow = rownum, ncol = colnum, i = inds[, 1], j = inds[, 2], 
        x = split( diff_mat[within], inds[, 1] ) )
  
}

time_delta_list <- function(regime_object) {
  
  regime_object$time_delta <- Map( time_delta_matrix, 
                                   series = regime_object$series, 
                                   response_index = regime_object$response_index,
                                   MoreArgs = list(time_window = regime_object$time_window) )
  regime_object
  
}
