##### A function that assigns proper names to the result of the MAP function
##### Daniel Dempsey

giveNames <- function(x) {
 
  ### Beta names
  parlist <- split(x$par, x$pinds)
  default_names <- names(parlist[[1]])
  len_names <- length(default_names)
  x_nms <- names(x$regime_object)
  levs <- lapply(x$regime_object, function(y) { levels(y$series[[1]][[2]]) }) 
  levs <- lapply(levs, "[", i = -1)
  beta_names <- unlist( mapply( function(x, y) { paste0(y, x, "_beta") }, x = levs, y = x_nms, SIMPLIFY = FALSE ) )
  extra_names <- default_names[ sequence(len_names - length(beta_names)) ]
  beta_names <- c(extra_names, beta_names)
  
  ### Theta names
  theta_inds <- x$pinds[x$pinds > 1] - 1
  x_list <- split(x_nms[theta_inds], theta_inds)
  theta_names <- unlist( mapply(function(x) { paste0( x, "_theta", seq(1, length(x)) ) }, x = x_list, SIMPLIFY = FALSE) )
  
  ### Return result
  c(beta_names, theta_names)
  
}

