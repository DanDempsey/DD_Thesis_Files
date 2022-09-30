##### Function to unpack formula objects containing MIDAS regimes
##### Daniel Dempsey

midas_formula_unpack <- function(formula, Zenv) {
  
  ### Unpack the supplied formula
  stringform <- as.character(formula)
  rhs <- strsplit(stringform, " ")[[3L]]
  lhs <- eval( formula[[2L]], Zenv )
  Zenv$response_vector <- do.call( c, lapply(lhs, "[[", i = "value") )
  
  ### Extract the midas_regime object and remake the formula
  n_terms <- length(rhs)
  regime_indices <- logical( n_terms )
  regime_terms <- rhs
  ind <- 1L
  for ( i in 1:n_terms ) {
    tc <- try(get(rhs[i], Zenv), silent = TRUE)
    if ( inherits(tc, c("midas_regime_list", "midas_regime")) ) {
      len <- length(tc)
      ind <- seq( ind, sum(ind, len, -1L) )
      regime_indices[i] <- TRUE
      rhs[i] <- paste0("midas_design_matrices(regime_object)", collapse = " + ")
      ind <- ind[len] + 1
    }
  }
  if ( all(!regime_indices) )
    stop("No MIDAS regime object supplied.")
  regime_object_function <- paste0( "c(", paste(regime_terms[regime_indices], collapse = ", "), ")" )
  regime_object <- eval( parse(text = regime_object_function), Zenv )  
  
  ### Insert the response times into the regime objects and compute the time window differences
  ro <- lapply( regime_object, "[[<-", i = "response_index", value = lapply(lhs, "[[", i = "time") )
  Zenv$regime_object <- lapply( ro, time_delta_list )
  
  ### Return new formula object
  as.formula( paste0("response_vector ~ ", rhs), env = Zenv )
  
}

  