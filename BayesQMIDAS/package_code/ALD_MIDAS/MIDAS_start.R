##### Function to extract starting values for MAP estimation and MCMC inference
##### Daniel Dempsey

MIDAS_start <- function(formula, Zenv) {
  
  startm <- model.matrix( formula, Zenv )
  if ( is.null(get("start", Zenv)) )
    start <- glm.fit(startm, Zenv$response_vector, family = get("family", Zenv))$coefficients
  n_DLF <- length( Zenv$regime_object )
  n_beta <- ncol(startm)
  if (length(start) > n_beta) {
    warning("Too many start slope parameters supplied. Ignoring extra terms.")
    start <- start[1L:n_beta]
  }
  if (length(start) < n_beta) {
    warning("Not enough start slope parameters supplied. Extra terms will be set to zero.")
    start <- c(start, rep(0, length(start) - n_beta))
  }
  theta_start <- lapply(Zenv$regime_object, "[[", i = "DLF_parameters")
  Zenv$pars <- c( start, unlist(theta_start) )
  Zenv$pinds <- c( rep(1L, n_beta), rep(seq(2L, n_DLF + 1L), lengths(theta_start)) )
  
}

MIDAS_start_ald <- function(formula, start, prior, Zenv) {
  
  n_y <- length(Zenv$yv)
  pinds <- c( rep(1, length(prior$beta)), rep(2, n_y), rep(3, n_y) )
  
  if ( length(start) < length(pinds) ) {
    if ( !is.null(start) )
      warning("Not enough starting values supplied; ignoring supplied vector and using default instead.")
    start <- c( prior$beta, rep(0.01, n_y), rep(1, n_y) )
  }
  
  if ( length(start) > length(pinds) ) {
    warning("Too many starting values supplied; truncating extra terms.")
    start <- start[ seq_along(pinds) ]
  }
  
  Zenv$pinds <- pinds
  start
  
}
