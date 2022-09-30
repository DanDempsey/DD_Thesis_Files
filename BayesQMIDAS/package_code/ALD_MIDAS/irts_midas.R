##### Infers the Maximum A Posteriori of the IRTS-MIDAS
##### Daniel Dempsey

irts_midas <- function( formula, data, start = NULL, family = "binomial", 
                        n_cores = 1, gr = midas_gradient, Ofunction = "optim", ... ) {
  
  ### Initialise environment to be passed into model.frame
  if (missing(data))
    Zenv <- new.env(parent = environment())
  else {
    Zenv <- as.environment(data)
    parent.env(Zenv) <- environment()
  }
  
  ### Check that the family parameter is valid
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    stop("Supplied family is not recognized.")
  }
  
  ### Check that Ofunction parameter is valid
  if ( is.function(Ofunction) )
    Ofunction <- as.character( substitute(Ofunction) )
  if( !is.character(Ofunction) )
    stop("Ofunction must be either a function or string.")
  if( !(Ofunction %in% c("optim", "spg", "optimx", "dry_run")) ) 
    stop("Supplied Ofunction is not in the supported optimization functions list.")
  
  ### Unpack the formula
  form <- midas_formula_unpack( formula, Zenv )
  
  ### Assign starting values and prepare parameters
  MIDAS_start( form, Zenv )
  
  ### Gather all relevant info and run the fitting function
  Ofunction_pars <- list( Ofunction = Ofunction, par = Zenv$pars, fn = midas_LOSS, 
                          family = family, pinds = Zenv$pinds, formula = form, 
                          gr = gr, ... )
  invisible( list2env(list(Ofunction = Ofunction, n_cores = n_cores), envir = Zenv) )
  irts_midas.fit( Ofunction_pars )
  
}

