##### Function that runs the optimization function to compute MAP
##### Daniel Dempsey

irts_midas.fit <- function(x) {
  
  ### Run the optimisation function if user did not request a dry run
  Ofunction <- x$Ofunction
  x[ which( x == "optim" ) ] <- NULL
  if (Ofunction != "dry_run") {
    #opt <- optim(par = x$par, fn = x$fn, family = x$family, prior = x$prior, lp = x$lp,
    #             pinds = x$pinds, formula = x$formula, method = x$method, hessian = x$hessian)
    opt <- try(do.call(Ofunction, x), silent = TRUE)
    if (inherits(opt, "try-error")) {
      stop("The optimisation algorithm failed with the following message:\n\n", 
           opt, "\nPlease try other starting values or a different optimisation function.")
    }
    
    ### Extract reuslts
    if (Ofunction == "optimx") {
      bmet <- which.min(opt$value)
      par <- as.numeric(opt[bmet, 1:length(args$par)])
    }
    else { par <- opt$par }
    x$convergence <- opt$convergence
  }
  
  ### Do nothing if user requested a dry run
  else { x$opt <- list(message = "Dry run: no optimisation performed", par = par) }
  
  ### Give informative names to the parameters and return the result
  x$regime_object <- evalq(regime_object, environment(x$formula))
  x$opt <- opt
  x$par <- par
  nms <- giveNames(x)
  names(x$par) <- names(x$opt$par) <- nms
  if (!is.null(x$opt$hessian)) {
    colnames(x$opt$hessian) <- rownames(x$opt$hessian) <- nms
  }
  
  ### Compute model diagnostics from the glm function
  #form <- x$formula
  #WX <- model.matrix( form, environment(form) )
  #colnames(WX) <- split(nms, x$pinds)[[1]]
  #x$glm_diagnostics <- dry_glm_fit(WX, environment(form)$response_vector, 
  #                                 start = split(x$par, x$pinds)[[1]], family = x$family)
  x
  
}

