##### Midas regime characteristics for a collection of covariates
##### Daniel Dempsey

mr <- function(value, time, group, time_window = 30, DLF = "irts_nealmon", 
               DLF_gradient = "nealmon_gradient_fun", DLF_parameters) {
  
  ### Checking input
  # value cannot be missing
  if ( missing(value) )
    stop("value must be specified.")
  # If time is missing, assume regularity
  if ( missing(time) )
    time <- seq(1, NROW(value))
  # value and time must be the same length
  if (NROW(value) != length(time))
    stop("value and time must be the same length.")
  # group must be NULL or the same length as value
  if (!missing(group)) {
    if (NROW(value) != length(group))
      stop("if not missing, group must be the same length as value.")
  }
  # The DLF must be either a function or the name of a function
  if (is.character(DLF))
    DLF <- try(get(DLF, mode = "function", envir = parent.frame()), silent = TRUE)
  if (!is.function(DLF))
    stop("DLF must be a function or a string that corresponds to the name of a function.")
  # Same as above for DLF gradient function
  if (is.character(DLF_gradient))
    DLF_gradient <- try(get(DLF_gradient, mode = "function", envir = parent.frame()), silent = TRUE)
  if (!is.function(DLF_gradient))
    stop("The DLF gradient must be a function or a string that corresponds to the name of a function.")
  # The time window must be greater than 0
  if (time_window <= 0)
    stop("time_window must be greater than 0.")
  # Check that the DLF and number of parameters are valid
  DLF_check <- try(DLF(1, DLF_parameters), silent = TRUE)
  if (inherits(DLF_check, "try-error"))
    stop(paste0("\nThe DLF or starting weights are not valid. \nThis may have occurred if the DLF you provided has no DLF_parameter argument.\n\nThe error message reads as follows:\n", DLF_check[1]))
  
  ### Compile all the characteristics into a list
  if ( missing(group) )
    group <- rep("All", nrow(series))
  series <- data.frame(time = time, value = value)
  ord <- order(series$time)
  series <- split(series[ord, ], group[ord])
  series <- lapply(series, na.omit)
  model_matrix_int <- Map(model.matrix, data = series, MoreArgs = list(~value))
  model_matrix <- lapply(model_matrix_int, subset.matrix, select = -1)
  res <- mget(c("series", "time_window", "DLF", "DLF_gradient", "DLF_parameters", "model_matrix"))
  class(res) <- c("midas_regime")
  res
  
}

midas_regimes <- function(value, time, group, time_window = 30, DLF = "irts_nealmon", 
                          DLF_gradient = "nealmon_gradient_fun", DLF_parameters) {
  
  # Value cannot be missing
  if (missing(value))
    stop("value must be specified.")
  # If time is missing then set it as a regular ts 
  if (missing(time))
    time <- mapply(seq_along, value)
  # Change DLF to a list
  if (is.function(DLF))
    DLF <- list(DLF)
  # Change DLF gradient to a list
  if (is.function(DLF_gradient))
    DLF_gradient <- list(DLF_gradient)
  # Change group
  if (missing(group))
    group <- rep("All", NROW(value))
  # Formatting
  format_fix <- function(x) {
    if (is.matrix(x))
      x <- as.data.frame(x)
    if (!is.list(x))
      x <- list(x)
    x
  }
  value <- format_fix(value)
  time <- format_fix(time)
  group <- format_fix(group)
  DLF_parameters <- format_fix(DLF_parameters)
  
  ### Run the mr function across all values
  res <- Map(mr, value = value, time = time, group = group, time_window = time_window, 
             DLF = DLF, DLF_gradient = DLF_gradient, DLF_parameters = DLF_parameters)
  
  # Make sure every element has a name
  n <- names(res)
  if (is.null(n))
    names(res) <- paste0("midas_covariate_", seq(1, length(res)))
  m <- which(n == "")
  names(res)[m] <- paste0("midas_covariate_", seq(1, length(m)))
  
  class(res) <- "midas_regime_list"
  res
  
}

response_regime <- function(value, time, group) {
  
  ### Sense Checks
  # value cannot be missing
  if ( missing(value) )
    stop("value must be specified.")
  # If time is missing, assume regularity
  if ( missing(time) )
    time <- seq(1, NROW(value))
  # value and time must be the same length
  if (NROW(value) != length(time))
    stop("value and time must be the same length.")
  # group must be NULL or the same length as value
  if (!missing(group)) {
    if (NROW(value) != length(group))
      stop("if not missing, group must be the same length as value.")
  }
  else
    group <- rep("All", NROW(value))
  
  ### Compile all the characteristics into a list
  series <- data.frame(time = time, value = value)
  ord <- order(series$time)
  series <- split(series[ord, ], group[ord])
  series <- lapply(series, na.omit)
  class(series) <- c("response_regime")
  series
  
}

