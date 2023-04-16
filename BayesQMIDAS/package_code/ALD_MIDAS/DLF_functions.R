##### Normalised exponential Almon lag weighting model for uneven time series
##### Daniel Dempsey

### Normalised Exponential Almon
#irts_nealmon <- function(delta_t, theta) {
#  
#  p <- length(theta)
#  delta_t_mat <- matrix(rep(delta_t, each = p) ^ (1:p), ncol = p, byrow = TRUE) # Matrix of delta_t polynomials
#  lw <- delta_t_mat %*% theta
#  exp(lw - matrixStats::logSumExp(lw)) # Normalising this way is more numerically stable
#  
#}

### Two degree Exponential Almon
irts_nealmon <- function(delta_t, theta) {
  
  lw <- matrix(rep(delta_t, each = 2) ^ (1:2), ncol = 2, byrow = TRUE) %*% theta
  exp(lw - matrixStats::logSumExp(lw))
  
}

nealmon_gradient_fun <- function( delta_t, theta, M ) {
  
  dat <- matrix( rep( delta_t, each = 2 ) ^ ( 1:2 ), ncol = 2, byrow = TRUE )
  w1 <- exp( dat %*% theta )
  
  dat_trim <- dat[, M ]
  sw <- sum( w1 )
  s1 <- sw * dat_trim
  s2 <- c( crossprod( w1, dat_trim ) )
  
  w2 <- s1 - s2

  w3 <- sw^2 
  
  ( w2 * w1 ) / w3
  
}

nealmon_prior <- function(pars, pinds, beta_sd = 10, theta_sd = 10, theta_shape = 2.5, theta_rate = 1/2.5) {
  
  par <- split(pars, pinds)
  
  # Beta has a normal prior
  beta_prior <- sum( dnorm(par[[1]], mean = 0, sd = beta_sd, log = TRUE) )
  
  # Extract thetas
  theta <- par[-1]
  theta_2_ind <- lengths(theta)
  theta_2 <- unlist( mapply("[", x = theta, i = theta_2_ind, SIMPLIFY = FALSE) )  
  theta_1 <- unlist( mapply("[", x = theta, i = -theta_2_ind, SIMPLIFY = FALSE) )
  
  # Theta 1 has a normal prior
  theta_1_prior <- sum( dnorm(theta_1, mean = 0, sd = theta_sd, log = TRUE) )
  
  # Minus theta 2 has a gamma prior
  theta_2_prior <- sum( dgamma(-theta_2, shape = theta_shape, rate = theta_rate, log = TRUE) )
  
  # Full prior
  sum( beta_prior, theta_1_prior, theta_2_prior )
  
}

### Beta
irts_beta <- function(delta_t, theta) {
  
  ### First apply safety
  eps <- .Machine$double.eps
  x <- delta_t / max(delta_t)
  x <- ifelse(x == 1, x - eps, x)
  x <- ifelse(x == 0, x + eps, x)
  
  ### Compute weights
  w <- dbeta( x, theta[1], theta[2] )
  w / sum( w )
  
}

beta_prior <- function(pars, pinds, slope_sd = 10, theta_shape = 2.5, theta_rate = 1/2.5) {
  
  par <- split(pars, pinds)
  
  # Slope beta have a normal prior
  slope_prior <- sum( dnorm(par[[1]], mean = 0, sd = slope_sd, log = TRUE) )
  
  # DLF thetas have a gamma prior
  theta_prior <- do.call( "sum", mapply(dgamma, x = par[-1], shape = theta_shape, 
                                        rate = theta_rate, MoreArgs = list(log = TRUE)) )
  # Full prior
  sum( slope_prior, theta_prior )
  
}

### Gompertz
irts_gompertz <- function(delta_t, theta) {
  
  s <- delta_t <- theta[3]
  z <- exp(theta[1] * delta_t)
  res <- z * exp(-theta[2] * z)
  res / sum(res)
  
}

irts_reverse_gompertz <- function(delta_t, theta) {
  
  s <- theta[3] - delta_t
  z <- exp(theta[1] * s)
  res <- z * exp(-theta[2] * z)
  res / sum(res)
  
}

### Log Cauchy
irts_logCauchy_DLF <- function(delta_t, theta) {
  
  res <- (1 / delta_t) * 1/(theta[2]^2 + (log(delta_t) - theta[1])^2)
  res / sum(res)
      
}

### Nakagami
irts_nakagami_DLF <- function(delta_t, theta) {
  
  res <- ( delta_t ^ (2*theta[1] - 1) ) * exp( -(theta[1]/theta[2]) * delta_t^2 )
  res / sum(res)
  
}
