##### Sampler for a truncated ALD
##### SOURCE PACKAGE

pALD <- function(q, mu = 0, sigma = 1, p = 0.5) {
  ifelse( test = q < mu, 
          yes = p * exp((1 - p) * (q - mu)/sigma), 
          no = 1 - (1 - p) * exp(-p * (q - mu)/sigma) )
}

qALD <- function (prob, mu = 0, sigma = 1, p = 0.5)  {
  ifelse( test = prob < p, 
          yes = mu + (sigma*log(prob/p))/(1 - p), 
          no = mu - sigma*log((1 - prob)/(1 - p))/p )
}

# Used elsewhere
rALD <- function(n, mu = 0, sigma = 1, p = 0.5) {
  u <- runif(n)
  mapply(qALD, prob = u, mu = mu, sigma = sigma, p = p)
}

# Not used anywhere, but included here for the sake of completeness (based on ald package)
dALD <- function (y, mu = 0, sigma = 1, p = 0.5) {
  ifelse(test = y < mu, 
         yes = (p * (1 - p)/sigma) * exp((1 - p) * (y - mu)/sigma), 
         no = (p * (1 - p)/sigma) * exp(-p * (y - mu)/sigma))
}

rTALD <- function(n, rtrunc, mu = 0, sigma = 1, p = 0.5) {
  bound <- pALD( 0, mu = mu, sigma = sigma, p = p )
  u <- runif( n, min = ifelse( rtrunc, bound, 0 ), max = ifelse( rtrunc, 1, bound ) )
  qALD( prob = u, mu = mu, sigma = sigma, p = p )
}

# Benoit 2017 version
truncALD <- function(y, mu, sigma, p) {
  trunc <- -mu/sigma
  if ((y == 0) & (trunc <= 0)) {
    u <- runif(n = 1)
    fn_val <- -(abs(trunc) - log(u)) / (1 - p)
  } else if ((y == 1) & (trunc >= 0)) {
    u <- runif(n = 1)
    fn_val <- trunc - log(u) / p
  } else {
    g <- pALD(q = trunc, mu = 0, sigma = 1, p = p)
    if (y == 1) {
      min <- g
      max <- 1
    } else {
      min = 0
      max = g
    }
    u <- runif(n = 1)
    u <- min + (max - min) * u
    fn_val <- qALD(prob = u, mu = 0, sigma = 1, p = p)
  }
  fn_val <- mu + fn_val * sigma
}
