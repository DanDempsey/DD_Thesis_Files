### Use R0 to get beta
getbeta <- function(R0t, pars, p_age, CONTACTMATRIX = contacts) {
  
  ### Extract Parameters
  h <- pars["h"]
  k <- pars["k"]
  L <- pars["L"]
  Cv <- pars["Cv"]
  Dv <- pars["Dv"]
  f <- pars["f"]
  q <- pars["q"]
  tv <- pars["tv"]
  TT <- pars["TT"]
  
  ### Prepare population contact matrix
  n <- dim(CONTACTMATRIX[[1]])[1]
  Csym <- lapply(CONTACTMATRIX, function(x, p_age) (x + t(x)*((p_age)%*%t(1/p_age)))/2, p_age) # make sure contacts are reciprocal
  C <- Csym[[1]] + Csym[[2]] + Csym[[3]] + Csym[[4]]
  
  ### Create the N matrix
  transmission_rates <- c(0, 1, h, k, 1, k, 1)
  N_vals <- matrix(0, nrow = n, ncol = n*7) # Initialise matrix of non-zero values
  col_ind <- seq(1, n*7, 7)
  for(i in 1:n) {
    count <- 1
    for(j in 1:n) {
      ci <- col_ind[count]
      N_vals[i,ci:(ci+6)] = transmission_rates * C[i,j] * p_age[i]/p_age[j] # See report for derivation
      count <- count + 1
    }
  }
  ### N is the F matrix with beta factored out
  N <- sparseMatrix( dims = rep(n*7, 2), i = rep(seq(1, 7*n, 7), n*7), j = rep(seq(1, n*7), each = n), x = as.vector(N_vals) )
  
  ### Create the inverse V matrix
  V_inv_vals <- c(L, (1-f)*(Cv-L), Cv-L, f*Dv, Dv, (1-f)*q*(Dv-Cv+L),
                  q*(Dv-Cv+L), Dv-Cv+L, TT*tv*(1-f), TT*tv, TT,
                  tv*(1-f)*(Dv-Cv+L-TT), tv*(Dv-Cv+L-TT), Dv-Cv+L-TT,
                  Dv-Cv+L-TT, (1-q-tv)*(1-f)*(Dv-Cv+L), (1-q-tv)*(Dv-Cv+L),
                  Dv-Cv+L)
  V_inv_i <- c(1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7)
  V_inv_j <- c(1, 1, 2, 1, 3, 1, 2, 4, 1, 2, 5, 1, 2, 5, 6, 1, 2, 7)
  V_inv_list <- rep( list(spMatrix(nrow = 7, ncol = 7, i = V_inv_i, j = V_inv_j, x = V_inv_vals)), n)
  
  V_inv <- bdiag(V_inv_list) # See report for why it's block diagonal
  
  ### Calculate and return beta
  eig <- eigen(N%*%V_inv)
  R0t / max(Re(eig$values))
  
}

