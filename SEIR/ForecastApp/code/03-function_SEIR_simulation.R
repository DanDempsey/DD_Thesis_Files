### The following code is primarily authored by Fatima-Zahra Jaouimaa

## Define the model for lsoda
SEIR_model_D <- function (t_ind, x, parms) {
  
  # Initialise the time-dependent variables
  S <- x[1:16]
  Ev <- x[17:32]
  Ip <- x[33:48]
  IA <- x[49:64]
  Ii <- x[65:80]
  It <- x[81:96]
  Iti <- x[97:112]
  Iq <- x[113:128]
  R <- x[129:144]
  
  # Define model parameter values
  L <- parms[["L"]]
  denom_1 <- parms[["denom_1"]]
  denom_2 <- parms[["denom_2"]]
  Dv <- parms[["Dv"]]
  f <- parms[["f"]]
  k <- parms[["k"]]
  tv <- parms[["tv"]]
  q <- parms[["q"]]
  TT <- parms[["TT"]]
  linfo <- parms[["linfo"]]
  
  C <- parms[["C1"]] * parms[["intervention_scales"]][(t_ind >= linfo[[1]]) & (t_ind < (linfo[[2]] + 1))]
  
  # calculate the number of infections and recoveries between time t_ind and t_ind + dt
  dSdt <- -S*parms[["beta"]]*C%*%(Ip + parms[["h"]]*IA + It + Iq + k*Ii + k*Iti)/parms[["N_age"]]
  dEvdt <- -Ev/L - dSdt
  dIpdt <- -Ip/denom_1 + (1 - f)*Ev/L
  dIAdt <- -IA/Dv + f*Ev/L
  dIidt <- -Ii/denom_2 + q*Ip/denom_1
  dItdt <- -It/TT + tv*Ip/denom_1
  dItidt <- -Iti/(denom_2 - TT) + It/TT
  dIqdt <- -Iq/denom_2 + (1 - q - tv)*Ip/denom_1
  dRdt <- IA/Dv + Ii/denom_2 + Iq/denom_2 + Iti/(denom_2 - TT)
  
  list(c(dSdt, dEvdt, dIpdt, dIAdt, dIidt, dItdt, dItidt, dIqdt, dRdt))
  
}

SEIR_model_simulation <- function(pars,
                                  dateStart = as.Date('2020-02-28'),
                                  lockdown_information = NULL,
                                  POP = population,
                                  contacts_ireland = contacts,
                                  beta = 0.1816126,
                                  startval,
                                  dt = 1,  
                                  tmax = 225)   
{
  
  ## Load population information
  N_age <- POP$popage
  
  ## Initialising the compartments
  groups <- dim(contacts_ireland[[1]])[2]
  
  ## Setting time scale                     
  numSteps <- tmax/dt
  times <- seq(from = 0, to = tmax, by = dt)
  dateEnd <- dateStart + (tmax - 1)
  
  ## Defining time points at which interventions come in
  linfo <- data.frame(c1 = difftime(lockdown_information[[1]], dateStart, units = "days"),
                      c2 = difftime(lockdown_information[[2]], dateStart, units = "days"))
  
  ## defining all parameters required for solving model equations
  L <- pars[1]
  Cv <- pars[2]
  Dv <- pars[3]
  parms <- list(L = L, denom_1 = Cv - L, Dv = Dv, h = pars[4],
                denom_2 = Dv - Cv + L, f = pars[5], tv = pars[6], 
                q = pars[7], k = pars[8], TT = pars[9], beta = beta, 
                N_age = N_age, C1 = contacts_ireland[[5]], linfo = linfo, 
                intervention_scales = lockdown_information[[3]])
  
  ## Solving the equations and returning result
  sol <- lsoda(startval, times, SEIR_model_D, parms)
  
  list( solution = as.data.frame(sol), N_age = N_age, beta = beta,
        dateStart = dateStart, dateEnd = dateEnd,
        lockdown_information = lockdown_information )
  
}

