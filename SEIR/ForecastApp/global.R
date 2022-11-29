####------------------ Covid Predictor: Global Shiny File ------------------####

####-- 1. Required Libraries -----------------------------------------------####
  
  ### List of all required packages for the application to run
  library('shiny')
  library('shinydashboard') 
  library('Matrix')
  library('deSolve')
  library('plotly')
  library('tidyverse')
  library('doParallel')
  library('shinycssloaders')

####-- 2. Required Functions -----------------------------------------------####
  
  ### Locate all files in the code directory
  func_files <- list.files(path = 'code', full.names = TRUE)
  
  ### Load all function files in the code directory
  sapply(func_files, source)
  
####-- 3. Required Data ----------------------------------------------------####
  
  ### Load Data Objects as a List
  data_sources <- load_data(data_path = 'data/', county = 'Dublin')
  
  ### Create Data Frames from the Data Object List
  list2env(data_sources, envir = .GlobalEnv)
  
  ### Remove the List containing the Data Objects
  rm(data_sources)
  
####-- 4. Extra App Data  --------------------------------------------------####

  n_cores <- max( detectCores() - 1, 1 )
  
  ### Create Default Age Groups
  def_age_groups <- c(paste0(seq(0, 70, 5), " - ", seq(4, 74, 5)), "75+")
  forecast_age_groups <- c("0 - 14", paste0(seq(15, 65, 10), " - ", seq(24, 74, 10)), "75+") 
  dub_population$agegroup <- def_age_groups
  
  ### Vector of compartment names
  comp_vec <- c("S", "Ev", "Ip", "IA", "Ii", "It", "Iti", "Iq", "R")
  
  ### Vector of lockdown measures
  lockdown_measures <- c('No Intervention', paste0('Lockdown Level ', 1:5))
  
  ### Today's date
  td <- as.Date(Sys.time())
  
  ### Starting values for ODE solvers
  dub_xstart <- c(dub_population$popage - 15/16,
                  rep(1, 16),
                  rep(1/16, 16),
                  rep(0, 16 * 6))
  names(dub_xstart) <- c(paste0('S_',1:16),
                         paste0('Ev_',1:16),
                         paste0('Ip_',1:16),
                         paste0('IA_',1:16),
                         paste0('Ii_',1:16),
                         paste0('It_',1:16),
                         paste0('Iti_',1:16),
                         paste0('Iq_',1:16),
                         paste0('R_',1:16))

  ### Default SEIR model parameters
  def_pars <- c(3.7, 5.8, 11.6, 0.55, 0.2, 0.8, 0.1, 0.05, 7)
  names(def_pars) <- c('L', 'Cv', 'Dv', 'h', 'f', 'tv', 'q', 'k', 'TT')
  
  ### Starting transmission rate
  dub_def_beta <- getbeta(3.7, pars = def_pars, p_age = dub_population$propage,
            CONTACTMATRIX = contacts)
  
  ### 'Known' and full dataset sizes for forecast tab 
  N_known <- 13
  N_full <- 68
  
  ### Average Multiple Lockdowns
  L0 <- mean(c(optim_res$optim_res[optim_res$policy == 'No Intervention']))
  L2 <- optim_res$optim_res[optim_res$policy == 'Lockdown Level 2']
  L3 <- optim_res$optim_res[optim_res$policy == 'Lockdown Level 3']
  L5 <- optim_res$optim_res[optim_res$policy == 'Lockdown Level 5']
  
  ### Linearly interpolate missing lockdown measures
  L1 <- mean( c( L0, L2 ) )
  L4 <- mean( c( L3, L5 ) )
  L_dat <- data.frame(policy = paste0('Lockdown Level ', c( 1, 4 )), 
                      optim_res = c(L1, L4))
  
  optim_res <- optim_res[ c( 1, 2 ) ]
  optim_res <- rbind( optim_res, L_dat )
  rm( L0, L1, L2, L3, L4, L5, L_dat )
  
  ### Include new interventions into the bootstrap dataset
  boot_lockdown_scalars['Lockdown Level 0'] <- boot_lockdown_scalars$`No Intervention`
  boot_lockdown_scalars['Lockdown Level 2'] <- boot_lockdown_scalars$`Lockdown Level 2`
  boot_lockdown_scalars['Lockdown Level 3'] <- boot_lockdown_scalars$`Lockdown Level 3`
  boot_lockdown_scalars['Lockdown Level 5'] <- boot_lockdown_scalars$`Lockdown Level 5`
  boot_lockdown_scalars['Lockdown Level 1'] <- (boot_lockdown_scalars$`Lockdown Level 0` + boot_lockdown_scalars$`Lockdown Level 2`)/2
  boot_lockdown_scalars['Lockdown Level 4'] <- (boot_lockdown_scalars$`Lockdown Level 3` + boot_lockdown_scalars$`Lockdown Level 5`)/2
  
  ### Use the transposed bootstrapped data
  boot_scales_t <- cbind( policy = names(boot_lockdown_scalars), 
                          as.data.frame(t(boot_lockdown_scalars)) )
  rm(boot_lockdown_scalars)
  
