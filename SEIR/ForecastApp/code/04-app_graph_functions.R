##### Extra functions for app, mainly for graphics

### A convenience function simply to extract appropriate compartment
comp_extract <- function(dat, comp) {
  dat[grepl(paste0(comp, "_"), names(dat))]
}

### Summary plot function for model dashboard
summary_plot <- function(dat, xval, comp_vec, group, I_type){
  
  # Extract Compartments
  comps <- Map(comp_extract, comp = comp_vec,
               MoreArgs = list(dat = dat))
  names(comps)[3:8] <- c("I_Pr", "I_As", "I_Im", "I_Aw", "I_Is", "I_No")
  comps$I_Sy <- comps$I_Im + comps$I_Aw + comps$I_Is + comps$I_No
  comps$I_Al <- comps$I_Pr + comps$I_As + comps$I_Sy 
  comps$I <- comps[[paste0("I_", substr(I_type, start = 1, stop = 2))]]
  
  # Sum across desired age groups
  comp_groups <- lapply(comps[c("S", "Ev", "I", "R")],
                        "[", i = group)
  comp_draw <- lapply(comp_groups, rowSums)
  
  # Draw the Plot
  plot_ly(x = ~xval, y = ~comp_draw$S, name = 'Susceptible', type = 'scatter', mode = 'lines',
          line = list(color = 'rgb(69, 95, 245)')) %>%
    add_trace(y = ~comp_draw$Ev, name = 'Exposed', mode = 'lines', line = list(color = 'rgb(214, 122, 17)')) %>% 
    add_trace(y = ~comp_draw$I, name = 'Infected', mode = 'lines', line = list(color = 'rgb(186, 24, 19)')) %>%
    add_trace(y = ~comp_draw$R, name = 'Removed', mode = 'lines', line = list(color = 'rgb(23, 191, 26)')) %>% 
    layout(xaxis = list(title = "Time"), 
           yaxis = list(title = "Compartment Size"), 
           title = list(text = "Overall Output"))
  
}

### Individual component plots for age breakdowns in model dashboard
comp_plot <- function(dat, xval, comp, group, inp){
  
  comp_dat <- dat[grepl(comp, names(dat))][group]
  
  # Draw the Plot
  p <- plot_ly(x = ~xval, y = NA, type = 'scatter', mode = 'lines')
  
  for ( i in seq_along(inp$age_sel) ) {
    p <- p %>% add_trace(y = comp_dat[[i]], name = inp$age_sel[i], mode = 'lines')
  }
  
  p 
  
}

### A function that cleans up compartments, gives them appropriate names,
### and adds a 'Symptomatic' and 'All' infectious compartments
comp_surmise <- function(x) {
  comps <- Map(comp_extract, comp = comp_vec, 
               MoreArgs = list(dat = x))
  names(comps) <- c("Su", "Ex", "Pr", "As", "I_Im", 
                    "I_Aw", "I_Is", "I_No", "Re")
  comps$Sy <- comps$I_Im + comps$I_Aw + comps$I_Is + comps$I_No
  comps$Al <- comps$Pr + comps$As + comps$Sy
  comps
}

### Create forecasting data
comp_sel <- function(x, y, z = 1:16) {
  comps <- comp_surmise(x)
  comps_sel <- comps[[substr(y, start = 1, stop = 2)]]
  N <- nrow(comps_sel)
  if (length(z) == 1) { return(comps_sel[(N-N_full+1):N, z]) }
  rowSums(comps_sel[(N-N_full+1):N, z])
}

### Makes it so that any date without an explicit intervention is
### assigned "No intervention"
intervention_adjust <- function(dat, x) {
  
  old_end <- dat$end[nrow(dat)]
  missing_interventions <- as.numeric(difftime(x, old_end, units = "days"))
  if (missing_interventions <= 1)
    return(dat)
  rbind(dat, data.frame(start = old_end + 1, 
                        end = x, 
                        policy = "No Intervention"))
  
}
  
### Extract quantiles for forecast tab to create upper and lower bounds 
df_quant <- function(x, p) {
  sapply(x, quantile, probs = p)
}

### Compute estimated deaths
comp_deaths <- function(x, forecaster = FALSE) {
  comps <- comp_surmise(x)
  denom_1 <- def_pars["Dv"] - def_pars["Cv"] + def_pars["L"]
  denom_2 <- denom_1 - def_pars["TT"]
  comps_sel <- (comps$I_Im + comps$I_No)/denom_1 + comps$I_Is/denom_2
  N <- nrow(comps_sel)
  forcast <- (N-N_full+N_known+1):N 
  
  if (!forecaster) {
    ags <- rep(0, 8)
    ags[1] <- sum(comps_sel[forcast, 1:3])
    ags[2] <- sum(comps_sel[forcast, 4:5])
    ags[3] <- sum(comps_sel[forcast, 6:7])
    ags[4] <- sum(comps_sel[forcast, 8:9])
    ags[5] <- sum(comps_sel[forcast, 10:11])
    ags[6] <- sum(comps_sel[forcast, 12:13])
    ags[7] <- sum(comps_sel[forcast, 14:15])
    ags[8] <- sum(comps_sel[forcast, 16])
    return(ags * est_deaths$Estimated_Deaths)
  }
  ags <- rep(0, 9)
  ags[1] <- sum(comps_sel[forcast, 1:3])
  ags[2] <- sum(comps_sel[forcast, 4:5])
  ags[3] <- sum(comps_sel[forcast, 6:7])
  ags[4] <- sum(comps_sel[forcast, 8:9])
  ags[5] <- sum(comps_sel[forcast, 10:11])
  ags[6] <- sum(comps_sel[forcast, 12:13])
  ags[7] <- sum(comps_sel[forcast, 14])
  ags[8] <- sum(comps_sel[forcast, 15])
  ags[9] <- sum(comps_sel[forcast, 16])
  est_alt <- est_deaths$Estimated_Deaths[c(1:8, 8)]
  prop_70 <- as.numeric( dub_population[15, 2] / sum(dub_population[14:15, 2]) )
  est_alt[7] <- est_alt[7] * (1 - prop_70)
  est_alt[8] <- est_alt[8] * prop_70
  ags * est_alt
}



