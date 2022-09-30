## load data: population age structure, contact matrices, daily accumulated cases,
##            and interventions (start and end dates)

load_data <- function(data_path = 'data/', county = 'Dublin', 
  load_county_pop_data = TRUE, load_IRL_pop_data = TRUE, 
  load_contact_matrices = TRUE, load_case_data = FALSE, 
  load_interventions = TRUE, load_bootstrapped_scalars = TRUE,
  load_optim_scalars = TRUE, load_deaths = TRUE, load_costs = TRUE) {
  
  # Empty list for output - allows for variable output size
  out <- vector("list", sum(load_county_pop_data, load_IRL_pop_data, 
                            load_contact_matrices, load_case_data, 
                            load_interventions, load_bootstrapped_scalars,
                            load_optim_scalars, load_deaths, load_costs))
  count <- 0
  
  # Load County Population Data
  if (isTRUE(load_county_pop_data)) {
    count <- count + 1
    pop_file <- paste0(data_path, county, '_pop_2019.csv')
    out[[count]] <- readr::read_csv(pop_file, col_types = cols())
    names(out)[[count]] <- 'dub_population'
  }
  
  # Load County Population Data
  if (isTRUE(load_IRL_pop_data)) {
    count <- count + 1
    pop_file <- paste0(data_path, 'Ireland_pop_2019.csv')
    out[[count]] <- readr::read_csv(pop_file, col_types = cols())[-1]
    names(out)[[count]] <- 'irl_population'
  }
  
  # Load (Projected) Contact Matrices Data 
  if (isTRUE(load_contact_matrices)) {
    count <- count + 1
    contacts_file <- paste0(data_path, 'contacts_IRL.Rdata')
    out[[count]] <- get(load(contacts_file))
    names(out)[[count]] <- 'contacts'
  }
  
  # Load Case Data
  if (isTRUE(load_case_data)) {
    count <- count + 1
    cases_file <- paste0(data_path, 'Covid19CountyStatisticsHPSCIreland.csv')
    out[[count]] <- read_csv(cases_file, col_types = cols()) %>% 
      filter(CountyName == county) %>%
      select(TimeStamp, ConfirmedCovidCases) %>% 
      rename(date = TimeStamp , cases = ConfirmedCovidCases) %>%
      mutate(date = as.Date(date)) 
    names(out)[[count]] <- 'cumulative_cases'
  }
  
  # Load Intervention Data
  if (isTRUE(load_interventions)) { 
    count <- count + 1
    interventions_file <- paste0(data_path, county, '_Interventions.csv')
    out[[count]] <- read_csv(
        interventions_file, col_types = cols()
      ) %>% 
      mutate(
        start = as.Date(start, "%d/%m/%Y"),
        end = as.Date(end, "%d/%m/%Y")
      ) 
    names(out)[[count]] <- 'interventions_info'
  }
  
  # Load Bootstrapped Data
  if (load_bootstrapped_scalars) {
    count <- count + 1
    out[[count]] <- read_csv(paste0(data_path, 'bootstrapped_scalars.csv'),
                                     col_types = paste(rep("d", 12), collapse = ''))
    names(out[[count]])[c(10, 12)] <- c("Pre-Christmas Level 3", "After-Christmas Level 5")
    names(out)[[count]] <- 'boot_lockdown_scalars'
  }
  
  if (load_optim_scalars) {
    count <- count + 1
    out[[count]] <- read_csv(paste0(data_path, 'optim_scalars.csv'),
                                     col_types = "cd")
    names(out)[[count]] <- 'optim_res'
  }
  
  # Load Estimated Costs and Deaths
  if (load_deaths) {
    count <- count + 1
    out[[count]] <- read_csv(paste0(data_path, 'Estimated_Deaths.csv'),
                                     col_types = "cd")
    names(out)[[count]] <- 'est_deaths'
  }
  
  if (load_costs) {
    count <- count + 1
    out[[count]] <- read_csv(paste0(data_path, 'Estimated_Costs.csv'),
                                     col_types = "cd")
    names(out)[[count]] <- 'est_costs'
  }
  
  return(out)
}


