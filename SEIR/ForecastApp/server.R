### Server for Shiny App

server <- function(input, output, session) {
  
  ####-- 1. Info ------------------------------------------------####
  
  pp_url <- a("Jaouimaa et. al. (2021)", href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0260632")
  output$pp_link <- renderUI({
    tagList("This app was created as supplemental material for the following research paper:", 
            pp_url)
  })
  
  sfi_url <- a("SFI Homepage", href="https://www.sfi.ie/")
  output$sfi_link <- renderUI({
    tagList(sfi_url)
  })
  
  output$General_Intro <- renderText({
    "This app allows users to create 8 week forecasts based on our age structured SEIR model for modelling COVID-19 case counts in Ireland.
    These models are very sensitive to our modelling assumptions, and so should only be considered an uncertain projection."
  })
  
  output$Acknowledgement <- renderText({
    "This research was funded by Science Foundation Ireland (SFI)."
  })
  
  ####-- 2. Forecast Settings ----------------------------------------------####
  
  # tmax
  tmax_forecast <- reactive({
    as.numeric( difftime(input$start_date + 55, as.Date('2020-02-29'), units = "days") )
  })
  
  # Intervention correction
  intinfo3 <- reactive({
    intervention_adjust(interventions_info, input$start_date - 1)
  })
  
  # Selecting ages
  age_sel_fc <- reactive({
    age_vec <- numeric(0)
    if ("0 - 14" %in% input$forecast_age_sel) {age_vec <- 1:3}
    if ("15 - 24" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 4:5)}
    if ("25 - 34" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 6:7)}
    if ("35 - 44" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 8:9)}
    if ("45 - 54" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 10:11)}
    if ("55 - 64" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 12:13)}
    if ("65 - 74" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 14:15)}
    if ("75+" %in% input$forecast_age_sel) {age_vec <- c(age_vec, 16)}
    age_small_vec <- which(forecast_age_groups %in% input$forecast_age_sel)
    list(large_age_vec = age_vec, small_age_vec = age_small_vec)
  })
  
  # Lockdown Info
  linfo_forecast <- reactive({
    date_start_seq <- seq.Date(input$start_date, input$start_date + 55, 14)
    chosen_levs <- c(input$res1, input$res2, input$res3, input$res4)
    lockdown_forecast <- rbind(intinfo3(),
                               data.frame(start = date_start_seq, 
                                          end = date_start_seq + 13,
                                          policy = chosen_levs))
    optim_dat <- merge(lockdown_forecast, optim_res, by = "policy")
    merge(optim_dat, boot_scales_t, by = "policy")
  })
  
  # Simulate SEIR Model
  dub_forecasters <- eventReactive(input$fit_forecast, {
    
    tmax <- tmax_forecast()
    linfo <- linfo_forecast()
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    clusterExport(
               cl, varlist=c("SEIR_model_simulation", "def_pars", 
                             "contacts", "dub_xstart", "dub_population",
                             "dub_def_beta", "lsoda", "SEIR_model_D", 
                             "tmax", "linfo"),
               envir=environment())
    
    All_Runs <- foreach(i = 5:ncol(linfo), .combine = rbind) %dopar% {
      SEIR_model_simulation(pars = def_pars,
                            contacts_ireland = contacts,
                            dateStart = as.Date('2020-02-29'),
                            startval = dub_xstart,
                            POP = dub_population,
                            beta = dub_def_beta,
                            tmax = tmax, 
                            lockdown_information = linfo[, c(2:3, i)])$solution
    }
    
    All_Runs <- split(All_Runs, All_Runs$time)
    UL <- foreach(x = All_Runs, .combine = rbind,
                  .final = as.data.frame) %dopar% {
      sapply(x, quantile, probs = 0.975, names = FALSE)
    }
    LL <- foreach(x = All_Runs, .combine = rbind,
                  .final = as.data.frame) %dopar% {
      sapply(x, quantile, probs = 0.025, names = FALSE)
    }
    rm(All_Runs)
    stopCluster(cl)
    
    MID <- SEIR_model_simulation(pars = def_pars,
                                 contacts_ireland = contacts,
                                 dateStart = as.Date('2020-02-29'),
                                 startval = dub_xstart,
                                 POP = dub_population,
                                 beta = dub_def_beta,
                                 tmax = tmax, 
                                 lockdown_information = linfo[, 2:4])$solution
    
    list(MID = MID, UL = UL, LL = LL)
    
  })

  sel_deaths <- eventReactive(input$fit_forecast, {
    UL <- comp_deaths(dub_forecasters()$UL)
    MID <- comp_deaths(dub_forecasters()$MID)
    LL <- comp_deaths(dub_forecasters()$LL) 
    data.frame(UL = UL, MID = MID, LL = LL, 
               row.names = forecast_age_groups)
  })
    
  sel_deaths_age_groups <- reactive({
    sel_deaths()[age_sel_fc()$small_age_vec, ]
  })
  
  ### Deaths
  output$deaths <- renderTable({
    deaths <- rep(0, 3)
    deaths[1] <- sum( sel_deaths_age_groups()$UL )
    deaths[2] <- sum( sel_deaths_age_groups()$MID )
    deaths[3] <- sum( sel_deaths_age_groups()$LL )
    data.frame(Type = c("Upper", "Forecast", "Lower"), Deaths = deaths)
  }, align = "l", digits = 0)
  
  ### Costs
  output$Cost_12 <- renderText({ 
    cost <- round( est_costs$Estimated_Costs[est_costs$Level == input$res1], 2 )
    paste(cost, "b€/day")
  })
  
  output$Cost_34 <- renderText({ 
    cost <- round( est_costs$Estimated_Costs[est_costs$Level == input$res2], 2 )
    paste(cost, "b€/day")
  })
  
  output$Cost_56 <- renderText({ 
    cost <- round( est_costs$Estimated_Costs[est_costs$Level == input$res3], 2 )
    paste(cost, "b€/day")
  })
  
  output$Cost_78 <- renderText({ 
    cost <- round( est_costs$Estimated_Costs[est_costs$Level == input$res4], 2 )
    paste(cost, "b€/day")
  })
  
  output$TotalCostBox <- renderInfoBox({
    sel_levs <- data.frame(Level = c(input$res1, input$res2, input$res3, input$res4))
    sel_costs <- left_join(sel_levs, est_costs, by = "Level")
    
    infoBox(title = "Expected Cost", color = "red", icon = icon("euro", verify_fa = FALSE),
            subtitle = "Billion Euros",
            value = round( sum(sel_costs$Estimated_Costs * 14), 2 ))
  })
  
  ### Output
  dub_out <- reactive({
    Map(comp_sel, x = dub_forecasters(),
        MoreArgs = list(y = input$Disp_comp, z = age_sel_fc()$large_age_vec))
  })
  
  date_place <- eventReactive(input$fit_forecast, {
    input$start_date + seq(1-N_known, 55)
  })
  
  segs <- eventReactive(input$fit_forecast, {
    seq(input$start_date, input$start_date + 42, 14)
  })
  
  output$forecast_plot_dub <- renderPlotly({
    
    plot_ly(x = ~date_place(), y = ~dub_out()$UL, name = 'Upper Limit', type = 'scatter', 
            mode = 'lines', line = list(color = 'blue', dash = "dash")) %>%
      add_trace(y = ~dub_out()$LL, name = "Lower Limit", mode = 'lines', 
                fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)',
                line = list(color = 'blue')) %>%
      add_trace(y = ~dub_out()$MID, name = "Forecast", mode = 'lines',
                line = list(color = 'darkorange')) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$MID[1:N_known], name = "Forecast", 
                mode = 'lines',
                line = list(color = 'darkorange', dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$LL[1:N_known], name = "Lower Limit", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$UL[1:N_known], name = "Upper Limit", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_segments(x = segs(), xend = segs(), y = min(dub_out()$LL), yend = max(dub_out()$UL),
                   line = list(color = 'black', dash = "solid", width = 0.5)) %>%
      layout(showlegend = FALSE, xaxis = list(title = "Date"), 
             yaxis = list(title = input$Disp_comp, range = c(0, max(dub_out()$UL))),
             paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)')
    
  })
 
  output$forecast_exp_deaths <- renderPlotly({
    
    deaths <- as.data.frame(sel_deaths())[age_sel_fc()$small_age_vec, ]
    xax <- forecast_age_groups[age_sel_fc()$small_age_vec]
    
    plot_ly(x = xax, y = deaths$MID, type = 'scatter', mode = 'markers',
            error_y = list(type = "data", symmetric = FALSE, 
                  arrayminus = deaths$MID - deaths$LL, 
                  array = deaths$UL - deaths$MID))
    
  })
   
}

