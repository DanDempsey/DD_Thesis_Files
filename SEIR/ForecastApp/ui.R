### UI for Shiny App

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Forecast Settings", tabName = "forecast", icon = icon("wrench")),
    menuItem("Info", tabName = "info", icon = icon("info"))
  )
)

body <- dashboardBody(
  tabItems(
    #Tab 1: Info
    tabItem(tabName = "info",
            box(title = "About This Application", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("General_Intro"),
              uiOutput("pp_link")
            ),
            box(title = "Acknowledgement", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("Acknowledgement"),
              uiOutput("sfi_link")
            )
    ),
    #Tab 2: Forecast Settings
    tabItem(tabName = "forecast",
            fluidRow(
              box( title = "8 Week Compartment Projection", width = 10, solidHeader = TRUE,
                   status = "primary",
                   plotlyOutput("forecast_plot_dub") %>%
                                withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
              box( title = "Display Options", width = 2, solidHeader = TRUE, status = "warning",
                 selectInput(inputId = "Disp_comp", label = "Compartment:",
                             choices = list("Susceptible", "Exposed", "All Infected", 
                                            "Symptomatic Infected", "Pre-symptomatic Infected", 
                                            "Asymptomatic Infected", "Removed"),
                             selected = "Infected:All"),
                 dateInput(inputId = "start_date", label = "Start Date:", 
                           value = as.Date('2021-02-01'), min = as.Date('2020-02-29'), 
                           max = as.Date('2021-02-01')),
                 checkboxGroupInput(inputId = "forecast_age_sel", label = "Age Groups:", 
                                      choices = forecast_age_groups, selected = forecast_age_groups)
                 )
            ),
            fluidRow(
              infoBox("", "Select Lockdown Levels:", 
                      icon = icon("chevron-right"), width = 2, fill = TRUE),
              box( title = "Weeks 1 & 2:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res1", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_12")
                   ),
              box( title = "Weeks 3 & 4:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res2", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_34")
              ),
              box( title = "Weeks 5 & 6:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res3", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_56")
              ),
              box( title = "Weeks 7 & 8:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res4", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_78")
              ),
              box( title = "Forecast Button", width = 2, solidHeader = TRUE, status = "success", 
                   actionButton("fit_forecast", label = "Create Forecast", icon = icon("play-circle"))
              )
            ),
            fluidRow(
              box( title = "8 Week Projected Deaths by Age Group", width = 10, solidHeader = TRUE,
                   status = "danger",
                   plotlyOutput("forecast_exp_deaths") %>%
                                withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
              box( title = "Total Projected Deaths", width = 2, solidHeader = TRUE, status = "danger",
                     tableOutput('deaths') ),
              infoBoxOutput('TotalCostBox', width = 2)
            )
    )
  )
)

# Put them together into a dashboardPage
dashboardPage( skin = "blue",
  dashboardHeader(title = "Ireland COVID-19 Incidence Modelling",
                  titleWidth = 380),
  sidebar,
  body
)

