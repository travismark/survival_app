# Travis Baer
# Fall 2016

library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram
#shinyUI(
dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
  
  # Application title
  titlePanel("Survival Example"),
  
  fluidPage(
    "Data Source: Klein and Moeschberger (1997) Survival Analysis Techniques for Censored and truncated data,
Springer",
    fluidRow(
      column(width = 12, 
           box(width = NULL, status = "primary", collapsible = TRUE, 
                 solidHeader = TRUE, title = "Data Selection",
              column(width=8, 
                     tags$div(selectInput("dataset_name", label = "Select a Dataset", 
                                 choices = "larynx", width = 250), 
                              uiOutput("split_by_selection",  
                                       style = "padding-left: 10px"), 
                              class = "well", style = "display: flex"),
                     #box(h4("Larynx Dataset"), width = NULL, status = "primary", 
                         verbatimTextOutput("data_summary")),
              column(width=4,
                     #box(width = NULL, status = "primary", 
                     uiOutput("summary_plot_select"),
                     plotOutput("summary_plot")
                     #)
                     )
             )
      )),
    fluidRow(column(width = 12,
                    box(width = NULL, status = "primary", title = "Distribution Fits", 
                        solidHeader = TRUE, collapsible = TRUE,
                        tableOutput("dist_table_one"))
                    )),
    fluidRow(
      column(width = 6,
      box(title = "Survival Plot", solidHeader = TRUE, 
          width = NULL, status = "primary",
      plotOutput("surv_plot"))
      ),
      column(width = 6,
      # Show a plot of the generated distribution
      box(title = "Distribution Plot", width = NULL, 
          status = "primary", solidHeader = TRUE,
            plotOutput("dist_plot")
  ))
))
))


