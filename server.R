# Travis Baer
# Fall 2016

library(shiny)
library(dplyr)
library(survival)
library(survminer)
library(KMsurv)
data("larynx")
source("reliability_functions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  event_data <- reactive({
    larynx <- larynx[larynx$time > 0, ] # remove any 0-time events
    larynx$age_interval <- cut(larynx$age, breaks = c(40,60,70,90)) # ~ 30 in each of three buckets
    larynx$stage <- ordered(larynx$stage)
    return(larynx)
  })
  
  surv_obj <- reactive({
    Surv(time = event_data()[["time"]], event = event_data()[["delta"]])
  })
  
  # weibull reactives
  surv_fit_weib <- reactive({
    survreg(surv_obj() ~ 1, dist = "weibull")
  })
  
  rank_pts_weib <- reactive({
    CalculateModKM(select(event_data(), time, delta), params = surv_fit_best()[["icoef"]], 
                   distribution = "weibull")
  })
  
  # exponential reactives
  surv_fit_exp <- reactive({
    survreg(surv_obj() ~ 1, dist = "exponential")
  })
  
  rank_pts_exp <- reactive({
    CalculateModKM(select(event_data(), time, delta), params = surv_fit_exp()[["icoef"]], 
                   distribution = "exponential")
  })
   
  # best fits
  surv_fit_best <- reactive({
    survreg(surv_obj() ~ 1, dist = "weibull")
  })
  
  rank_pts_best <- reactive({
    CalculateModKM(select(event_data(), time, delta), params = surv_fit_best()[["icoef"]], 
                 distribution = "weibull")
  })
  
  # split data
  dists <- reactive({
    # splits data and fits distributions
    to_return <- event_data()
    if(input$split_by == "None"){
      # only fits a single distribution per distribution type if split_by is None
      to_return <- tibble("Data Subset" = "All", 
                          "data" = list(to_return)) %>% mutate(weib = purrr::map(data, ~ FitADist(.)), 
                                                               exp = purrr::map(data, ~ FitADist(., "exp")))
      
      
    } else {
      
    }
  })
  
  # Outputs -----------------------------------------------------------------
  
  output$data_summary <- renderPrint({
    # TODO: replace with table(s) - summary() returns a table to pass to renderTable
    print(paste(nrow(event_data()), "x", ncol(event_data())))  
    print(c("stage - Stage of disease (1=stage 1, 2=stage2, 3=stage 3, 4=stage 4)",
      "time - Time to death or on-study time, months",
"age - Age at diagnosis of larynx cancer",
"diagyr - Year of diagnosis of larynx cancer",
"delta - Death indicator (0=alive, 1=dead)"))
    summary(event_data())
  })
  
  output$split_by_selection <- renderUI({
    selectInput("split_by", label = "Select Factor to Split Data By", width = 250,
                choices = c("None", names(event_data())[c(1,6)]))
  })
  
  output$dist_table_one <- renderTable({
    pval <- BetaTest(surv_fit_weib()$loglik[1], surv_fit_exp()$loglik[1])
    weib_param <- TidySurv(surv_fit_weib()) # parameter_name, estimate, se
    exp_param  <- TidySurv(surv_fit_exp())
    data.frame("Distribution Name" = c("Weibull", "Exponential"),
               "parameter_one"     = c(weib_param$estimate[1], exp_param$estimate[1]),
               "parameter_one_se"  = c(weib_param$se[1], exp_param$se[1]),
               "parameter_two"     = c(weib_param$estimate[2], NA),
               "parameter_two_se"  = c(weib_param$se[2], NA),
               "log_likelihood"    = c(surv_fit_weib()$loglik[1], surv_fit_exp()$loglik[1]),
               "anderson_darling"  = c(CalculateADA(rank_pts_weib()), CalculateADA(rank_pts_exp())),
               "beta_eq_one_pval"  = c(pval, NA))
  })
  
  output$summary_plot_select <- renderUI({
    selectInput("summary_plot_selection", label = "Select a Field to Plot", choices = names(event_data()))
  })
  
  output$distribution_info_table <- renderTable({
    to_return <- surv_fit_weib() %>% ReturnParameters()
  })
  
  output$summary_plot <- renderPlot({
    s <- input$summary_plot_selection
    if(!is.null(s)){
      if(any(class(event_data()[[s]]) %in% c("numeric","integer"))){
        hist(event_data()[[s]], 
             xlab = s, main = "Histogram")
      } else {
        barplot(table(event_data()[[s]]), xlab = s)
      }
    }
  })
  
  output$surv_plot <- renderPlot({
    ggsurvplot(survfit(surv_obj() ~ 1))
  })
  
  output$dist_plot <- renderPlot({
    ranked_points <- CalculateModKM(select(event_data(), time, delta), params = surv_fit_best()[["icoef"]], 
                                    distribution = "weibull")
    parameter_names <- data.frame(value = "All", parameter_name = "All")
    DisplayWeibullPlot(ranked_points, parameter_names, surv_fit_best()[["icoef"]], surv_fit_best()[["dist"]])
  })
  
})

