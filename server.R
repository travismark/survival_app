# Travis Baer
# Fall 2016

library(shiny)
library(tidyverse)
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
    names(larynx)[names(larynx)=="delta"] <- "causal"
    return(larynx)
  })
  
  surv_obj <- reactive({
    Surv(time = event_data()[["time"]], event = event_data()[["causal"]])
  })
  
  # weibull reactives
  surv_fit_weib <- reactive({
    survreg(surv_obj() ~ 1, dist = "weibull")
  })
  
  rank_pts_weib <- reactive({
    CalculateModKM(select(event_data(), time, causal), params = surv_fit_best()[["icoef"]], 
                   distribution = "weibull")
  })
  
  # exponential reactives
  surv_fit_exp <- reactive({
    survreg(surv_obj() ~ 1, dist = "exponential")
  })
  
  rank_pts_exp <- reactive({
    CalculateModKM(select(event_data(), time, causal), params = surv_fit_exp()[["icoef"]], 
                   distribution = "exponential")
  })
   
  # best fits
  surv_fit_best <- reactive({
    survreg(surv_obj() ~ 1, dist = "weibull")
  })
  
  rank_pts_best <- reactive({
    CalculateModKM(select(event_data(), time, causal), params = surv_fit_best()[["icoef"]], 
                 distribution = "weibull")
  })
  
  # split data
  dists <- reactive({
    if(is.null(input$split_by)){
      return(NULL)
    }
    # splits data and fits distributions
    to_return <- event_data()
    if(input$split_by == "None"){
      # only fits a single distribution per distribution type if split_by is None
      to_return <- tibble("subset_name" = "All Data", "data" = list(to_return)) %>% 
        mutate(weib = purrr::map(data, ~ FitADist(.)), 
               exp = purrr::map(data, ~ FitADist(., d = "exp"))) %>%
        mutate(weib_p_val = purrr::map2_dbl(.x = weib, .y = exp, 
                                            .f = ~ BetaTest(-.x$loglik[1], -.y$loglik[1]))) %>%
        mutate(ev_cause = map_int(data, ~sum(.$causal)), ev_susp = map_int(data, ~sum(.$causal==0)))
    } else {
      to_return <- to_return %>% group_by_(input$split_by) %>% nest() %>%
        mutate(weib = purrr::map(data, ~ FitADist(.)), 
               exp = purrr::map(data, ~ FitADist(., d = "exp"))) %>%
        mutate(weib_p_val = purrr::map2_dbl(.x = weib, .y = exp, 
                                            .f = ~ BetaTest(-.x$loglik[1], -.y$loglik[1]))) %>%
        mutate(ev_cause = map_int(data, ~sum(.$causal)), ev_susp = map_int(data, ~sum(.$causal==0)))
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
    if(is.null(dists())){
      return(NULL)
    }
    #pval <- BetaTest(surv_fit_weib()$loglik[1], surv_fit_exp()$loglik[1])
    to_return_wei <- dists() %>% 
      mutate(dist_name = "weibull", rest = purrr::map(weib, ~ AntiTidySurv(.)),
             ranked_points = map2(.x = data, .y = weib, 
                                       .f = ~ CalculateModKM(select(.x, time, causal) %>% 
                                                               as.data.frame(), params = .y[["icoef"]], 
                                                             distribution = "weibull")), 
             anderson_darling = map_dbl(ranked_points, ~ CalculateADA(.))) %>% 
                             unnest(rest)
    to_return_exp <- dists() %>% 
      mutate(dist_name = "exponential", rest = purrr::map(exp, ~ AntiTidySurv(.)), 
             ranked_points = map2(.x = data, .y = exp, 
                                       .f = ~ CalculateModKM(select(.x, time, causal) %>% 
                                                               as.data.frame(), params = .y[["icoef"]], 
                                                             distribution = "exponential")),
             anderson_darling = map_dbl(ranked_points, ~ CalculateADA(.))) %>% unnest(rest)
    
    to_return <- bind_rows(to_return_wei, to_return_exp) %>% arrange_(names(to_return_wei)[1]) %>% 
      select(1, 8, 6, 7, 11:ncol(.), 10, 5)
    
    # weib_param <- TidySurv(surv_fit_weib()) # parameter_name, estimate, se
    # exp_param  <- TidySurv(surv_fit_exp())1
    # data.frame("dist_name" = c("Weibull", "Exponential"),
    #            "param1"     = c(weib_param$estimate[1], exp_param$estimate[1]),
    #            "param1_one_se"  = c(weib_param$se[1], exp_param$se[1]),
    #            "param2"     = c(weib_param$estimate[2], NA),
    #            "param2_se"  = c(weib_param$se[2], NA),
    #            "log_likelihood"    = c(surv_fit_weib()$loglik[1], surv_fit_exp()$loglik[1]),
    #            "anderson_darling"  = c(CalculateADA(rank_pts_weib()), CalculateADA(rank_pts_exp())),
    #            "beta_eq_one_pval"  = c(pval, NA))
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
      if(any(class(event_data()[[s]]) %in% c("numeric", "integer"))){
        #hist(event_data()[[s]], 
        #     xlab = s, main = "Histogram")
        ggplot(data = event_data(), aes(x = s)) + geom_histogram()
      } else {
        barplot(table(event_data()[[s]]), xlab = s)
      }
    }
  })
  
  output$surv_options <- renderUI({
    if(is.null(input$split_by) || input$split_by == "None"){
      return(NULL)
    }
    return(
      tags$div(checkboxInput("surv_split_by", "Split Survival Plot by Factor", value = TRUE),
               conditionalPanel(
                 condition = "input.surv_split_by == true",
                 tags$div(selectInput("surv_split_by_choice", multiple = TRUE,
                             label = paste("Select a", MakeTitleCase(input$split_by)),
                             choices = levels(event_data()[[input$split_by]]),
                             selected = levels(event_data()[[input$split_by]])),
                          style = "padding-left: 20px")
               ),
               conditionalPanel(
                 condition = "input.surv_split_by == true",
                 checkboxInput("surv_split_confint", "Conf. Intervals", value = FALSE),
                 style = "padding-left: 20px"),
        class = "well", style = "display:inline-flex; padding-bottom: 0px; padding-top: 14px")
    )
  })
  
  output$dist_options <- renderUI({
    if(is.null(input$split_by) || input$split_by == "None"){
      return(NULL)
    }
    return(
      tags$div(checkboxInput("dist_split_by", "Split Distribution Plot by Factor", value = TRUE),
               conditionalPanel(
                 condition = "input.dist_split_by == true",
                 tags$div(selectInput("dist_split_by_choice", multiple = TRUE,
                                      label = paste("Select a", MakeTitleCase(input$split_by)),
                                      choices = levels(event_data()[[input$split_by]]),
                                      selected = levels(event_data()[[input$split_by]])),
                          style = "padding-left: 20px")
               ),
               class = "well", style = "display:inline-flex; padding-bottom: 0px; padding-top: 14px")
    )
  })
  
  output$surv_plot <- renderPlot({
    if(is.null(input$surv_split_by) || !input$surv_split_by){
      # if null then data hasn't been split by factor. provide single plot of all data
      # if false user wants to continue to group data for a single surv
      ggsurvplot(survfit(surv_obj() ~ 1))
    } else if (input$surv_split_by){
      if(is.null(input$surv_split_by_choice)){
        return(NULL)
      }
      # if false then user wants to view a line per selected factor level
      
      filtered_data <- event_data()[event_data()[[input$split_by]] %in% input$surv_split_by_choice, ]
      filtered_data[[input$split_by]] <- factor(filtered_data[[input$split_by]])
      split_data <- filtered_data[[input$split_by]]
      surv_obj <- Surv(time = filtered_data[["time"]], event = filtered_data[["causal"]])
      
      ggsurvplot(survfit(surv_obj ~ split_data), 
                 legend.labs = MakeTitleCase(paste(input$split_by,
                                                   levels(split_data))),
                 conf.int = input$surv_split_confint
      )
    }
  })
  
  output$dist_plot <- renderPlot({
    ranked_points <- CalculateModKM(select(event_data(), time, causal), params = surv_fit_best()[["icoef"]], 
                                    distribution = "weibull")
    parameter_names <- data.frame(value = "All", parameter_name = "All")
    DisplayWeibullPlot(ranked_points, parameter_names, surv_fit_best()[["icoef"]], surv_fit_best()[["dist"]])
  })
  
  output$interactive_plot <- renderggiraph({
    to_return <- event_data() %>% ggplot(aes(y = as.numeric(row.names(larynx)), x = time, colour = as.logical(causal))) +
      geom_point_interactive(aes(tooltip = paste0("Diagnosis Year: ", diagyr, "<br>Age: ", age))) + 
      facet_grid(paste("Stage", stage) ~ ., scales = "free_y") + 
      labs(x = "Months until Event", y = "Patient Id") +
      theme_bw() + theme(axis.text.y = element_blank(), strip.background = element_rect("#3c8dbc"),
                         panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                         legend.position = "top") +
      scale_color_viridis(discrete = TRUE, name = "Recorded Death")
    ggiraph(code = print(to_return), width = 1)
  })
  
  output$parametric_regression_select <- renderUI({
    choices = names(event_data())[!names(event_data()) %in% c("time", "causal")]
    selectInput("parametric_regression_selection", label = "Select Factor to Regress By",
                choices = choices)
  })
  
  output$parameteric_regression <- renderPrint({
    if(is.null(input$parametric_regression_selection)){
      return(NULL)
    }
    survreg(surv_obj() ~ event_data()[[input$parametric_regression_selection]], dist = "weibull") %>% 
      summary() %>% print()
  })
  
})

