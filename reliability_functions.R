library(tibble)

FitADist <- function(x, time = "time", causal = "causal", d = "weibull"){
  # fit a distribution for a dataset
  so <- Surv(time = x[[time]], event = x[[causal]])
  return(survreg(so ~ 1, dist = d))
}

BetaTest <- function(wll, ell){
  # test if beta different than one
  # https://en.wikipedia.org/wiki/Likelihood-ratio_test
  test_stat <- 2*(abs(wll)-abs(ell))
  return(1 - pchisq(test_stat, 1))
}

TidySurv <- function(x){
  # like broom package function - to get distribution parameter stats
  # Args:
  #  x: survfit object
  if(x$dist == "weibull"){
    scale = exp(x$icoef[1]) %>% unname()
    shape = 1/exp(x$icoef[2]) %>% unname()
    se = diag(x$var)
    return(tibble::tibble(
      "parameter_name" = c("scale","shape"),
      "estimate" = c(scale, shape),
      "se" = c(sqrt(se[1]) * scale, sqrt(se[2]) * shape) %>% unname()
    ))
  } else {
    mu = exp(x$icoef[1]) %>% unname()
    return(tibble::tibble(
      "parameter_name" = "mean",
      "estimate" = mu,
      "se" = sqrt(x$var) * mu
    ))
  }
}

AntiTidySurv <- function(x){
  # like broom package function - to get distribution parameter stats
  # Args:
  #  x: survfit object
  if(x$dist == "weibull"){
    scale = exp(x$icoef[1]) %>% unname()
    shape = 1/exp(x$icoef[2]) %>% unname()
    se = diag(x$var)
    return(tibble::tibble(
      "param1" = scale,
      "param1_se" = sqrt(se[1]) * scale %>% unname(),
      "param2" = shape,
      "param2_se" = sqrt(se[2]) * shape %>% unname(),
      "dist_mean" = scale*gamma(1+1/shape),
      "loglik" = -x$loglik[1]
    ))
  } else {
    mu <- exp(x$icoef[1]) %>% unname()
    return(tibble::tibble(
      "param1" = mu,
      "param1_se" = (sqrt(x$var) * mu) %>% as.numeric(),
      "param2" = NA,
      "param2_se" = NA,
      "dist_mean" = mu,
      "loglik" = -x$loglik[1]
    ))
  }
}

ReturnParameters <- function(xs, se = FALSE){
  # return parameters as list
  if(xs[["dist"]] == "weibull"){
    to_return <- list(
      scale = exp(xs$icoef[1]) %>% unname(),
      shape = 1/exp(xs$icoef[2]) %>% unname()
    )
    if(se){
      return()
    } 
  } else {
    to_return <- list(
      mean = exp(xs$icoef[1]) %>% unname()
    )
    if(se){
      return()
    }
  }
  return(to_return)
}


CalculateModKM <- function(interval_value, params, distribution, source_package = "survreg") {
  # Calculates Modified Kaplan Meier and distribution values for each event in the distribution
  #  Also re-parameterizes survreg distribution parameters
  # Args:
  #  interval_value: simplified interval data for the distribution in a 2-column data frame
  #   > interval value
  #   > causal or non-causal (1/0)
  #  params: parameters from the fit
  #  distribution: string of distribution type, like weibull
  #  source_pacakge: string of survreg or fitdistcens
  # Returns:
  #  a data frame with the modified values, rank, and ages
  
  # TODO: Use classes to automate these calls, i.e. don't have to check source_package
  if (source_package %in% "survreg") {
    if (distribution %in% "weibull") {
      scale <- unname(exp(params[1]))
      shape <- unname(1/exp(params[2]))
    } else if (distribution %in% "exponential") {
      scale <- unname(exp(params))
      shape <- 1
    }
  } else if (source_package %in% "fitdistcens") {
    if (distribution %in% "weibull") {
      scale <- unname(params[2])
      shape <- unname(params[1])
    } else if (distribution %in% "exponential") {
      scale <- unname(params)
      shape <- 1
    }
  } else {
    stop("Failed attempting to calculate kaplan meier statistics for a fit 
         of type other than fitdistcens or survreg")
  }
  
  # Extract and sort interval data into two data frames with the distribution probability
  causal_data <- data.frame("interval_value" = sort(interval_value[, 1][interval_value[, 2]==1]),
                            "fit_probability" = pweibull(sort(interval_value[, 1][interval_value[, 2]==1]),shape,scale),
                            "causal" = 1)
  # If some suspension data exists, create the data frame and combine with causal (order matters)
  if (sum(interval_value[, 2]==0)>0) {
    suspension_data <- data.frame("interval_value" = sort(interval_value[, 1][interval_value[, 2]==0]),
                                  "fit_probability" = pweibull(sort(interval_value[, 1][interval_value[, 2]==0]),shape,scale),
                                  "causal" = 0)
    interval_value <- bind_rows(causal_data, suspension_data)
  } else {
    interval_value <- causal_data
  }
  # Order by event time and causal b/f suspension, then add row name as rank
  interval_value <- arrange(interval_value, interval_value, -causal) %>% rownames_to_column(var = "rank") 
  interval_value$rank <- as.integer(interval_value$rank)
  
  # Calculate Modified Kaplan-Meier rank (just for the non-censored events)
  #  using nrow(interval_value) allows calculating reverse-rank from rank
  max_rank <- nrow(interval_value)
  
  # Calculate:
  # p_prime - Normal Kaplan Meier statistic (matches JMP and fitsurv from survival package)
  # p       - Modified Kaplan Meier statistic
  interval_value <- mutate(interval_value, 
                           p_intmd =   (max_rank - row_number())^causal / 
                             (max_rank - row_number() + 1)^causal,
                           p_prime = 1 - cumprod(p_intmd), 
                           prev_p_prime = lag(p_prime, 1L, default = 0),
                           p = (1 - ((1 - p_prime) + (1 - prev_p_prime)) / 2)) %>%
    select(-prev_p_prime,-p_intmd)
  # using mutate with cumprod is much faster than vapply and prod(1:rank)
  # other uses of mutate are mostly for elegance
  
  return(interval_value)
  } # end CalculateModKM

CalculateADA <- function(ranked_points) {
  # Calculate's a distribution's anderson darling adjusted statistic
  # Args
  #  ranked_points: interval data with modified Kaplan Meier ranks (empirical dist) and distribution fit
  #   p is modified kaplan meier statistic, p_prime is standard kaplan meier statistic
  # Returns
  #  the anderson darling adjusted statistic (a double)
  
  # Only take the plotted (noncensored) points
  ranked_points      <- ranked_points[ranked_points$causal==1,] 
  lgth               <- nrow(ranked_points) # number of noncensored points
  ranked_points$rank <- seq_len(lgth) # recalculate the rank
  
  ranked_points <- mutate(ranked_points, 
                          prev_fit_prob = lag(fit_probability, 1L, default = 0),
                          prev_p = lag(p, 1L, default = 0),
                          ada_contribution = -fit_probability - log(1-fit_probability) + 
                            prev_fit_prob + log(1-prev_fit_prob) + 
                            2*log(1-fit_probability)*prev_p - 2*log(1-prev_fit_prob)*prev_p +
                            log(fit_probability)*prev_p^2 - log(1-fit_probability)*prev_p^2 - 
                            # log(0) is infinity, so must use ifelse here:
                            ifelse(prev_fit_prob==0,0,log(prev_fit_prob)*prev_p^2) +
                            log(1-prev_fit_prob)*prev_p^2
  )
  
  # the n+1th row value doesn't belong in the data frame
  final_p <- ranked_points$p[lgth]
  final_fit <- ranked_points$fit_probability[lgth]
  
  final_point <- -(1-1E-12) - log(1-(1-1E-12)) + 
    ranked_points$fit_probability[lgth] + log(1 - final_fit) + 
    2*log(1-(1-1E-12))*final_p - 2*log(1-final_fit)*final_p + 
    log(1-(1E-12))*final_p^2 - log(1-(1-(1E-12)))*final_p^2 -
    log(final_fit)*final_p^2 + log(1-final_fit)*final_p^2
  
  return(lgth*(sum(ranked_points$ada_contribution,final_point)))
} # end CalculateADA

MakeTitlePretty <- function(plot_title, parameter_separator, max_char) {
  # Recursive function to add line breaks to a title for weibull plots
  # Args:
  #   plot_title: original or modified plot title
  #   parameter_separator: character string with which parameters are spaced apart in the title
  #   max_char: ideal number of characters to fit on a plot row
  # Returns
  #   a plot title or part of a plot title
  
  ends_vector <- gregexpr(parameter_separator, plot_title)[[1]] # where parameters end on the same line
  # Four cases:
  #  1) plot title is less than the character limit: return title w/o recursion
  #  2) there's only one parameter on this so: return title w/o recursion
  #  3) the first parameter is over the limit - split it into its own line and recurse
  #  4) some other parameter is over the limit - make everything before this parameter its own line and recurse
  
  # case 1 - short line
  if (nchar(plot_title) <= max_char) { 
    return(plot_title)
  }
  
  # case 2 - one param
  if (ends_vector[1] == -1) { 
    return(plot_title)
  }
  
  ends_vector <- c(ends_vector,nchar(plot_title))
  too_long_param <- which(ends_vector > max_char)[1] # first parameter after max_char characters
  # case 3 - split is on first param
  if(too_long_param == 1) {
    split_this_many_before <- 0 # split the first (only) separator that's over the limit
  } else {
    # case 4
    split_this_many_before <- 1
  }
  
  plot_line       <- substr(plot_title, start = 1, stop = ends_vector[too_long_param-split_this_many_before]-1)
  remaining_title <- substr(plot_title, 
                            start = ends_vector[too_long_param-split_this_many_before] + nchar(parameter_separator), 
                            stop = nchar(plot_title))
  # recurse
  plot_title <- paste(plot_line, MakeTitlePretty(remaining_title, parameter_separator, max_char),sep = "\n")
  
  return(plot_title)
}

DisplayWeibullPlot <- function(ranked_points, dist_param_names, params, distribution_type = "weibull") {
  # Save a probability plot of the events and their fit
  #   overlayed with a histogram of the censored values
  #   uses weibull-scale log-log axes:
  #     x is log(time)
  #     y is log(log(1/(1-Unreliability(t))))
  # Args:
  #   ranked_points: modified Kaplan Meier ranked points for all intervals
  #   parameter_names: sample of original data frame with classifiers
  # params: two parameters for the weibull distribution
  # catgs: grouping categories passed to ddply to define which factors to not use "ALL"
  # uses Modified Kaplan-Meier method for plotting uncensored and censored data
  
  # Extract Parameters
  if (distribution_type %in% "weibull") {
    params[1] <- (exp(params[1])) # scale
    params[2] <- (1/exp(params[2])) # shape
    params[3] <- params[1]*gamma(1+1/params[2])
    names(params) <- c("scale", "shape", "dist_mean")
  } else if (distribution_type %in% "exponential") {
    params <- exp(params)
    names(params) <- "mean"
  } else (
    stop("Trying to make a weibull plot from a distribution other than weibull or exponential")
  )
  
  # Build the plot title as a string
  #   hopefully there are no
  ## TODO: find a better way to order these parameter names
  dist_param_names    <- arrange(dist_param_names, parameter_name) # order alphabetically
  parameter_separator <- "; "
  plot_title <- paste(dist_param_names$parameter_name,dist_param_names$value,sep=":",collapse=parameter_separator)
  # try to split into new lines about every 40 characters
  plot_title <- MakeTitlePretty(plot_title, parameter_separator, max_char = 40)
  distribution_type_title <-paste0(toupper(substr(distribution_type,1,1)),
                                   substr(distribution_type,2,nchar(distribution_type)))
  
  # Plot
  op <- par() # save the plot parameters
  # dots on plot
  ## TODO: make the filename more informative
  #png(file = paste(plot_dir,"WeibullPlot_", distribution_type, "_", 
  #                 dist_param_names$interval_parameter_set_id[1], ".png", sep=""), 
  #    width=5.96, height=4.54, units="in", res=144)
  par(mar=c(3.6, 4, 4, 2) + 0.1) # push the plot down a little
  weibplot(ranked_points[ranked_points$causal==1, ]$interval_value,
           ranked_points[ranked_points$causal==1, ]$p,
           #forcexlim=c(.9, log10(max(ranked_points$interval_value)) + .02),
           #forceylim=c(-7, 0), 
           xlab="", ylab="Event Probability", col=rgb(255/255, 0/255, 0/255), pch=16, cex = 1.2,
           main=paste0("Removal Plot with ", distribution_type_title, " Fit\n", plot_title))
  title(sub="Time", line=2)
  abline(h=c(log(-log(1-c(.001, .01, .1)))), v=c(c(seq(-2,3))), col=rgb(1, 0, 0)) # Red
  abline(h=c(log(-log(1-c(.003, .02, .05, .25, .5, .75, .9, .96, .99, .999)))), col=rgb(38/255, 201/255, 38/255))
  abline(v=c(log10(seq(0.1, 0.9, 0.1)), log10(seq(1, 9)), log10(seq(20, 90, 10)), log10(seq(200, 900, 100)), 
             log10(seq(2000, 9000, 1000))), col=rgb(38/255, 201/255, 38/255)) # Green
  # fitted line
  if (distribution_type %in% "weibull") {
    line.41 <- qweibull(seq(0.0001, .99, .0005), scale = params[1], shape = params[2], log.p=F)
  } else { # exponential
    line.41 <- qweibull(seq(0.0001, .99, .0005), scale = params[1], shape = 1, log.p=F)
  }
  
  lines(x=log10(line.41), y=log(-log(1-seq(.0001, .99, .0005))))
  # legend
  if (distribution_type %in% "weibull") {
    legend_text <- c(paste("Shape: ", round(params[2], digits=2)), paste("Scale: ",round(params[1], digits=2)), 
                     paste("Mean: ", round(params[3], digits=2)), 
                     paste("Observed: ", nrow(ranked_points[ranked_points$causal==1, ])), 
                     paste("Censored: ", nrow(ranked_points[ranked_points$causal==0, ])))
  } else { 
    legend_text <- c(paste("Exponential\nMean: ", round(params[1], digits=3)), 
                     paste("Observed: ", nrow(ranked_points[ranked_points$causal==1, ])), 
                     paste("Censored: ", nrow(ranked_points[ranked_points$causal==0,])))
  }
  legend(x="topleft", legend=legend_text, cex=.85, ncol=1)#, inset=c(-0.05,-0.03))
  
  #histogram (if there are censored values)
  if(nrow(ranked_points[ranked_points$causal==0, ])==0) { # just plot a zero at median of non-censored hours
    text(x=log10(median(ranked_points[ranked_points$causal==1, ]$interval_value)), y=-7, labels="0")
  } else { # plot the histogram
    par(new=T, mar=c(2.9, 4.1, 4.1, 2.1))
    hist_holder <- hist(log10(ranked_points[ranked_points$causal==0, ]$interval_value), 
                        breaks="Sturges",plot=F) # to get the bin heights and counts
    hist_plot   <- hist(log10(ranked_points[ranked_points$causal==0, ]$interval_value), 
                        breaks="Sturges", #xlim=c(.9, log10(max(ranked_points$interval_value)) + .02), 
                        ylim=c(0, max(hist_holder$counts)*2.4), border=1, 
                        col=rgb(149/255, 184/255, 251/255), ylab=NULL, xlab=NULL, 
                        main=NULL, labels=T, axes=F, plot=T)
  } # end if there is enough data to plot the histogram
  # redraw the dots
  par(new=T, mar=c(3.6, 4, 4, 2) + 0.1)
  weibplot(ranked_points[ranked_points$causal==1, ]$interval_value,
           ranked_points[ranked_points$causal==1, ]$p,
           points = TRUE, xlab="", ylab="",
           col=rgb(255/255, 0/255, 0/255), pch=16, cex = 1.2)
  options(warn=-1)
  par(op)
  options(warn=0)
  #invisible(dev.off())
} # end plot Function

###plotting functions###
# fix the scale to Weibull:
# draw the plot:
weibplot <- function(x, y, points = FALSE, log='xy', ..., forceylim=c(0,0), forcexlim=c(0,0)) {
  x <- log(x,10)
  y <- log(-log(1-y))
  xlg <- TRUE # hard-coded for now
  ylg <- TRUE
  yl <- ifelse(forceylim==c(0,0),range(y),forceylim)
  xl <- ifelse(forcexlim==c(0,0),range(x),forcexlim)
  plot(x,y,...,axes=FALSE,ylim=yl,xlim=xl)
  if(!points){
    if(xlg){drawlogaxis(1,xl)}else{axis(1, at=pretty(xl), labels=pretty(xl))}
    if(ylg){drawweibaxis()}else{axis(2, at=pretty(yl), labels=pretty(yl))}
  }
  graphics::box()
}

# draw the axes:
drawlogaxis <- function(side,range)  {
  par(tck=0.02)
  mlog <- floor(min(range))
  Mlog <- ceiling(max(range))
  SeqLog <- c(mlog:Mlog)
  Nlog <- (Mlog-mlog)+1
  axis(side,at=SeqLog,labels=10^SeqLog)
  ats <- log(seq(from=2,to=9,by=1),10)
  mod <- NULL
  for(i in SeqLog)
  {
    mod <- c(mod,rep(i,length(ats)))
  }
  ats <- rep(ats,Nlog)
  ats <- ats+mod
  par(tck=0.02/3)
  axis(side,at=ats,labels=NA)
}

drawweibaxis <- function()  {
  par(tck=0.02)
  SeqWeib <- c(.001,.003,.01,.02,.05,.1,.25,.5,.75,.9,.96,.99,.999)
  axis(2, labels=SeqWeib, at=(log(-log(1-SeqWeib))), las=2)
}
### end plotting functions ###