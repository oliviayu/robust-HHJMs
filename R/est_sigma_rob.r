#' Robust estimate of standard error of residual term
#'
est_sigma_rob <- function(param, lower = 0, upper = 100,
                          fixed_param, other_param = NULL,
                          model_object, long.data, 
                          Bis, Silent){
  
  linear_pred <- model_object$reg_equation
  resForm <- paste("(", model_object$response, "-(", linear_pred, "))/", names(param))
  Tpt <- model_object$robust$Tpt
  delim_val <- model_object$left_censoring$limit_value
  uncensored_idx <- (long.data[, model_object$left_censoring$indicator] == 0)
  sub_data <- long.data[uncensored_idx, ]

  # robust score function
  gr <- function(xx, subB){
    # all parameter values
    par_val <- c(Vassign(names(param), xx), fixed_param, other_param)
    subB <- subB[uncensored_idx, ]
    # evaluate the residual values
    resVal <- with(sub_data, with(par_val, with(subB, eval(parse(text = resForm)))))
    OKSign <- as.numeric(abs(resVal) <= Tpt)
    negSign <- as.numeric(resVal < -Tpt)
    posSign <- as.numeric(resVal > Tpt)
    # evaluate the linear predictor
    expts <- with(sub_data, with(par_val, with(subB, eval(parse(text = linear_pred)))))
    uSign <- as.numeric(delim_val > expts + Tpt*xx)
    bSign <- as.numeric(delim_val < expts - Tpt*xx)
    mSign <- as.numeric(delim_val >= expts - Tpt*xx)*
      as.numeric(delim_val <= expts + Tpt*xx)
    
    # the derivative part
    pt1 <- sum(sapply(resVal, function(x){max(min(x, Tpt), -Tpt)^2}))
    # the expectation part/ correction term
    d_tilde <- as.numeric((delim_val - expts)/xx)
    Phi_d <- pnorm(d_tilde)
    Phi_c <- pnorm(Tpt)

    term1 <- Tpt^2
    term2 <- (2*Tpt^2*(1 - Phi_c) - Tpt^2*Phi_d + 2*Phi_c -
                2*Tpt*dnorm(Tpt) - 1)/(1 - Phi_d)
    term3 <- (Tpt^2*(1 - Phi_c) + Phi_c - Tpt*dnorm(Tpt) - Phi_d +
                d_tilde*dnorm(d_tilde))/(1 - Phi_d)
    pt2 <- term1*sum(uSign) + sum(term2[which(bSign == 1)]) + 
      sum(term3[which(mSign == 1)])

    pt1 - pt2  
  }
  
  
  mcmc_size <- length(Bis)
  # approximate score function with MCMC samples
  subject_id <- names(Bis[[1]])[1]
  gr_sum <- function(xx){
    fy <- numeric(1)
    for(i in 1:mcmc_size){
      subB <- dplyr::select(long.data, subject_id) %>%
        dplyr::left_join(Bis[[i]], by = subject_id)
      fy <- fy+ gr(xx, subB)
    }
    fy/mcmc_size
  }

  # Start iteration of finding optimal solution
  message <- 'try-error'
  M <- 1
  while(message == 'try-error' & M < 10){
    result <- try(uniroot(gr_sum, c(lower, upper)), silent=T)
    # TODO: investigate if this is necessary
    lower <- lower + 0.02
    upper <- upper - 0.02
    M <- M + 1
    message <- class(result)
    # cat("lower =", lower, '\n')
    # cat("upper =", upper, '\n')
  }

  if(message != 'try-error'){
    as.list(Vassign(names(param), result$root))
  } else {
    print(result)
  }
}