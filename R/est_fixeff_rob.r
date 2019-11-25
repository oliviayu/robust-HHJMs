#' Estimate fixed parameters using robust inference method
#' 
#' @importFrom nleqslv nleqslv
est_fixeff_rob <- function(param, sigma, other_param,
                           model_object,
                           long.data,
                           Bis, invSIGMA,
                           Silent = TRUE){ 

  p <- length(param)
  linear_pred <- model_object$reg_equation
  resForm <- paste("(", model_object$response, "-(", linear_pred, "))/", names(sigma))
  Tpt <- model_object$robust$Tpt
  delim_val <- model_object$left_censoring$limit_value
  uncensored_idx <- (long.data[, model_object$left_censoring$indicator] == 0)
  sub_data <- long.data[uncensored_idx, ]
  
  # score function/nonlinear equation to be solved
  ff <- function(xx, subB){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param, sigma)
    subB <- subB[uncensored_idx, ]
    # evaluate the residual values
    resVal <- with(sub_data, with(par_val, with(subB, eval(parse(text = resForm)))))
    OKSign <- as.numeric(abs(resVal) <= Tpt)
    negSign <- as.numeric(resVal < -Tpt)
    posSign <- as.numeric(resVal > Tpt)
    # evaluate the linear predictor
    expts <- with(sub_data, with(par_val, with(subB, eval(parse(text = linear_pred)))))
    uSign <- as.numeric(delim_val > expts + Tpt*sigma[[1]])
    bSign <- as.numeric(delim_val < expts - Tpt*sigma[[1]])
    mSign <- as.numeric(delim_val >= expts - Tpt*sigma[[1]])* 
             as.numeric(delim_val <= expts + Tpt*sigma[[1]])
    Phi_d <- pnorm(as.numeric((delim_val - expts)/sigma[[1]]))
    phi_d <- dnorm(as.numeric((delim_val - expts)/sigma[[1]]))
    
    lapply(1:p, function(i){
      term0 <- Vderiv(resForm, names(param)[i])
      term0Val <- - with(sub_data, with(par_val, with(subB, eval(parse(text = term0)))))
      pt1Val <- sum( term0Val*(resVal*OKSign+Tpt*posSign+ (-Tpt)*negSign))

      term1 <- Tpt*Phi_d/(1-Phi_d)*term0Val*bSign
      term2 <- Tpt*term0Val*uSign
      term3 <- (Tpt*pnorm(-Tpt)+phi_d-dnorm(Tpt))/(1-Phi_d)*term0Val*mSign
      term1[bSign==0] <- 0
      term2[uSign==0] <- 0 
      term3[mSign==0] <- 0 
      pt2Val <- sum( term1+term2+term3)
      pt1Val - pt2Val
    }) %>% unlist()
  }
  
  mcmc_size <- length(Bis)
  # approximate score function with MCMC samples
  subject_id <- names(Bis[[1]])[1]
  fy_sum <- function(xx){
    fy <- numeric(p)
    for(i in 1:mcmc_size){
      subB_i <- dplyr::select(long.data, subject_id) %>%
        dplyr::left_join(Bis[[i]], by = subject_id)
      fy <- fy + ff(xx, subB_i)
    }
    fy/mcmc_size
  }
  
  # Start iteration of finding optimal solution
  message <- -1
  M <- 1
  while( !message%in%c(1,2) & M<100){
    str_val0 <- sapply(param, function(x){
      x + rnorm(1,0, min(1, abs(x/5)))
    })
    
    result <- try(nleqslv::nleqslv(str_val0, fy_sum, control = list(maxit = 20),
                                   method = "Newton", jacobian = TRUE), silent=T)
    error_mess <- attr(result, "class")
    
    if(length(error_mess) != 0){ 
      message <- -1
    } else {
      message <- result$termcd
    }
    
    if(Silent == F){cat("M=", M,"\n"); print(result)}
    M <- M +1
  }
  
  if(message %in% c(1,2)){
    gamma <- result$x
    names(gamma) <- names(param)
    fval <- result$fvec
    jac <- result$jac
  } else {
    stop("Iteration stops because param parameters can not be successfully estimated.")  
  }

  list(gamma = as.list(gamma), fval = fval, jac = jac)
}