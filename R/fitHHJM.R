#' Fit robust joint model based on h-likelihood
#' 
#' This function is for robust joint modelling of multiple (generalized) 
#' linear mixed-effect models and a survival model (Cox PH or Weibull).
#' The continuous longitudinal data may be left-censored due to lower 
#' limit of quantification.
#' @import matrixcalc
#' @import magrittr
#' @param glmeObject a list of sublists. Each sublist provides the necessary
#' information of an LME or GLME model.
#' @param survObject a list that provides the necessary information of 
#' a survival model.
#' @param long.data longitudinal data
#' @param surv.data survival data
#' @param subject_id
#' @param event_time
#' @export
fitHHJM <- function(glmeObjects, survObject, long.data, surv.data, subject_id,
                    randeff_info = list(distribution = "t-dist", degree = 3),
                    mcmc_size = 100, burn_in = 100, thin = 20,
                    itertol = 1e-3, Ptol = 1e-2, iterMax = 10, Silent = T){
  # get log hlikelihood
  long_config <- lapply(glmeObjects, get_log_density_glm)
  long_log_hlike <- lapply(long_config, function(x){x$log_density}) %>% paste(collapse = "+")
  surv_config <- get_log_density_ph(survObject)
  surv_log_hlike <- surv_config$log_density
  RespLog <- list(long_log_hlike, surv_log_hlike)
  
  survFit <- survObject$initial_fit
  if(class(survFit) == "coxph"){
    hazard_est <- basehaz(survFit) %>%
      tibble::add_column(step_hazard = c(1e-16, pmax(1e-16, diff(.$hazard, 1))), .before = 1)
    names(hazard_est) <- c(surv_config$additional_param, survObject$response)
    surv.data <- merge(surv.data, hazard_est, by.all = survObject$response) %>%
      dplyr::arrange_at(subject_id)
    fixed_param <- c()
  } else if(class(survFit) == "survreg" & survObject$dist == "weibull"){
    weibPar = c(-summary(survFit)$coeff[1]/summary(survFit)$scale,
                1/summary(survFit)$scale) %>%
      Vassign(surv_config$additional_param, .)
    fixed_param <- Vassign(surv_config$additional_param, weibPar)
  }
  
  random_effects <- sapply(glmeObjects, function(x){x$random_effects}) %>% unlist %>% unique
  fixed_param <- lapply(c(glmeObjects, list(survObject)), function(x){
    Vassign(x$fixed_param$names, x$fixed_param$start_values)
  }) %>% c(., fixed_param) %>% do.call(c, .)
  sigma_param <- lapply(long_config, function(x){x$additional_param}) %>% 
    unlist() %>% Vassign(rep(0.5, length(.)))
  disp_param <- lapply(glmeObjects, function(x){
    Vassign(x$disp_param$names, x$disp_param$start_values)
  })
  disp_param <- c(disp_param, sigma_param) %>% do.call(c, .)
  
  p <- length(fixed_param)  #  dimension of fixed parameters 
  q <- length(random_effects) # diemnsion of random effects
  group <- long.data[ , subject_id]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repreated measurements for each subject 
  N <- nrow(long.data)
  
  if(randeff_info$distribution == 'normal'){
    invSIGMA <- SIGMA <- diag(1, q, q)
  } else if(randeff_info$distribution == 't-dist'){
    SIGMA <- diag(1, q, q)*(randeff_info$degree - 2)/randeff_info$degree
    invSIGMA <- solve(SIGMA)
  }
  
  # add bounds for the parameters
  if(is.null(survObject$distribution)){
    fixed_lower <- rep(-Inf, p)
  }  else if(survObject$distribution == "weibull"){
    fixed_lower <- c(rep(-Inf, p-1), 0)
  }
  names(fixed_lower) <- names(fixed_param)
  disp_lower <- lapply(glmeObjects, function(x){x$disp_param$lower}) %>% 
    unlist %>% c(., rep(0, length(sigma_param)))
  disp_upper <- lapply(glmeObjects, function(x){x$disp_param$upper}) %>% 
    unlist %>% c(., rep(Inf, length(sigma_param)))
  names(disp_lower) <- names(disp_upper) <- names(disp_param)

  likDiff <- Diff <- 1
  convergence <- 1
  m <- 1
  Lval0 <- NULL
  
  while(likDiff > itertol & Diff > Ptol & m < iterMax){
    ################################################
    ########### Estimate random effects ############
    ################################################
    # estimate random effects by max log h-likelihood
    Bi <- est_randeff(RespLog, long.data, surv.data, subject_id,
                      random_effects, invSIGMA,
                      par_val = c(disp_param, fixed_param),
                      distribution = randeff_info$distribution,
                      degree = randeff_info$degree,
                      Silent = T, scale = T)
    # print(cov(Bi[, -1]))
    # generate MCMC samples of random effects (adaptive Metro-Hastings)
    par_list <- list(fixed_param = fixed_param,
                     SIGMA = solve(invSIGMA),
                     disp_param = disp_param,
                     RespLog = RespLog,
                     random_effects = random_effects,
                     distribution = randeff_info$distribution,
                     degree = randeff_info$degree)
    Bi_samples <- sample_randeff_MCMC(long.data, surv.data, subject_id,
                                      Bi, par_list, mcmc_size, burn_in, thin)

    ################################################
    ############# Estimate parameters ##############
    ################################################
    # Robust estimate of fixed parameters
    rob_fixed_param <- rob_fixed_param_sd <- c()
    rob_disp_param <- c()
    for(model_object in glmeObjects){
      if(!is.null(model_object$robust) & !is.null(model_object$left_censoring)){
        if(model_object$distribution == "normal"){
          param <- c(fixed_param, disp_param) %>%
            .[names(.) %in% c(model_object$fixed_param$names, model_object$disp_param$names)]
          sigma <- disp_param[grepl(model_object$response, names(disp_param), fixed = TRUE)]
          robest <- est_fixeff_rob(param = param,
                                   sigma = sigma,
                                   other_param = fixed_param[!names(fixed_param) %in% names(param)],
                                   model_object,
                                   long.data,
                                   Bis = Bi_samples, invSIGMA,
                                   Silent)
          rob_fixed_param <- c(rob_fixed_param, robest$gamma)
          rob_fixed_param_sd <- c(rob_fixed_param_sd, sqrt(-diag(solve(robest$jac))))
          
          # Robust estimate of dispersion parameters
          new_sigma <- est_sigma_rob(param = sigma,
                                     fixed_param = robest$gamma,
                                     model_object = model_object,
                                     long.data = long.data,
                                     Bis = Bi_samples, Silent = Silent)
          rob_disp_param <- c(rob_disp_param, new_sigma)
        }
      }
    }
    # cat("robust estimates done ... \n")
    
    # Non-robust estimate of fixed parameters
    param <- fixed_param[!names(fixed_param) %in% names(rob_fixed_param)]
    other_param <- c(rob_fixed_param, rob_disp_param,
                     disp_param[!names(disp_param) %in% names(rob_disp_param)])
    res_fixeff <-  est_fixeff(param = param,
                              other_param = other_param,
                              RespLog, long.data, surv.data,
                              Bi, invSIGMA,
                              lower = fixed_lower[names(fixed_lower) %in% names(param)],
                              distribution = randeff_info$distribution, 
                              degree = randeff_info$degree,
                              Silent)
    new_fixed_param <- c(rob_fixed_param, res_fixeff$gamma)

    # Non-robust estimate of dispersion parameters
    other_param <- c(new_fixed_param, rob_disp_param)
    param <- disp_param[!names(disp_param) %in% names(other_param)]
    res_disp <- est_disp(param, invSIGMA, other_param, RespLog,
                         deriv_vars = c(random_effects, names(fixed_param)),
                         long.data, surv.data, Bi, Lval0,
                         lower = disp_lower[names(disp_lower) %in% names(param)],
                         upper = disp_upper[names(disp_upper) %in% names(param)],
                         distribution = randeff_info$distribution,
                         degree = randeff_info$degree, 
                         Silent)
    new_disp_param <- c(rob_disp_param, res_disp$sigma)
    new_invSIGMA <- res_disp$invSIGMA
    new_Lval0 <- res_disp$Lval

    # Estimate baseline hazard function in Cox model
    if(is.null(survObject$distribution)){
      new_h <- est_base_hazard(surv.data, survObject, 
                               c(new_fixed_param, new_disp_param), Bi)
      surv.data <- new_h$surv.data
    }
    
    ####################################################    
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    new_loglike_value <- approx_log_like(RespLog, c(new_fixed_param, new_disp_param), 
                                         long.data, surv.data, Bi, new_invSIGMA, 
                                         distribution = randeff_info$distribution, 
                                         degree = randeff_info$degree)
    
    if(m == 1){
      likDiff <- 1
    } else{
      likDiff <- abs((new_loglike_value - loglike_value)/loglike_value)
    }
    
    # calcuate relative changes in mean parameters
    old_fixed_param <- c(fixed_param, disp_param)[names(new_fixed_param)]
    old_disp_param <- c(fixed_param, disp_param)[names(new_disp_param)]
    relative_change <- abs((unlist(new_fixed_param) - unlist(old_fixed_param))/(unlist(old_fixed_param) + 1e-16))
    Diff <- mean(relative_change)
    
    # print iterating result
    cat("############## Iteration:", m, "###############","\n")
    cat("Approximate log likelihood:", new_loglike_value, "\n")
    cat("Estimates of fixed parameters: \n")
    cat(round(unlist(new_fixed_param), 2), "\n")
    cat("Estimates of dispersion parameters:\n")
    cat(round(unlist(new_disp_param), 2), "\n")
    cat("Average relative changes in fixed parameters:",
        round(Diff, 3), "\n")
    cat("Relative change in log likelihood:", likDiff, "\n")
    cat("##########################################","\n")
    
    fixed_param <- new_fixed_param
    invSIGMA <- new_invSIGMA
    disp_param <- new_disp_param
    Lval0 <- new_Lval0
    loglike_value <- new_loglike_value
    m <- m + 1
  }
  
  # messages about convergence success or failure
  if((likDiff > itertol  & Diff > Ptol)){
    warning("Iteration limit reached without covergence.")
    convergence <- 1
  }
  if(likDiff <= itertol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= itertol.")
    convergence <- 0
  }
  if(Diff <= Ptol){
    message("Successful convergence. Iteration stops because FixedParDiff <= Ptol.")
    convergence <- 0
  }
  
  fixed_param_sd <- NULL
  if(convergence == 0){
    # Estimate s.e. of non-robust fixed parameter estimates
    nonrob_fixed_param_sd <- c()
    for(model_object in c(glmeObjects, list(survObject))){
      if(is.null(model_object$robust)){
        param <- fixed_param[model_object$fixed_param$names]
        other_param <- c(fixed_param, disp_param) %>%
          .[!names(.) %in% names(param)]
        sd_est <- get_sd_MCEM(param, RespLog, long.data, surv.data,
                              Bi_samples, random_effects, subject_id, 
                              other_param, 
                              distribution = randeff_info$distribution, 
                              degree = randeff_info$degree,
                              Silent)
        nonrob_fixed_param_sd <- c(nonrob_fixed_param_sd, sd_est)
      }
    }
    fixed_param_sd <- c(rob_fixed_param_sd, nonrob_fixed_param_sd)[names(fixed_param)]
  }
  
  output <- list(fixed_est = unlist(fixed_param), 
                 fixed_sd = fixed_param_sd,
                 disp_est = unlist(disp_param),
                 Bi = Bi,
                 Bi_samples = Bi_samples,
                 SIGMA = solve(invSIGMA),
                 convergence = convergence,
                 loglike_value = loglike_value,
                 long.data = long.data,
                 surv.data = surv.data,
                 RespLog = RespLog,
                 rob_param = names(rob_fixed_param))
  class(output) <- "hhjm"
}


