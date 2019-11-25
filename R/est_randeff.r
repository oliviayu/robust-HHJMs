#' Estimate random effects by maximizing h-likelihood function,
#' given fixed and dispersion parameters
#' 
est_randeff_by_hlike <- function(RespLog, sub_long, sub_surv, random_effects, 
                                 invSIGMA, par_val, distribution,
                                 degree, Silent = T){
  
  # Negative h-likelihood with random effects as input
  ff <- function(xx){
    fy <- numeric(1)
    Bi <- data.frame(id = NA, Vassign(random_effects, xx))
    B <- do.call("rbind", replicate(nrow(sub_long), Bi, simplify = FALSE))
    fy <- eval_log_hlike(RespLog, par_val, sub_long, sub_surv, B, Bi,
                         invSIGMA, distribution, degree)
    -fy
  }
  
  # Gradient of negative h-likelihood w.r.t. random effects
  q <- length(random_effects)
  gr_long <- deriv(formula(paste("~", RespLog[[1]])), random_effects)
  gr_surv <- deriv(formula(paste("~", RespLog[[2]])), random_effects)
  gr <- function(xx){
    fy <- numeric(q)
    # all parameter values
    par_val <- Vassign(random_effects, xx) %>% 
      c(., par_val)
    # gradient values from longitudinal and survival data
    val_long <-  with(par_val, with(sub_long, attr(eval(gr_long), "gradient"))) %>%
      apply(., 2, sum) %>% as.vector()
    #cat("long gradient:", val_long, '\n')
    val_surv <-  with(par_val, with(sub_surv, attr(eval(gr_surv), "gradient"))) %>%
      as.vector()
    #cat("surv gradient:", val_surv,'\n')
    # gradient values from the distribution of random effect
    if(distribution == "normal"){
      val_raneff <-  as.vector(-invSIGMA%*%xx)
    } else if(distribution == "t-dist"){
      denom = as.vector(degree+xx%*%invSIGMA%*%xx)
      val_raneff <- as.vector(-(degree + q)*(invSIGMA%*%xx)/denom)
    }
    fy <- -(val_long + val_surv + val_raneff)
    return(fy)
  }
  
  # Start iteration of finding optimal solution
  message <- -1
  M <- 0
  while(message != 0 & M < 50){
    bval <- rnorm(q, 0, 1)
    result <- try(optim(par = bval, fn = ff, gr = gr, method = "L-BFGS-B",
                        control = list(maxit = 2000, trace=0)), silent = T)
    error_mess <- attr(result, "class")
    
    if(length(error_mess) == 1){
      message <- -1
    } else {
      message <- result$convergence
    }
    if(Silent == FALSE){
      cat(paste("M = ", M, ", message = ", message, ".", sep = ""), '\n')
    }
    M <- M + 1
  }
  
  if(message == 0){
    result$par
  } else {
    stop("Iteration stops because random effects can not be successfully estimated.")  
  }
}

#' Obtain estimates of random effects 
#' 
est_randeff <- function(RespLog, long.data, surv.data, subject_id,
                        random_effects, invSIGMA, par_val,
                        distribution = "normal", degree = 0,
                        Silent = T, scale = T){
  
  uniqueID <- unique(long.data[, subject_id])
  n <- length(uniqueID)
  ni <- table(long.data[, subject_id])
  N <- nrow(long.data)
  q <- length(random_effects)
  
  Bi_est <- lapply(uniqueID, function(id){
    sub_long <- subset(long.data, long.data[, subject_id] == id)
    sub_surv <- subset(surv.data, surv.data[, subject_id] == id)
    bi <- est_randeff_by_hlike(RespLog, sub_long, sub_surv, random_effects,
                               invSIGMA, par_val, distribution,
                               degree, Silent)
  }) %>% do.call(rbind, .)

  if(scale == TRUE){ Bi_est <- scale(Bi_est) }
  
  Bi_df <- cbind(uniqueID, data.frame(Bi_est))
  names(Bi_df) <- c(subject_id, random_effects)
  # print(str(Bi_df))
  Bi_df
}