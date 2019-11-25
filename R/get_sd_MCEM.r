#' Estimate s.e. of non-robust fixed parameters,
#' using MCEM approximation
#' 
get_sd_MCEM <- function(param, RespLog, long.data, surv.data,
                        Bis, random_effects, subject_id, other_param,
                        distribution, degree, Silent){
  
  q <- length(random_effects)
  p <- length(param)
  n <- nrow(surv.data)
  gr.long <- deriv(formula(paste("~", RespLog[[1]])), names(param), hessian = T)
  gr.surv <- deriv(formula(paste("~", RespLog[[2]])), names(param), hessian = T)
  
  gr <- function(xx, subB_i){
    fy <- numeric(p)
    # assign values to parameters
    par_val <- c(Vassign(names(param), xx), other_param)
    subB <- dplyr::select(long.data, subject_id) %>%
      dplyr::left_join(subB_i, by = subject_id)
    fy1 <- with(long.data, with(par_val, with(subB, attr(eval(gr.long), "gradient"))))
    fy2 <- with(surv.data, with(par_val, with(subB_i, attr(eval(gr.surv), "gradient"))))
    score_i <- apply(fy1, 2, function(x){
      tapply(x, long.data[, subject_id], sum)
    }) + fy2
    
    jac1 <- with(long.data, with(par_val, with(subB, attr(eval(gr.long), "hessian"))))
    jac2 <- with(surv.data, with(par_val, with(subB_i, attr(eval(gr.surv), "hessian"))))
    jac_i <- apply(jac1, MARGIN = c(2, 3), sum) + apply(jac2, MARGIN = c(2, 3), sum)
    
    # score_i%*%score_i
    return(list(score = score_i, jac = jac_i))
  }
  
  mcmc_size <- length(Bis)
  mat1 <- mat2 <- matrix(0, p, p)
  score <- rep(0, p)
  for(i in 1:mcmc_size){
    out <- gr(unlist(param), Bis[[i]])
    mat1 <- mat1 + t(out$score)%*%out$score
    mat2 <- mat2 + out$jac
    score <- score + out$score
  }
  
  result <- (mat1 + mat2 - t(score)%*%score/mcmc_size)/mcmc_size
  sqrt(diag(-solve(result)))
}