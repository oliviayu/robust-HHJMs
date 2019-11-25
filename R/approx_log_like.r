#' Evaluate the approximate log marginal likelihood function
#' 
approx_log_like <- function(RespLog, par_val, long.data, surv.data, Bi,
                          invSIGMA, distribution, degree){
  par_val <- as.list(par_val)
  random_effects <- names(Bi)[-1]
  q <- length(random_effects)
  n <- nrow(Bi)
  subject_id <- names(Bi)[1]
  B <- dplyr::select(long.data, subject_id) %>% dplyr::left_join(Bi, by = subject_id)
  
  hloglike_value <- eval_log_hlike(RespLog, par_val, long.data, surv.data, 
                                   B, Bi, invSIGMA, distribution, degree)

  Hmats <- getHmat(RespLog, pars = random_effects)
  nH1 <- evalMat(Hmats[[1]], q, long.data, par.val = par_val, raneff.val = B)
  nH2 <- evalMat(Hmats[[2]], q, surv.data, par.val = par_val, raneff.val = Bi)  # survival part
  nH3 <- -invSIGMA*n
  Hval <- -(nH1 + nH2 + nH3)
  
  hloglike_value - 0.5*log(det(Hval/2/pi))
}

#' Evaluate log h-likelihood
#' 
eval_log_hlike <- function(RespLog, par_val, long.data, surv.data, B, Bi,
                           invSIGMA, distribution, degree){
  q <- ncol(Bi) - 1 
  lhlike1 <- with(long.data, with(par_val, with(B, eval(parse(text = RespLog[[1]])))))
  lhlike2 <- with(surv.data, with(par_val, with(Bi, eval(parse(text = RespLog[[2]])))))
  deno <- diag(as.matrix(Bi[, -1])%*%invSIGMA%*%t(as.matrix(Bi[, -1])))

  if(distribution == "normal"){
    lhlike3 <- 0.5*log(det(invSIGMA)) - diag(0.5*deno)
  }else if (distribution == "t-dist"){
    lhlike3 <- 0.5*log(det(invSIGMA)) - (degree + q)/2*log(1 + deno/degree)
  }
  sum(lhlike1) + sum(lhlike2) + sum(lhlike3)
}