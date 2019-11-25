#' Estimate fixed parameters by maximizing the profile h-likelihood function
#' 
#' @importFrom dplyr select left_join
#' @importFrom Deriv Simplify
est_fixeff <- function(param, other_param,
                       RespLog, long.data, surv.data,
                       Bi, invSIGMA,
                       lower = NULL,
                       distribution = "normal", degree = 0,
                       Silent = TRUE){ 
  
  subject_id <- names(Bi)[1]
  random_effects <- names(Bi)[-1]
  q <- length(random_effects)
  n <- nrow(Bi)
  p <- length(param)
  B <- dplyr::select(long.data, subject_id) %>%
    dplyr::left_join(Bi, by = subject_id)
  
  # Derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars = random_effects)
  deno <- diag(as.matrix(Bi[, -1])%*%invSIGMA%*%t(Bi[, -1]))
  numo <- as.matrix(Bi[, -1])%*%invSIGMA
  # evaluate the -H matrix for the random effects part
  if(distribution == "normal"){
    negH3 <- -invSIGMA*n
    lhlike3 <- -diag(0.5*as.matrix(Bi[, -1])%*%invSIGMA%*%t(Bi[, -1]))
  } else if (distribution == "t-dist"){
    negH3 <- lapply(1:n, function(i){
      numo1 <- (degree + q)*invSIGMA*(degree + deno[i]) - 
        (degree + q)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      -numo1/((degree + deno[i])^2)
    }) %>% Reduce('+', .)
  }

  # Negative adjusted profile h-likelihood with param as input
  ff <- function(xx){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)

    # evaluate the -H matrix
    negH1 <- evalMat(Hmats[[1]], q, long.data, par_val, raneff.val = B)   # longitudinal part
    negH2 <- evalMat(Hmats[[2]], q, surv.data, par_val, raneff.val = Bi)  # survival part
    H <-  -(negH1 + negH2 + negH3)

    # evaluate h-likelihood
    lhlike <- eval_log_hlike(RespLog, par_val, long.data, surv.data, 
                             B, Bi, invSIGMA, distribution, degree)
    # evaluate profile h-likelihood
    fy <- lhlike - 0.5*log(det(H/2/pi))
    return(-fy)
  }

  # Gradient of ff()
  gr <- function(xx){
    # all parameter values
    par_val <- Vassign(names(param), xx) %>% c(., other_param)

    # evaluate inverse -H matrix
    negH1 <- evalMat(Hmats[[1]], q, long.data, par_val, raneff.val = B)   # longitudinal part
    negH2 <- evalMat(Hmats[[2]], q, surv.data, par_val, raneff.val = Bi)  # survival part
    H <-  -as.matrix(negH1 + negH2 + negH3)
    invH <- solve(H)

    # calculate and evaluate dH/dbeta, where beta are fixed parameters
    lapply(1:p, function(i){
      # from longitudinal models
      dH <- lapply(1:(q^2), function(j){
        Simplify(Vderiv(Hmats[[1]][j], names(param)[i]))
      }) %>% unlist()
      # from survival model
      dH2 <- lapply(1:(q^2), function(j){
        Vderiv(Hmats[[2]][j], names(param)[i])
      }) %>% unlist()

      dH_val1 <- evalMat(dH, q, long.data, par_val, raneff.val = B)
      dH_val2 <- evalMat(dH2, q, surv.data, par_val, raneff.val = Bi)
      dH_val <- -(dH_val1 +dH_val2)  # Note: Hmats returns -H matrix

      # calculate dh/dbeta: derivative of log h-likelihood  
      df1 <- Vderiv(RespLog[[1]], names(param)[i])
      df1_val <- with(long.data, with(par_val, with(B, eval(parse(text = df1)))))
      df2 <- Vderiv(RespLog[[2]], names(param)[i])
      df2_val <- with(surv.data, with(par_val, with(Bi, eval(parse(text = df2)))))
      
      -(sum(df1_val) + sum(df2_val) - 0.5*matrix.trace(as.matrix(invH%*%dH_val)))
    }) %>% unlist
  }
  
  # Start iteration of finding optimal solution
  if(Silent == F) check <- 0  else check <- 1
  if(is.null(lower)){
    lower <- -Inf
    method <- "BFGS"
  } else {
    method <- "L-BFGS-B"
  }
  message <- -1
  M <- 0
  str_val <- param

  while(message != 0 & M < 20){
    str_val <- sapply(str_val, function(x){
      x + rnorm(1, 0, min(1, abs(x/5)))
    })
    result <- try(optim(par = str_val, fn = ff, gr = gr,
                        method = method, lower = lower,
                        control = list(trace = 1 - check, maxit = 1000)),
                  silent = T)
    error_mess <- attr(result, "class")
    
    if(length(error_mess) != 0 ){ 
      message <- -1
    } else {
      message <- result$convergence 
      str_val <- result$par
    }
    
    if(Silent == F){ print(message); print(M); print(result) }
    M <- M +1
  }
  
  if(message == 0){
    gamma <- result$par
    names(gamma) <- names(param)
    fval <- result$value
  } else {
    stop("Iteration stops because fixed parameters can not be successfully estimated.")  
  }
  
  return(list(gamma = gamma, fval = fval))
}