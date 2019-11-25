#' Evaluate H matrix
#'
eval_Hmat_for_disp <- function(Hmats, par_val, long.data, surv.data, 
                      B, Bi, invSIGMA, distribution, degree){
  
  q1 <- sqrt(length(Hmats$negH_long))
  q2 <- ncol(Bi) - 1
  n <- nrow(Bi)

  # part 1: longitudinal data
  nD1 <- evalMat(Hmats[[1]], q1, long.data, par_val, raneff.val = B)
  # part 2: survival data
  nD2 <- evalMat(Hmats[[2]], q1, surv.data, par_val, raneff.val = Bi)
  # part 3: random effect
  if(distribution == 'normal'){
    nD3 <- bdiag( -invSIGMA*n, diag(0, q1 - q2, q1 - q2))
  } else if(distribution == 't-dist'){
    deno <- diag(as.matrix(Bi[, -1])%*%invSIGMA%*%t(Bi[, -1]))
    numo <- as.matrix(Bi[, -1])%*%invSIGMA
    nD3 <- matrix(0, nrow = q2, ncol = q2)
    for(i in 1:n){
      numo1 <- (degree + q2)*invSIGMA*(degree + deno[i]) - (degree + q2)*as.matrix(numo[i, ])%*%(2*numo[i, ])
      nD3_i <- -numo1/((degree + deno[i])^2)
      nD3 <- nD3 + nD3_i
    }
    nD3 <- bdiag(nD3, diag(0, q1 - q2, q1 - q2))
  }
  -as.matrix(nD1+nD2+nD3)
}

#' Estimates the dispersion parameters in the joint models, 
#' by maximizing the adjusted profile h-likelihood function.
#'
est_disp <- function(param, invSIGMA, other_param,
                     RespLog, deriv_vars,
                     long.data, surv.data,
                     Bi, Lval0 = NULL, lower, upper,
                     distribution = 'normal', degree = 0,
                     Silent = T){ 
  
  q <- length(deriv_vars)
  q1 <- length(param)
  q2 <- ncol(Bi) - 1
  n <- nrow(Bi)
  L2  <- strMat(q2)   # Return matrices of strings for cov(bi) 
  q0 <- length(L2$Mpar)
  subject_id <- names(Bi)[1]
  B <- dplyr::select(long.data, subject_id) %>% dplyr::left_join(Bi, by = subject_id)
  
  # derive -H matrix, where H is defined in the adjusted profile h-likelihood
  Hmats <- getHmat(RespLog, pars = deriv_vars)
  # derive derivative of -H
  dhlike1 <- deriv(formula(paste("~", RespLog[[1]])), names(param))
  dhlike2 <- deriv(formula(paste("~", RespLog[[2]])), names(param))
  dH1 <- dH2 <- as.list(rep(NA, q^2))
  for(i in 1:q^2){
    dH1[[i]] <- deriv(formula(paste("~", Hmats[[1]][i])), names(param))
    dH2[[i]] <- deriv(formula(paste("~", Hmats[[2]][i])), names(param))
  }
  
  # Return the negative value of the adjusted profile h-likelihood
  ff <- function(xx){
    if(!is.null(param)){
      # The first q1 are dispersion parameters in the models,
      # and the rests are for the cov. matrix of random effects
      par_val <- c(Vassign(names(param), xx[1:q1]), other_param)
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    } else { 
      # No other dispersion parameters except for the covariance 
      # matrix of random effects
      par_val <- other_param
      Lval <- Vassign(L2$Mpar, xx)
    }
    if(distribution == 'normal'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval) # i.e. SIGMA
      mat <- solve(invmat)
    } else if(distribution == 't-dist'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval)*(degree - 2)/degree
      mat <- solve(invmat)
    }
    # evaluate the H matrix
    H <- eval_Hmat_for_disp(Hmats, par_val, long.data, surv.data, B, Bi, 
                   mat, distribution, degree)
    # evaluate the h-likelihood
    lhlike <- eval_log_hlike(RespLog, par_val, long.data, surv.data, B, Bi,
                             mat, distribution, degree)
    # evaluate the negative adjusted profile h-likelihood
    -(lhlike - 0.5*log(det(H/2/pi)))
  }
  
  # Return the gradient of ff
  gr <- function(xx){
    # assign values to the parameters
    if(!is.null(param)){
      # The first q1 are dispersion parameters in the models,
      # and the rests are for the cov. matrix of random effects
      par_val <- c(Vassign(names(param), xx[1:q1]), other_param)
      Lval <- Vassign(L2$Mpar, xx[-(1:q1)])
    } else { 
      # No other dispersion parameters except for the covariance 
      # matrix of random effects
      par_val <- other_param
      Lval <- Vassign(L2$Mpar, xx) 
    }
    if(distribution == 'normal'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval)
    } else if(distribution=='t-dist'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval)*(degree-2)/degree
    }
    mat <- solve(invmat)

    # evaluate the H matrix
    H <- eval_Hmat_for_disp(Hmats, par_val, long.data, surv.data, B, Bi, 
                   mat, distribution, degree)
    invH <- solve(H)
    
    # derivative w.r.t. dispersion parameters other than cov(bi)
    dh_val1 <- with(par_val, with(long.data, with(B, attr(eval(dhlike1), "gradient")))) %>%
      apply(., 2, sum)
    dh_val2 <- with(par_val, with(surv.data, with(Bi, attr(eval(dhlike2), "gradient")))) %>%
      apply(., 2, sum)
    myDmat <- myDmat2 <- matrix(NA, nrow = q^2, ncol = q1)
    for(i in 1:(q^2)){
      myDmat[i, ] <- with(par_val, with(long.data, with(B, attr(eval(dH1[[i]]), "gradient")))) %>%
        apply(., 2, sum)
      myDmat2[i, ] <- with(par_val, with(surv.data, with(Bi, attr(eval(dH2[[i]]), "gradient")))) %>%
        apply(., 2, sum)
    }
    traces <- sapply(1:q1, function(i){
      dH1_val <- matrix(myDmat[, i], q, q)
      dH2_val <- matrix(myDmat2[, i], q, q)
      matrix.trace(invH%*%(-dH1_val-dH2_val))
    })
    fy1 <- -(dh_val1 + dh_val2 - 0.5*traces)
    
    # derivative w.r.t. cov(bi)
    fy2 <- sapply(1:q0, function(i){
      dM <- L2$dM[, , i]
      dM_val <- evalMat(dM, q2, par.val = Lval)
      Dmat <- -mat%*%dM_val%*%mat # d invSIGMA/d xi
      XDmat <- diag(as.matrix(Bi[, -1])%*%Dmat%*%t(as.matrix(Bi[, -1])))
      if(distribution == "normal"){
        # from l_h
        gh1 <- sum(-0.5 * matrix.trace(mat%*%dM_val) - 0.5*XDmat)
        # from correction term log(|H|)
        DH <- Dmat*n
      } else if (distribution == "t-dist"){
        deno <- diag(as.matrix(Bi[, -1])%*%mat%*%t(Bi[, -1]))
        numo <- as.matrix(Bi[, -1])%*%mat

        gh1 <- sum( 0.5* matrix.trace(invmat%*%Dmat) - (degree + q2)/(degree + deno)*0.5*XDmat)
        DH <- matrix(0, nrow=q2, ncol=q2)
        for(j in 1:n){
          XBi <- t(as.matrix(Bi[j, -1]))%*%as.matrix(Bi[j, -1])
          fg1 <- (degree + q2)*mat*(degree + deno[j]) - (degree + q2)*as.matrix(numo[j, ])%*%(2*numo[j, ]) # == numo1
          fg2 <- (degree+deno[j])^2
          Dg1 <- (degree*(degree + q2)+(degree + q2)*deno[j])*Dmat + (degree + q2)*mat*XDmat[j]-
            2*(degree + q2)*Dmat%*%XBi%*%mat - 2*(degree + q2)*mat%*%XBi%*%Dmat
          Dg2 <- as.numeric(2*(degree + deno[j])*XDmat[j])
          DH_j <- -(Dg1*fg2 - fg1*Dg2)/(fg2^2)
          DH <- DH + DH_j
        }
      }
      DH2 <- -bdiag(DH, diag(rep(0, q - q2)))
      gh2 <- -0.5*matrix.trace(as.matrix(invH%*%DH2))
      -(gh1+gh2)
    })
    
    c(fy1, fy2)
  }
  
  # Start iteration of finding optimal solution
  message <- -1
  M <- 0
  if(Silent==F) check=0  else check=1
  if((q1+q2) == 1){
    method1 <- method2 <- "Brent"
  } else{
    method <- "L-BFGS-B"
  }
  while(message != 0 & M < 10){
    if(is.null(Lval0)){
      Lval00 <- runif( q0, 0, pi)
    } else{
      Lval00 <- Lval0 + rnorm(length(Lval0), 0, 0.01)
    }
    str_val0 <-  c(unlist(param) + rnorm(length(param), 0, 0.3), Lval00)
    result <- try(optim(par = str_val0, fn = ff, gr = gr, method = method,
                        lower = c(lower, rep(0, q0)), upper = c(upper, rep(pi, q0)),
                        control = list(trace = 1 - check, maxit = 2000)), silent = T)
    error_mess <- attr(result, "class")
    
    if(length(error_mess) != 0 ){
      message <- -1
    } else {
      message <- result$convergence
      str_val <- result$par
    }
    
    if(Silent == F){ print(message); print(M); print(result) }
    M <- M + 1
  }
  
  if(message == 0){
    if(q1 > 0){
      output <- result$par[1:q1]
      names(output) <- names(param)
      Lval <- Vassign(L2$Mpar, result$par[-(1:q1)])
    } else {
      output <- NULL
      Lval <- Vassign(L2$Mpar, result$par)
    }
    if(distribution == 'normal'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval)
    } else if (distribution == 't-dist'){
      invmat <- evalMat(as.list(L2$M), q2, par.val = Lval)*(degree - 2)/degree
    }
    mat <- solve(invmat)
  } else {
    stop("Iteration stops because dispersion parameters can not be successfully estimated.
         The Hessian matrix of random effects and fixed parameters may not be invertible.
         It's probably because the random effects are highly correlated with each other.")
  }
  
  list(sigma = output, invSIGMA = mat, Lval = Lval)
}


