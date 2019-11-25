#' Derive H matrix in profile h-likelihood functions
#' 
#' Density of random effects is excluded here and will 
#' be included separately while estimating parameters.
#' Note that this function returns -H.
#' 
getHmat <- function(RespLog, pars){
  p <- length(pars)
  negH_long <-  getHessian(RespLog[[1]], pars)
  negH_surv <-  getHessian(RespLog[[2]], pars)
  list(negH_long = negH_long, negH_surv = negH_surv)
}

#' Assign names to a vector
#'
Vassign <- function(name, value){
  dat <- data.frame(as.list(value))
  names(dat) <- name
  return(dat)
}

#' Calculate Hessian matrix 
#' 
getHessian <- function(lik1, pars){
  lik <- parse(text = lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q^2))
  for(i in 1:q){
    for(j in 1:q){
      k <- (i - 1)*q + j
      result[[k]] <- D(D(lik, pars[i]), pars[j])
      names(result)[k] <- paste(pars[i], pars[j], sep = ",")
    }
  }
  return(result)
}

#' This function returns a matrix with elements in string.
#' Moreover, it returns L(l)'L(l) = SIGMA, which is a spherical
#' parameterization with diagonal elements equal to 1.
#' 
strMat <- function(q2){
  L <- c()
  L <- rbind(L, c(1, rep(0, q2 - 1)))
  Mpar <- c()
  
  for(i in 2:q2){
    l0 <- 1
    Li <- c()
    for(j in 1:i){
      m0 <- paste('L', i, 2:min((j+1), i), sep = "")
      Mpar <- c(Mpar, m0)
      
      if(j < i){
        sin0 <- rep(c("sin(", "cos("), c(length(m0) - 1, 1))
      } else {
        sin0 <- rep("sin(",  length(m0))
      }
      l1 <- paste(sin0, m0, rep(")", length(m0)), sep = "")
      Li <- c(Li, paste(c(l0, l1), collapse = "*"))
    }
    L <- rbind(L, c(Li, rep(0, q2 - i)))
  }
  L <- t(L)
  
  M <- matrix(NA, q2, q2)  # M=L'L
  for(i in 1:q2){
    for(j in 1:q2){
      M[i, j] <- Simplify(paste(L[, i], L[, j], sep = "*", collapse = "+"))
    }
  }
  diag(M) <- "1"
  
  Mpar <- unique(Mpar)
  Mp <- length(Mpar)
  dM <- array(NA, dim = c(q2, q2, Mp))
  for(i in 1:Mp){
    ttz <- unlist(lapply(M, function(x){ Vderiv(x, Mpar[i]) }))
    dM[,,i] <- matrix(as.character(ttz), q2, q2)
  }
  
  return(list(M=M, Mpar=Mpar, dM=dM))
}

#' Compute derivatives of an expression w.r.t 
#' a vector of parameters respectively
#' 
Vderiv <- function(lik1, pars){
  lik <- parse(text = lik1)
  q <- length(pars)
  result <- as.list(rep(NA, q))
  
  for(i in 1:q){
    result[[i]] <- D(lik, pars[i])
    names(result)[i] <- pars[i]
  }
  return(result)
}  

#' Evaluate a matrix
#' 
evalMat <- function(mat, q, data = NULL, par.val = NULL, 
                    raneff.val = NULL){
  D <- matrix(0, q, q)
  if(!is.data.frame(par.val)){
    par.val <- data.frame(as.list(par.val))
  }
  for(i in 1:q){
    for(j in 1:q){
      kk <- (i-1)*q+j
      Di <- with(data, with(par.val, with(raneff.val, eval(parse(text = mat[kk])))))
      if(!is.null(data)){
        if(nrow(data) > 1 & length(Di) == 1){ 
          D[i, j] <- Di*nrow(data) 
        } else { 
          D[i,j] <- sum(Di) 
        }
      } else {
        D[i,j] <- sum(Di)
      } 
    }
  }
  return(D)
}