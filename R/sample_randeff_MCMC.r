#' Log h-likelihood of longitudinal and survival data
likelihood <- function(xx, par_list, sub_long, sub_surv){
  par_val <- Vassign(par_list$random_effects, xx) %>% 
    c(., par_list$fixed_param, par_list$disp_param)
  # values from longitudinal data
  val_long <- with(par_val, with(sub_long, eval(parse(text = par_list$RespLog[[1]]))))
  # values from survival data
  val_surv <- with(par_val, with(sub_surv, eval(parse(text = par_list$RespLog[[2]]))))
  sum(val_long) + sum(val_surv)
}

#' Log prior of random effects
prior <- function(xx, par_list){
  invSIGMA <- par_list$SIGMA
  q <- length(par_list$random_effects)
  xx <- matrix(xx, nrow = 1)
  
  if(par_list$distribution == "normal"){
    val_raneff <-  -0.5*xx%*%invSIGMA%*%t(xx) + 0.5*log(det(invSIGMA)) - q*0.5*log(2*pi)
  } else if(par_list$distribution == "t-dist"){
    degree <- par_list$degree
    C_gamma <- gamma((degree + q)/2)/(gamma(degree/2)*degree^{q/2}*pi^{q/2})
    val_raneff <- -(degree + q)/2*log(1+xx%*%invSIGMA%*%t(xx)/degree) + 0.5*log(det(invSIGMA)) + log(C_gamma)
  }
  val_raneff
}

#' Un-normalized log posterior of random effects
posterior <- function(param, par_list, sub_long, sub_surv){
  likelihood(param, par_list, sub_long, sub_surv) + prior(param, par_list)
}

#' Generate MCMC samples using adaptive Metro_Hastings
#'
#' @importFrom MHadaptive Metro_Hastings mcmc_thin
#' @importFrom dplyr select left_join
sample_randeff_MCMC <- function(long.data, surv.data, subject_id,
                                Bi, par_list, mcmc_size, burn_in, thin){
  
  uniqueID <- unique(long.data[, subject_id])

  # The random effects are assumed independent across subjects.
  # Therefore, we can draw samples for each subject independently.
  # TODO: parallelize code to speed up
  Bi_samples <- lapply(uniqueID, function(id){
    sub_long <- subset(long.data, long.data[, subject_id] == id)
    sub_surv <- subset(surv.data, surv.data[, subject_id] == id)
    startvalue <- subset(Bi, Bi[, subject_id] == id) %>% 
      dplyr::select(-subject_id) %>% as.matrix()
    
    samples <- MHadaptive::Metro_Hastings(posterior, pars = startvalue, 
                                          prop_sigma = par_list$SIGMA, 
                                          par_names = par_list$random_effects,
                                          par_list = par_list,
                                          sub_long = sub_long, 
                                          sub_surv = sub_surv,
                                          iterations = mcmc_size*thin + burn_in - 1, 
                                          burn_in = burn_in, quiet = T) %>%
      MHadaptive::mcmc_thin(thin = thin) %>%
      .[["trace"]]
    output <- data.frame(samples)
    names(output) <- par_list$random_effects
    output
  })
  
  output <- lapply(Bi_samples, as.matrix) %>% simplify2array(., higher = T)
  res <- lapply(1:mcmc_size, function(k){
    Bi_k <- t(output[k, , ]) %>% scale %>% as.data.frame()
    df <- cbind(uniqueID, Bi_k)
    names(df)[1] <- subject_id
    df
  })
  names(res) <- paste0("sample", 1:mcmc_size)
  res
}

#################################################################
# Alternative implementation of MH sampling 
# (currently slower than Metro_Hastings and the results seem not consistent)

# proposalfunction <- function(param, SIGMA, par_list){
#   param <- as.numeric(param)
#   if(par_list$distribution == 'normal'){
#     mvtnorm::mvrnorm(1, mu = param, Sigma = SIGMA)
#   } else if(par_list$distribution == 't-dist'){
#     degree <- par_list$degree
#     mvtnorm::rmvt(1, delta = param, sigma = SIGMA*(degree-2)/degree,
#                 df = degree)
#   }
# }
# 
# run_metropolis_MCMC <- function(startvalue, SIGMA, par_list, 
#                                 sub_long, sub_surv,
#                                 mcmc_size, burn_in, thin){
#   samples <- array(dim = c(mcmc_size*thin + burn_in, length(startvalue)))
#   i <- 1
#   while(i <= nrow(samples)){
#     proposal <- proposalfunction(startvalue, SIGMA, par_list)
#     probab <- exp(posterior(proposal, par_list, sub_long, sub_surv) - 
#                    posterior(startvalue, par_list, sub_long, sub_surv))
#     if(!is.nan(probab)){
#       if(runif(1) < probab){
#         startvalue <- proposal
#         samples[i, ] <- proposal
#         i <- i+1
#       }
#     }
#     # cat("i = ", i, '\n')
#   }
# 
#   samples <- samples[seq(burn_in + 1, nrow(samples), by = thin), ]
#   samples <- data.frame(samples)
#   names(samples) <- par_list$random_effects
#   samples
# }
