
#' Get log density function of (generalized) linear model
#' 
get_log_density_glm <- function(model_object){
  additional_param <- NULL
  left_censored <- !is.null(model_object$left_censoring)

  if(!left_censored){

    if(model_object$distribution == "binomial"){
      
      res <- paste(model_object$response, "*(", 
                   model_object$reg_equation, ")-", "log(1+exp(", 
                   model_object$reg_equation, "))")
      
    } else if(model_object$distribution == "normal"){
      
      sigma <- paste0(model_object$response, "_sigma")
      res <- paste("- 0.5*(", model_object$response, "- (", 
                   model_object$reg_equation, "))^2/", 
                   sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
      additional_param <- c(additional_param, sigma)

    } else if(model_object$distribution == "poisson"){
      
      res <- paste(model_object$response, "*(", 
                   model_object$reg_equation, ")-exp(",
                   model_object$reg_equation, ")-log(factorial(",
                   model_object$response, "))")
    }
    
  } else {

    sigma <- paste0(model_object$response, "_sigma")
    log_density <- paste("- 0.5*(", model_object$response, 
          "- (", model_object$reg_equation, "))^2/", 
          sigma, "^2-log(", sigma, ")-0.5*log(2*pi)")
    Cresp <- model_object$left_censoring$indicator
    limit_val <- model_object$left_censoring$limit_value
    CstdNmu <- paste("(", limit_val, "-(",
                     model_object$reg_equation, "))/", sigma)
    leftint <- paste("1/(1+exp(-1.702*", CstdNmu, "))")

    if(model_object$left_censoring$method == "tobit"){
      # A Tobit model, assuming left-censored data follow normal distribution
      res <- paste( "(", log_density, ")*(1-", Cresp, ") +log(", leftint, ")*",Cresp )

    } else if(model_object$left_censoring$method == "truncated"){
      # Model uncensored data only, assumed following truncted normal distribution
      res <- paste("(", log_density, "-log(1-", leftint, "))*(1-", Cresp, ")")
    }
    
    additional_param <- c(additional_param, sigma)
  }
  
  list(log_density = paste("(", res, ")"),
       additional_param = additional_param)
}
