#' Get log density function of propotional hazard model
#' 
get_log_density_ph <- function(model_object){
  if(is.null(model_object$distribution)){
    # Cox PH model
    log_density <- paste(model_object$event, "*( log(hazard0) +",
                         model_object$reg_equation, ") - exp(", 
                         model_object$reg_equation, ")*cum_hazard0")
    
    additional_param <- c("hazard0", "cum_hazard0")
    
  } else if(model_object$distribution == "weibull"){
    # Weibull PH model
    log_density <- paste(model_object$event, "*( Wlogscale+log(Wshape)+log(",
                         model_object$response, ")*(Wshape-1) +", 
                         model_object$reg_equation, ") - exp(Wlogscale+", 
                         model_object$reg_equation, ")*",
                         model_object$response, "^Wshape")
    
    additional_param <- c("Wlogscale", "Wshape")
  }
  
  list(log_density = log_density,
       additional_param = additional_param)
}

