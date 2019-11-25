#' Summarize fitting results
#' 
#' @export
summary.hhjm <- function(output, newSD = NULL, digits = 3){
  est <- output$fixed_est
  if(is.null(newSD)){
    sd <- output$fixed_sd
  } else {
    sd <- newSD
  }
  Zvalue <- est/sd
  Pvalue <- (1 - pnorm(abs(Zvalue), 0, 1))*2
  Coeff <- data.frame(Estimate = est, 
                      Std.Error = sd, 
                      Zvalue, Pvalue)
  print(round(Coeff, digits = digits))
}