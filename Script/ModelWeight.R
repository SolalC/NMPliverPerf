calculate_model_weights <- function(fits, n_samples, n_param) {
  # bic <- sapply(fits, function(f) {
  #   n_param * log(n_samples) - 2 * f$loglik
  # })
  bic <- c()
  for(i in 1:length(fits)){
    bicCurrent <- n_param[i] * log(n_samples) - 2 * as.numeric(fits[[i]]$loglik)
    bic <- c(bic, bicCurrent)
  }
  names(bic) <- names(fits)
  min_bic <- min(bic)
  delta_bic <- bic - min_bic
  
  bicw <- exp(-0.5 * delta_bic) / sum(exp(-0.5 * delta_bic))
  names(bicw) <- names(fits)
  return(list(BIC = bic, BICW = bicw))
}