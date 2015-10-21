rmse <- function(error, na.rm = T) {
  sqrt(mean(error^2, na.rm = T))
}