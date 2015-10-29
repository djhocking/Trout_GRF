std <- function(data, var) {
  var_std <- (data[ , c(var)] - mean(data[ , c(var)], na.rm = TRUE)) / sd(data[ , c(var)], na.rm = TRUE)
}

std_vec <- function(x) {
  var_std <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  return(var_std)
}