#'@title Relative root mean squared error
#'@export
rrmse <- function(Bhat, B){
  sqrt(sum((Bhat-B)^2)/sum(B^2))
}
