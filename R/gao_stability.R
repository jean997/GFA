
#'@export
gao_stability_sparse <- function(L1, L2, rows_match=TRUE){
  if(!rows_match){
    L1 <- t(L1)
    L2 <- t(L2)
  }

  L1 <- L1[,colSums(L1) != 0]
  L2 <- L2[,colSums(L2) != 0]

  k1 <- ncol(L1)
  k2 <- ncol(L2)
  sigma <- abs(cor(L1, L2))

  #cat(dim(sigma), "\n")
  part1 <- apply(sigma, 1, function(x){
    max(x) - sum((((x > mean(x)))*x))/(k2-1)
  })
  part1 <- sum(part1)/(2*k1)

  part2 <- apply(sigma, 2, function(x){
    max(x) - sum((((x > mean(x)))*x))/(k1-1)
  })
  part2 <- sum(part2)/(2*k2)
  return(part1 + part2)
}
