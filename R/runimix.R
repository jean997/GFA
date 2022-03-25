runimix <- function(n, g){
  stopifnot(class(g)=="unimix")
  K <- length(g$pi)
  Z <- sample(1:K, size=n, replace = TRUE, prob=g$pi)
  beta <- rep(NA, n)
  for(k in 1:K){
    nk <- sum(Z==k)
    if(nk==0) next
    beta[Z==k] <- runif(n=nk , min=g$a[k], max = g$b[k])
  }
  return(beta)
}


rnormalmix <- function(n, g){
  stopifnot(class(g)=="normalmix")
  K <- length(g$pi)
  Z <- sample(1:K, size=n, replace = TRUE, prob=g$pi)
  beta <- rep(NA, n)
  for(k in 1:K){
    nk <- sum(Z==k)
    if(nk==0) next
    beta[Z==k] <- rnorm(n=nk , mean=g$mean[k], sd = g$sd[k])
  }
  return(beta)
}
