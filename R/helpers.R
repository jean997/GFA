
gfa_set_data <- function(Y, scale, R = NULL, params = NULL){
  if(!inherits(Y, "matrix")){
    stop("Data must have class matrix.\n")
  }

  dat <- list(n = nrow(Y), p = ncol(Y))
  R <- check_R(R, dat$p, params)

  stopifnot(inherits(scale, "numeric"))
  if(!length(scale) == dat$p){
      stop(paste0("Length of ", deparse(substitute(scale)), " does not match number of columns of data."))
  }

  dat$Y <- Y
  dat$scale <- scale
  dat$S <- 1

  if(!is.null(R)){
    dat$R <- R
  }
  dat$params <- params
  return(dat)
}

## jean notes on scale:
## if Y is z-scores, then we estimate F on z-score scale
## and then convert back by dividing each column of F by sqrt(N)
## For binary traits, we need to divide each column of F by
## sqrt(N*P*(1-P))/C(K) where P = N_case/N and
## C(K) = K*(1-K)/dnorm(qnorm(K, lower.tail = F))
## and K = population prevalence.
binary_const <- function(N, N_case, pop_prev){
  if(any(is.na(N))){
    stop("Illegal missing values in N")
  }
  if(any(is.na(N_case))){
    stop("Illegal missing values in N_case")
  }
  if(any(is.na(pop_prev))){
    stop("Illegal missing values in pop_prev")
  }
  P <- N_case/N
  if(!all(P < 1)){
    stop("Some of N_case are bigger than N\n")
  }
  const <- sqrt(N*P*(1-P))*dnorm(qnorm(pop_prev, lower.tail = TRUE))/(pop_prev*(1-pop_prev))
  return(const)
}


## The scale object will always contain the vector we
## need to divide cols of F by to get to standardized trait scale/liability scale.
## If S is provided, find the scale by looking at ratios of sds
## sp/s1 \approx sqrt(N_1)sqrt(Var Y_p)/sqrt(N_p)*sqrt(Var Y_1) for continuous trait or
## sp/s1 \approx sqrt(N_1)/sqrt(N_p*P_p*(1-P_p)) for p binary 1 cont
## sp/s1 \qpprox sqrt(N_)
## We fit with z-scores and use scale
## (1, sqrt(N2)/sqrt(N1), ...) \approx (1, median(s1)/median(s2), ...)
## Adjustment for binary traits is done in gfa_fit
get_scale_from_S <- function(S){
  if(ncol(S) == 1) return(1)
  scale <- apply(S[,-1, drop = FALSE], 2, function(s){
    median(S[,1]/s) ## relative to sample size
  })
  scale <- c(1, scale)
  #variant_scale <- apply( t(t(S)/trait_scale)  , 1, median)
  return(scale)
}

check_R <- function(R, p, params, tol = 1e-8){
  if(is.null(R)) return(R)
  if(!Matrix::isSymmetric(R)){
    stop(paste0(deparse(substitute(R)), " must be symmetric.\n"))
  }
  if(!ncol(R) == p){
    stop("Dimension of R not compatible with data.")
  }
  if(!all(abs(diag(R)-1) < tol)){
    #stop(paste0(deparse(substitute(R)), " should be a correlation matrix."))
    warning("R has values different from 1 on the diagonal. This is ok.")
  }

  vals <- eigen(R, only.values = TRUE)$values
  if(!all(vals >= -1*tol)){
    stop(paste0(deparse(substitute(R)), " is not positive definite. Project
                to the nearest wellconditioned positive definite matrix using
                condition(R)\n"))
  }
  if(max(vals)/min(vals) > params$cond_num){

    stop( paste0(deparse(substitute(R)), " is illconditioned with condition number", max(vals)/min(vals),
                 " and maximum allowable condition number ", params$cond_num,
                 ". Either project to nearest well conditioned matrix using Matrix::nearPD(.., posd.tol = 1/condnum) or choose a larger maximum condition
         number using params = list(cond_num = x)."))
  }

  v <- sum((vals - min(vals))^2)/sum(vals^2)
  if(v < params$lr_zero_thresh){
    message("R appears to be the identity or very close to the identity, so I am leaving it out.")
    return(NULL)
  }

  return(R)
}

check_args_znbs <- function(is_null_Z, is_null_N, is_null_B, is_null_S){
  if(!sum(is_null_Z, is_null_B) == 1){
    stop("Please supply exactly one of Z_hat and B_hat")
  }
  if(!is_null_B & is_null_S){
    stop("If B_hat is supplied, S must also be supplied.")
  }
  if(!is_null_Z & !is_null_S){
    stop("If Z_hat is supplied do not supply S")
  }
  if(!is_null_B & !is_null_N){
    stop("N is not needed if B_hat is supplied")
  }
  if(!is_null_Z & is_null_N){
    warning("Assuming equal sample sizes (factors will be returned on z-score scale.")
  }
}
