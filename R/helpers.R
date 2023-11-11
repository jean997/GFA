
gfa_set_data <- function(Y, scale = NULL, S = NULL, R = NULL, params = NULL){
  if(!inherits(Y, "matrix")){
    stop("Data must have class matrix.\n")
  }
  if(!is.null(scale) & !is.null(S)) stop("Something wrong, this shouldn't happen.")
  dat <- list(n = nrow(Y), p = ncol(Y))
  R <- check_R(R, dat$p, params)
  if(!is.null(S)){
    if(is.null(R)){
      ## effect scale no R
      dat$S <- S
      dat$scale <- rep(1, dat$p)
    }else{
      ## effect scale yes R
      message("For now, with R and on effect scale,
              all I know how to do is convert to z-scores
              and store the scale parameter..\n")
      dat$scale <- get_scale_from_S(S)
      dat$Y <- Y/S
      dat$R <- R
    }
  }else{
    ## we are on z-score scale
    dat$Y <- Y
    if(is.null(scale)){
      dat$scale <- rep(1, dat$p)
    }else{
      stopifnot(inherits(scale, "numeric"))
      if(!length(scale) == dat$p){
        stop(paste0("Length of ", deparse(substitute(scale)), " does not match number of columns of data."))
      }
      dat$scale <- scale
    }
    if(is.null(R)){
      dat$S <- 1
    }else{
      dat$R <- R
    }
  }
  dat$params <- params
  return(dat)
}

## jean notes on scale:
## if Y is z-scores, then we estimate F on z-score scale
## and then convert back by multiplying each column of F
## by sqrt(N)
## The scale object will always contain the vector we
## need to multiply cols of F by to get to trait scale.
## If S is provided, find the scale by looking at ratios of sds
## sp/s1 \apprix sqrt(N_1)/sqrt(N_p). One option is to convert
## B_hat to z-scores and store the scale parametrs (1, sqrt(N2)/sqrt(N1), ...)

get_scale_from_S <- function(S){
  trait_scale <- apply(S[,-1], 2, function(s){
    median(s/S[,1])
  })
  trait_scale <- sqrt(c(1, scale))
  #variant_scale <- apply( t(t(S)/trait_scale)  , 1, median)
  return(trait_scale)
}

check_R <- function(R, p, params, tol = 1e-8){
  if(is.null(R)) return(R)
  if(!Matrix::isSymmetric(R)){
    stop(paste0(deparse(substitute(R)), " must be symmetric.\n"))
  }
  if(!ncol(R) == p){
    stop("Dimension of R not compatible with data.")
  }
  eMvals <- eigen(R, only.values = TRUE)$values
  if(!all(eMvals >= -1*tol)){
    stop(paste0(deparse(substitute(R)), " is not positive semi-definite.\n"))
  }
  if(!all(abs(diag(R)-1) < tol)){
    stop(paste0(deparse(substitute(R)), " should be a correlation matrix."))
  }

  vals <- eigen(R, only.values = TRUE)$values
  if(min(vals)/max(vals) < params$min_ev){
    stop("Ratio of smallest eigenvalue to largest eigenvalue is smaller than ", params$min_ev)
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
