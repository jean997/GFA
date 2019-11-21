#library(flashier)

#'@export
run_flashier <- function(mats, var_type  =  c("constant", "by_row", "by_col", "kronecker",
                                              "zero", "noisy_constant", "noisy_byrow",
                                              "noisy_bycol") ,
                         init_type = c("flashier", "soft_impute", "from_data")){

  var_type <- match.arg(var_type)
  init_type <- match.arg(init_type)
  mats$beta_hat[is.na(mats$se_hat)] <- NA
  args <- list(data = mats$beta_hat, output.lvl = 1, tol = 0.01)

  if(var_type=="constant"){
    var_args <- list(var.type=0)
  }else if(var_type == "by_row"){
    var_args <- list(var.type=1)
  }else if(var_type == "by_col"){
    var_args <- list(var.type=2)
  }else if(var_type == "kronecker"){
    var_args <- list(var.type=c(1, 2))
  }else if(var_type == "zero"){
    var_args <- list(S = mats$se_hat, var.type=NULL)
  }else if(var_type=="noisy_constant"){
    var_args <- list(S = mats$se_hat, var.type=0)
  }else if(var_type == "noisy_byrow"){
    var_args <- list(S = mats$se_hat, var.type=1)
  }else if(var_type == "noisy_bycol"){
    var_args <- list(S = mats$se_hat, var.type=2)
  }

  if(init_type == "flashier"){
    init_args <- list(fit = "full")
  }else if(init_type == "soft_impute"){
    init_args <- list(fit = "full", init.fn = init.fn.softImpute)
  }else if(init_type == "from_data"){
    si_res <- softImpute::softImpute(mats$beta_hat, rank.max = 30, lambda = 0, type = "als")
    init_args <- list(fit = "backfit.only", EF.init = si_res)
  }

  all_args <- c(args, var_args, init_args)
  res <- do.call(flashier, all_args)
  res$types <- c(init_type, var_type)
  return(res)

}

