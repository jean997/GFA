
#'@export
generate_F_simple <- function(nblocks, type=c("nested", "difference",
                                               "checkers1", "checkers2")){

  type <- match.arg(type)
  if(type == "nested"){
    m <-  cbind(c(1,1,1),
                c(1, 1, 0),
                c(0, 0, 1))
  }else if(type == "difference"){
    m <- cbind(c(1, 1, 1), c(1, -1, 0), c(0, 0, 1))
  }else if(type == "checkers1"){
    m <- cbind(c(1, 1, 1),
               c(1, 1, 0 ),
               c(0, 1, 1))
  }else if(type = "checkers2"){
    m <-  cbind(c(1, 0, 1),
                c(1, 1, 0 ),
                c(0, 1, 1))
  }
  F_mat <- kronecker(diag(nblocks), m)
  return(F_mat)
}
