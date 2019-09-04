calc_lfsr <- function(flash_fit) {
  return(lapply(1:get.dim(flash),
                function(n) sapply(1:get.n.factors(flash),
                                   function(k) lfsr.one.n(flash, k, n))))

}

lfsr_one_n <- function(flash_fit, k, n) {
  factor <- extract.factor(flash, k)
  if (is.zero(factor) || all.fixed(factor, n)) {
    lfsr <- rep(NA, get.dims(flash)[n])
  } else {
    ebnm.res <- solve.ebnm(factor, n, flash, output = "lfsr")
    if (!is.null(ebnm.res$posterior) && !is.null(ebnm.res$posterior$lfsr)) {
      lfsr <- ebnm.res$posterior$lfsr
      fix.dim <- get.fix.dim(factor)
      if (!is.null(fix.dim) && (fix.dim == n))
        lfsr[get.fix.idx(factor)] <- NA
    } else {
      lfsr <- NULL
    }
  }
  return(lfsr)
}
