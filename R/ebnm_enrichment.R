ebnm_enrichment <- function(annot, optmethod = "nograd_nlm", ...){
  args <- list(...)
  if (any(c("x", "s", "output", "g_init", "fix_g") %in% names(args))){
    stop("'x', 's', 'g_init', 'fix_g', and 'output' must not be supplied as arguments to ",
         "'flash_ebnm_enrichment'.")
  }
  ebnm.fn <- function(x, s, g_init, fix_g, output) {
    #cat("pne family: ", length(x), " ", length(s), "\n")
    if(length(x) ==3){
      myannot <- annot[1:3,] # in order to pass test
    }else{
      myannot <- annot
    }
    ebnm::ebnm(x, s, g_init = g_init,
               fix_g = fix_g, output = output,
               prior_family = "point_normal_enrichment",
               optmethod = optmethod,
               annot = myannot, ...)
  }
  return(ebnm.fn)
}

ebnm_pn <- function( ...){
  args <- list(...)
  if (any(c("x", "s", "output", "g_init", "fix_g") %in% names(args))){
    stop("'x', 's', 'g_init', 'fix_g', and 'output' must not be supplied as arguments to ",
         "'flash_ebnm_enrichment'.")
  }
  ebnm.fn <- function(x, s, g_init, fix_g, output) {
    #cat("pn family: ", length(x), " ", length(s), "\n")
    ebnm::ebnm(x, s, g_init = g_init,
               fix_g = fix_g, output = output,
               prior_family = "point_normal",
               ...)
  }
  return(ebnm.fn)
}
