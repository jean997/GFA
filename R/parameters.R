#'@export
gfa_default_parameters <- function(){
  list(
  kmax = NULL,
  cond_num = 1000,
  max_lr_percent = 0.99,
  lr_zero_thresh = 1e-3,
  max_iter = 1000,
  extrapolate = TRUE,
  ebnm_fn_F = flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
  ebnm_fn_L = flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
  init_fn = flash_greedy_init_default,
  duplicate_check_thresh = 0.7,
  var_type = 2
  )

}
