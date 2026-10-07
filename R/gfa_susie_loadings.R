# Extract per-variant PIPs and posterior means of the loadings L (variant level,
# original LD basis) from a gfa_fit result that used susie_ebnm_fn for the
# loadings prior.  Requires susie_ebnm_fn.R (for susie_state_summary).
#
# Where things live
#   gfa_res$fit            the flashier "flash" object
#   gfa_res$fit$L_ghat[[k]] the fitted loadings prior for factor k = the
#                          susie_ebnm_state returned by susie_ebnm_fn
#                          (SuSiE alpha / mu / nu for every block, plus the
#                          shared mixture prior)
#   gfa_res$fit$L_pm       posterior mean of the WHITENED loadings X %*% ell
#                          (one row per whitened row, not per variant)
#
# What the PIPs mean
#   PIP[i, k] = P(variant i has a nonzero effect on factor k), i.e.
#   1 - prod_m (1 - alpha_mi * P(beta_m != 0 | gamma_m = i)).  It is
#   conditional (mean-field) on the other factors' posterior means.
#
# Posterior mean of L
#   ell_pm is E[ell] at the variant level, taken straight from the SuSiE state.
#   This is what t(W_b) %*% L_pm_b (i.e. R^{-1/2} times the whitened posterior
#   mean) gives when the whitening kept every direction (keep_frac small enough
#   that nothing was truncated).  If directions were truncated, t(W_b) %*% L_pm_b
#   is only the projection of E[ell] onto the retained directions, so ell_pm is
#   the better estimate.  Pass W = w$W to get both and compare.
#
# Arguments
#   gfa_res   output of gfa_fit (or a flashier fit object with $L_ghat)
#   var_index integer vector of the ORIGINAL variant indices in the order the
#             blocks were passed to susie_ebnm_fn, i.e.
#               var_index <- unlist(lapply(unique(blocks), function(b) ix[blocks == b]))
#   W         optional: w$W from make_whitening (adds L_via_Rinvhalf)
#   variant_names optional names for the rows
#
# Returns a list
#   pip          variants x K matrix of PIPs
#   ell_pm       variants x K posterior mean of L on flashier's internal scale
#   L_hat        ell_pm rescaled the way gfa_wrapup rescales L_hat (see below)
#   row_scale    the per-factor scale used for L_hat
#   L_via_Rinvhalf  (only if W given) t(W_b) %*% L_pm_b, same scale as ell_pm
#   var_index    the variant indices for the rows

gfa_susie_loadings <- function(gfa_res, W = NULL, variant_names = NULL) {
  fl <- if (!is.null(gfa_res$fit) && !is.null(gfa_res$fit$L_ghat)) gfa_res$fit else gfa_res
  if (is.null(fl$L_ghat))
    stop("No $L_ghat found. If gfa_fit hit max_iter the fit is still in gfa_res$fit; ",
         "if it was wrapped, pass gfa_res$fit directly.")
  K <- length(fl$L_ghat)
  states <- fl$L_ghat
  ok <- vapply(states, inherits, TRUE, what = "susie_ebnm_state")
  if (!all(ok))
    stop("L_ghat[[", which(!ok)[1], "]] is not a susie_ebnm_state; was susie_ebnm_fn ",
         "used for the loadings of every factor?")

  J <- nrow(fl$L_pm)
  pip <- ell <- matrix(NA_real_, J, K)
  for (k in seq_len(K)) {
    sm <- susie_state_summary(states[[k]])
    pk <- unlist(lapply(sm, `[[`, "pip"), use.names = FALSE)
    ek <- unlist(lapply(sm, `[[`, "ell_pm"), use.names = FALSE)
    pip[, k] <- pk; ell[, k] <- ek
  }
  rn <- if (!is.null(variant_names)) variant_names else 1:J
  dimnames(pip) <- dimnames(ell) <- list(as.character(rn), paste0("factor", seq_len(K)))

  # Scale used by gfa_wrapup: F is divided by `scale`, normalised to unit length,
  # and L is multiplied by the corresponding norm.
  F_hat <- fl$F_pm
  if (!is.null(gfa_res$scale)) F_hat <- F_hat / gfa_res$scale
  row_scale <- sqrt(colSums(F_hat^2))
  L_hat <- t(t(ell) * row_scale)
  dimnames(L_hat) <- dimnames(ell)


  out <- list(pip = pip, ell_pm = ell, L_hat = L_hat, row_scale = row_scale)

  if (!is.null(W)) {
    L_est <- t(t(fl$L_pm)*row_scale)
    rows <- vapply(W, nrow, 1L); cols <- vapply(W, ncol, 1L)
    if (sum(rows) != nrow(fl$L_pm)) stop("W does not match the rows of L_pm.")
    r_end <- cumsum(rows)
    r_start <- r_end - rows + 1L
    Lw <- lapply(seq_along(W), function(b){
      t(W[[b]]) %*% L_est[r_start[b]:r_end[b], , drop = FALSE]   # n_b x K
    }) %>% Reduce(rbind, .)
    dimnames(Lw) <- dimnames(ell)
    out$L_via_Rinvhalf <- Lw
  }
  out
}
