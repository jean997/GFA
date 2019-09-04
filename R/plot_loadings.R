#'@title Plot Loadings for flash object
#'@export
plot_loadings <- function(flash_fit){
  vals <- flash_fit$ldf$l
  dat <- reshape2::melt(vals) %>%
         rename(SNP = Var1, Factor = Var2) %>%
         mutate(Factor = as.factor(Factor))
  cols <- gtex.colors[1:ncol(vals)]
  min_val <- min(0, min(vals))
  max_val <- max(0, max(vals))

  ggplot(dat, aes(x=SNP, y=value)) +
    geom_point(size=0.5) +
    scale_x_discrete(labels = NULL) +
    theme_grey() +
    theme(legend.position="right",
          legend.text = element_text(size = 3.25),
          legend.title = element_blank()) +
    labs(y = "", x = "") +
    facet_wrap(~Factor, ncol = 2) +
    theme_bw()
}

#'@title Plot Loadings for flashier object
#'@export
plot_loadings_flashier <- function(flashier_fit){
  vals <- flashier_fit$loadings.pm[[1]]
  dat <- reshape2::melt(vals) %>%
    rename(SNP = Var1, Factor = Var2) %>%
    mutate(Factor = as.factor(Factor))
  cols <- gtex.colors[1:ncol(vals)]
  min_val <- min(0, min(vals))
  max_val <- max(0, max(vals))

  ggplot(dat, aes(x=SNP, y=value)) +
    geom_point(size=0.5) +
    scale_x_discrete(labels = NULL) +
    theme_grey() +
    theme(legend.position="right",
          legend.text = element_text(size = 3.25),
          legend.title = element_blank()) +
    labs(y = "", x = "") +
    facet_wrap(~Factor, ncol = 2) +
    theme_bw()
}
