
#'@title Plot Factors for flashier object
#'@export
plot_factors_flashier <- function(flashier_fit, names){
  x <-flashier_fit$loadings.pm[[2]]
  meltx <- melt(x) %>%
    rename(Trait = Var1, Factor = Var2) %>%
    mutate( Trait = names[Trait],
            Factor = as.factor(Factor))
  ggplot(data = meltx, aes(x=Factor, y=Trait, fill=value)) +
    geom_tile() +
    scale_fill_gradient2()
}
