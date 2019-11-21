
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


#'@title Plot Factors from a matrix
#'@export
plot_factors <- function(x, names, factor_names, trait_order){
  if(missing(factor_names)) factor_names <- seq(ncol(x))
  meltx <- melt(x) %>%
    rename(Trait = Var1, Factor = Var2) %>%
    mutate( Trait = names[Trait],
            Factor = as.factor(factor_names[Factor]))
  meltx$Trait = factor(meltx$Trait, levels = names[trait_order])
  ggplot(data = meltx, aes(x=Factor, y=Trait, fill=value)) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle=90))
}
