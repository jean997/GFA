#'@title Plot Factors from a matrix
#'@export
plot_factors <- function(x, names, factor_names, trait_order){
  if(is.null(x)) return(NULL)
  if(missing(names)) names <- trait_order <-  seq(nrow(x))
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

#'@title Plot factors as barplots from a matrix
#'@export
plot_factors_bars <- function(x, trait_names, factor_names, trait_order, which_factors = "all"){
  if(is.null(x)) return(NULL)
  if(missing(names)) names <- trait_order <-  seq(nrow(x))
  if(missing(factor_names)) factor_names <- paset0("Factor ", seq(ncol(x)))
  if(which_factors == "all") which_factors = seq(ncol(x))
  stopifnot(is.numeric(which_factors))
  X <- data.frame(x[,which_factors])
  names(X) <- factor_names[which_factors]
  X$name <- names
  trait_order <- trait_order[trait_order %in% X$name]
  X <- mutate(X, name = factor(name, levels = trait_order))
  plt <- X %>%
    gather("factor", "value", -name) %>%
    ggplot(.) +
    geom_bar(aes(x = name, y  = value), stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~ factor, nrow = length(which_factors), scales = "free_y")
  return(plt)
}

