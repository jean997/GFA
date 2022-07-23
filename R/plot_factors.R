#'@title Plot Factors from a matrix
#'@param row_names Labels for rows
#'@params col_names Labels for cols
#'@row_order
#'@col_order
#'@export
plot_factors <- function(x, row_names, col_names, row_order, col_order,
                         row_title = "Trait", col_title = "Factor"){
  if(is.null(x)) return(NULL)
  if(missing(row_names)) row_names <- row_order <-  seq(nrow(x))
  if(missing(col_names)) col_names <- col_order <- seq(ncol(x))
  if(missing(row_order)) row_order <- seq(nrow(x))
  if(missing(col_order)) col_order <- seq(ncol(x))
  meltx <- melt(x)
  names(meltx)[1:2] <- c("R", "C")
  meltx <- meltx %>%
    mutate( R = row_names[R] %>%
                factor(., levels = row_names[row_order]),
            C = col_names[C] %>%
                factor(., levels = col_names[col_order]))
  ggplot(data = meltx, aes(x=C, y=R, fill=value)) +
    ylab(row_title) +
    xlab(col_title) +
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

