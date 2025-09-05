#'@title Plot Factors from a matrix
#'@param row_names Labels for rows
#'@param col_names Labels for cols
#'@param row_order Order of rows
#'@param col_order Order of columns
#'@param row_title Title for horizontal axis
#'@param col_title Title for vertical axis
#'@return A ggplot object
#'@export
plot_factors <- function(x, row_names, col_names, row_order, col_order,
                         row_title = "Trait", col_title = "Factor"){
  if(is.null(x)) return(NULL)
  if(missing(row_names)) row_names <-  seq(nrow(x))
  if(missing(col_names)) col_names <-  seq(ncol(x))
  if(missing(row_order)) row_order <- seq(nrow(x))
  if(missing(col_order)) col_order <- seq(ncol(x))
  rownames(x) <- row_names
  colnames(x) <- col_names
  x <- x[row_order, col_order]
  meltx <- melt(x)
  names(meltx)[1:2] <- c("R", "C")
  meltx <- meltx %>%
     mutate( R = factor(R, levels = row_names[row_order]),
             C = factor(C, levels = col_names[col_order]))
  ggplot(data = meltx, aes(x=C, y=R, fill=value)) +
    ylab(row_title) +
    xlab(col_title) +
    geom_tile() +
    scale_fill_gradient2() +
    theme(axis.text.x = element_text(angle=90))
}

#'@title Plot factors as barplots from a matrix
#'@export
plot_factors_bars <- function(x, trait_names, factor_names, trait_order, which_factors = seq(ncol(x))){
  if(is.null(x)) return(NULL)
  if(missing(trait_names)) trait_names <- seq(nrow(x))
  if(missing(trait_order)) trait_order <- trait_names
  if(missing(factor_names)) factor_names <- paste0("Factor ", seq(ncol(x)))
  stopifnot(is.numeric(which_factors))
  X <- data.frame(x)
  names(X) <- factor_names
  X$name <- trait_names
  trait_order <- trait_order[trait_order %in% X$name]
  X <- mutate(X, name = factor(name, levels = trait_order))
  plt <- X %>%
    pivot_longer(names_to = "factor", values_to = "value", cols =  all_of(factor_names[which_factors])) %>%
    ggplot(.) +
    geom_bar(aes(x = name, y  = value), stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~ factor, nrow = length(which_factors), scales = "free_y")
  return(plt)
}

