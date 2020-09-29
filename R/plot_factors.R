
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

if(FALSE){
  mats <- readRDS("analysis_data/bcai_gwas_mats_order1.RDS")
  traits <- str_split(mats$traits, "/") %>% map(., 2) %>%
    unlist(.) %>%
    str_replace(., ".top_summary_statistics.tsv.gz", "") %>%
    str_replace(., "bcai_astle_", "") %>%
    str_replace(., "bcai_", "")
  fit <- readRDS("analysis_data/bcai_order1__soft_impute__noisy_bycol.RDS")
x <- fit$loadings.pm[[2]]
names <- traits
factor_names <- seq(ncol(x))
new_names <- c(rep("", 36), "SLE", "RA", "CD", "IBD", "UC", "Asth", "PSC", "Allg")
trait_order <- c(1:10, 16:17 , 21:22, 25:28, 24,36, 23,29,30,31, 11:15, 18:20,32:35, 37, 38, 39, 40, 41, 43, 42, 44)
meltx <- melt(x) %>%
  rename(Trait = Var1, Factor = Var2) %>%
  mutate( Trait = names[Trait],
          Factor = as.factor(factor_names[Factor]))
meltx$Trait = factor(meltx$Trait, levels = names[trait_order])

p <-  ggplot(data = meltx, aes(x=Factor, y=Trait, fill=value)) +
      geom_tile() +
      scale_fill_gradient2(high="purple", low="darkorange", name="Factor effect\non trait") +
      scale_y_discrete(labels = new_names[trait_order]) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.ticks = element_blank(),
            axis.title = element_text(size=16),
            legend.title = element_text(size=12))
ggsave(p, file="~/Dropbox/Apps/Overleaf/job_applications/biostatistics/ucdavis/job_talk_ucd/fig/blood_cell_factors1.png",
       height=8, width=7, units="in", dpi=300)

p <-  ggplot(data = meltx, aes(x=Factor, y=Trait, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(high="purple", low="darkorange", name="") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=20),
        legend.position = "none",
        panel.border = element_rect(fill=NA, size=1.5))
ggsave(p, file="~/Dropbox/Apps/Overleaf/job_applications/biostatistics/ucdavis/job_talk_ucd/fig/blood_cell_factors0.png",
       height=8, width=7, units="in", dpi=300)

}
