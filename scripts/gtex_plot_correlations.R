#' Purpose: Plot target-surrogate correlations by chromosome.
#' Updated: 2021-11-25
library(cowplot)
library(dplyr)
library(data.table)
library(ggplot2)

# -----------------------------------------------------------------------------
# Plot correlations.
# -----------------------------------------------------------------------------

rho_blood <- readRDS(file = "correlations/ssn_blood_rho_master.rds") %>%
  dplyr::mutate(
    chr_factor = factor(chr),
    r2 = pearson^2,
    tissue = "Blood"
  ) %>%
  dplyr::select(tissue, chr_factor, gene, r2)

rho_cere <- readRDS(file = "correlations/ssn_cere_rho_master.rds") %>%
  dplyr::mutate(
    chr_factor = factor(chr),
    r2 = pearson^2,
    tissue = "Cerebellum"
  ) %>%
  dplyr::select(tissue, chr_factor, gene, r2)

rho_muscle <- readRDS(file = "correlations/ssn_muscle_rho_master.rds") %>%
  dplyr::mutate(
    chr_factor = factor(chr),
    r2 = pearson^2,
    tissue = "Muscle"
  ) %>%
  dplyr::select(tissue, chr_factor, gene, r2)

rho <- rbind(rho_blood, rho_cere, rho_muscle)

rho %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(
    n_gene = dplyr::n()
  )

q <- ggplot2::ggplot(data = rho) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top",
  ) +
  ggplot2::geom_boxplot(
    aes(x = chr_factor, y = r2, fill = tissue)
  ) +
  ggplot2::scale_x_discrete(
    name = "Chromosome"
  ) + 
  ggplot2::scale_y_continuous(
    name = expression(r^2)
  ) + 
  ggplot2::scale_fill_manual(
    name = "Tissue",
    values = c("#DB4437", "#4285F4", "#F4B400")
  )

ggsave(
  filename = "figures/gtex_correlations.pdf",
  device = "pdf",
  plot = q,
  width = 7.5,
  height = 3.0,
  units = "in",
  dpi = 360
)