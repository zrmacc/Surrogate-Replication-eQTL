#' Purpose: Count the number of Bonferroni significant associations for each
#' target-surrogate pair.
#' Updated: 2021-11-25
library(dplyr)

setwd("~/Documents/Lab/Projects/Surrogates/Cluster/")

#' Tabulate Hits and NCP
#' 
#' @param data Data.frame
#' @return Data.frame

TabHits <- function(data) {
  data = data[complete.cases(data), ]
  
  # Bonferroni threshold.
  threshold <- 0.05 / nrow(data)
  
  # Hits.
  n_uni <- sum(data$uni_p <= threshold)
  n_biv <- sum(data$biv_p <= threshold)
  
  # Partition by detection method.
  uni_only <- sum(data$uni_p <= threshold & data$biv_p > threshold)
  biv_only <- sum(data$uni_p > threshold & data$biv_p <= threshold)
  both <- sum(data$uni_p <= threshold & data$biv_p <= threshold)
  
  # NCP.
  df <- data %>%
    dplyr::filter(uni_p <= threshold | biv_p <= threshold) %>%
    dplyr::mutate(
      x2_uni = (uni_beta / uni_se)^2,
      x2_biv = (biv_beta / biv_se)^2,
      var_uni = uni_se^2,
      var_biv = biv_se^2,
      re = var_uni / var_biv,
    ) %>%
    dplyr::summarise(
      ncp_uni = mean(x2_uni),
      ncp_biv = mean(x2_biv),
      mean_re = mean(re)
    )
  out <- data.frame(
    n_uni = n_uni,
    n_biv = n_biv,
    uni_only = uni_only,
    biv_only = biv_only,
    both = both,
    df
  )
}


# --------------------------------------------------------------------------------------
# Significant hit counts.
# --------------------------------------------------------------------------------------

blood <- readRDS(file = "results/ssn_blood_master.rds")
blood_tab <- TabHits(blood)
show(blood_tab)

mus <- readRDS(file = "results/ssn_muscle_master.rds")
mus_tab <- TabHits(mus)
show(mus_tab)

cere <- readRDS(file = "results/ssn_cere_master.rds")
cere_tab <- TabHits(cere)
show(cere_tab)

# --------------------------------------------------------------------------------------
# Transcripts with significant eQTL in cerebellum not expressed in blood.
# --------------------------------------------------------------------------------------
threshold_cere <- 0.05 / nrow(cere)
threshold_mus <- 0.05 / nrow(mus)

sig_cere <- cere %>%
  dplyr::filter(biv_p <= threshold_cere)
sig_mus <- mus %>%
  dplyr::filter(biv_p <= threshold_mus)

gene_cere_not_blood <- setdiff(sig_cere$gene, sig_mus$gene)
sig_cere %>%
  dplyr::filter(gene %in% transcript_cere_not_blood)
