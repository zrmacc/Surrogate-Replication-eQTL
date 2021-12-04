#' Purpose: Tabulate type I error and non-centrality for large-sample simulations.
#' Updated: 2021-11-10
library(dplyr)
library(xtable)

# -----------------------------------------------------------------------------

# 5K.
data_5k <- data.table::fread(file = "results/large/5k/Size.txt") %>%
  dplyr::filter(alpha == 0.05) %>%
  dplyr::mutate(power_5k = power * 1e2, ncp_5k = ncp) %>%
  dplyr::select(rho, mt, ms, power_5k, ncp_5k)

# 10k.
data_10k <- data.table::fread(file = "results/large/10k/Size.txt") %>%
  dplyr::filter(alpha == 0.05) %>%
  dplyr::mutate(power_10k = power * 1e2, ncp_10k = ncp) %>%
  dplyr::select(rho, mt, ms, power_10k, ncp_10k)

# 20k.
data_20k <- data.table::fread(file = "results/large/20k/Size.txt") %>%
  dplyr::filter(alpha == 0.05) %>%
  dplyr::mutate(power_20k = power * 1e2, ncp_20k = ncp) %>%
  dplyr::select(rho, mt, ms, power_20k, ncp_20k)


# Format table.
out <- merge(x = data_5k, y = data_10k, by = c("rho", "mt", "ms"))
out <- merge(x = out, y = data_20k, by = c("rho", "mt", "ms"))
out$rho <- out$rho / 100
out$mt <- out$mt / 100
out$ms <- out$ms / 100

out <- out[order(out$ms, out$mt, out$rho), ]

print(xtable(out, digits = 3), include.rownames = FALSE)

# -----------------------------------------------------------------------------

# 5K.
data_5k <- data.table::fread(file = "results/large/5k/Size.txt") %>%
  dplyr::filter(alpha == 1e-05) %>%
  dplyr::mutate(power_5k = power * 1e5, ncp_5k = ncp) %>%
  dplyr::select(rho, mt, ms, power_5k, ncp_5k)

# 10k.
data_10k <- data.table::fread(file = "results/large/10k/Size.txt") %>%
  dplyr::filter(alpha == 1e-05) %>%
  dplyr::mutate(power_10k = power * 1e5, ncp_10k = ncp) %>%
  dplyr::select(rho, mt, ms, power_10k, ncp_10k)

# 20k.
data_20k <- data.table::fread(file = "results/large/20k/Size.txt") %>%
  dplyr::filter(alpha == 1e-05) %>%
  dplyr::mutate(power_20k = power * 1e5, ncp_20k = ncp) %>%
  dplyr::select(rho, mt, ms, power_20k, ncp_20k)


# Format table.
out <- merge(x = data_5k, y = data_10k, by = c("rho", "mt", "ms"))
out <- merge(x = out, y = data_20k, by = c("rho", "mt", "ms"))
out$rho <- out$rho / 100
out$mt <- out$mt / 100
out$ms <- out$ms / 100

out <- out[order(out$ms, out$mt, out$rho), ]

print(xtable(out, digits = 3), include.rownames = FALSE)
