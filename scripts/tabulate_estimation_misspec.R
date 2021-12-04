#' Purpose: Tabulate model misspecification simulations.
#' Updated: 2021-11-09

library(dplyr)
library(xtable)

FormatSE <- function(x, y) {
  out <- sprintf("%.3f (%.3f)", sqrt(x), sqrt(y))
  return(out)
}

# -----------------------------------------------------------------------------

# Normal.
data_norm <- data.table::fread(file = "results/misspec/norm/EstH0.txt") %>%
  dplyr::filter(Parameter == "b.biv") %>%
  dplyr::filter(ms == 0) %>%
  dplyr::mutate(
    bias = round(Bias, digits = 4),
    dist = "norm",
    se = FormatSE(ModelVar, EmpiricalVar)
  ) %>%
  dplyr::select(dist, rho, mt, bias, se)


# Exponential.
data_exp <- data.table::fread(file = "results/misspec/exp/EstH0.txt") %>%
  dplyr::filter(Parameter == "b.biv") %>%
  dplyr::filter(ms == 0) %>%
  dplyr::mutate(
    bias = round(Bias, digits = 4),
    dist = "exp",
    se = FormatSE(ModelVar, EmpiricalVar)
  ) %>%
  dplyr::select(dist, rho, mt, bias, se)

# Log-normal.
data_lognorm <- data.table::fread(file = "results/misspec/lognorm/EstH0.txt") %>%
  dplyr::filter(Parameter == "b.biv") %>%
  dplyr::filter(ms == 0) %>%
  dplyr::mutate(
    bias = round(Bias, digits = 4),
    dist = "lognorm",
    se = FormatSE(ModelVar, EmpiricalVar)
  ) %>%
  dplyr::select(dist, rho, mt, bias, se)

# Student.
data_student <- data.table::fread(file = "results/misspec/student/EstH0.txt") %>%
  dplyr::filter(Parameter == "b.biv") %>%
  dplyr::filter(ms == 0) %>%
  dplyr::mutate(
    bias = round(Bias, digits = 4),
    dist = "student",
    se = FormatSE(ModelVar, EmpiricalVar)
  ) %>%
  dplyr::select(dist, rho, mt, bias, se)

# Format table.
data <- rbind(data_norm, data_exp, data_lognorm, data_student) %>%
  tidyr::pivot_wider(
    id_cols = c("rho", "mt"), 
    names_from = "dist",
    values_from = c("bias", "se")
  ) %>%
  dplyr::mutate(
    rho = rho / 100,
    mt = mt / 100,
  ) %>%
  dplyr::select(
    rho, mt, 
    bias_norm, se_norm, 
    bias_exp, se_exp, 
    bias_lognorm, se_lognorm,
    bias_student, se_student
  )

print(xtable::xtable(data, digits = 3), include.rownames = FALSE)

