#' Purpose: Tabulate type I error and non-centrality for misspecification simulations.
#' Updated: 2021-11-09
library(dplyr)
library(xtable)

#' Load Rejection Probability
#' 
#' @param dist Distribution.
#' @param file File path.
#' @param is_t1e Is this rejection probability a t1e?
#' @param is_uni Return unilateral results?
#' @return Data.frame.

LoadRejection <- function(dist, file, is_t1e = TRUE, is_uni = TRUE) {
  out <- data.table::fread(file = file) %>%
    dplyr::filter(alpha == 5e-2) %>%
    dplyr::mutate(
      dist = dist, 
      power = power * 1e2,
      rho = rho / 100,
      mt = mt / 100,
      ms = ms / 100
    ) 
  
  if (is_uni) {
    out <- out %>%
      dplyr::filter(ms == 0) 
  } else {
    out <- out %>%
      dplyr::filter(ms > 0)
  }
    
  if (is_t1e) {
    out <- out %>% 
      dplyr::rename(
        t1e = power, 
        t1e_se = se_power,
        null_ncp = ncp,
        null_ncp_se = se_ncp
      ) 
    if (is_uni) {
      out <- out %>%
        dplyr::select(dist, rho, mt, t1e, t1e_se, null_ncp, null_ncp_se)
    } else {
      out <- out %>%
        dplyr::select(dist, rho, mt, ms, t1e, t1e_se, null_ncp, null_ncp_se)
    }
    
  } else {
    out <- out %>% 
      dplyr::filter(pve == 0.5) %>%
      dplyr::rename(
        power_se = se_power,
        alt_ncp = ncp,
        alt_ncp_se = se_ncp
      ) 
    if (is_uni) {
      out <- out %>%
        dplyr::select(dist, rho, mt, power, power_se, alt_ncp, alt_ncp_se)
    } else {
      out <- out %>%
        dplyr::select(dist, rho, mt, ms, power, power_se, alt_ncp, alt_ncp_se)  
    }
  }
  return(out)
}


#' Load Data
#' 
#' @param dist Distribution.
#' @param is_uni Return unilateral results? 
#' @return Data.frame.

LoadData <- function(dist, is_uni = TRUE) {
  t1e_file <- file.path("results/misspec", dist, "Size.txt")
  pwr_file <- file.path("results/misspec", dist, "Power.txt")
  
  t1e_data <- LoadRejection(dist, t1e_file, is_t1e = TRUE, is_uni)
  pwr_data <- LoadRejection(dist, pwr_file, is_t1e = FALSE, is_uni)
  
  if (is_uni) {
    data <- t1e_data %>%
      dplyr::inner_join(pwr_data, by = c("dist", "rho", "mt")) %>%
      dplyr::select(dist, rho, mt, t1e, null_ncp, power, alt_ncp)
  } else {
    data <- t1e_data %>%
      dplyr::inner_join(pwr_data, by = c("dist", "rho", "mt", "ms")) %>%
      dplyr::select(dist, rho, mt, ms, t1e, null_ncp, power, alt_ncp) 
  }
  return(data)
}

# Normal.
data_norm <- LoadData(dist = "norm") %>% dplyr::select(-dist)
data_norm_bi <- LoadData(dist = "norm", is_uni = FALSE) %>% dplyr::select(-dist)

# Exponential.
data_exp <- LoadData(dist = "exp") %>% dplyr::select(-dist)
print(xtable::xtable(data_exp), include.rownames = FALSE)

# Log-normal.
data_lognorm <- LoadData(dist = "lognorm") %>% dplyr::select(-dist)
print(xtable::xtable(data_lognorm), include.rownames = FALSE)

# Student.
data_student <- LoadData(dist = "student") %>% dplyr::select(-dist)
print(xtable::xtable(data_student), include.rownames = FALSE)

# Norm table.
print(xtable::xtable(data_norm), include.rownames = FALSE)
print(xtable::xtable(data_norm_bi), include.rownames = FALSE)

