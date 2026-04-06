# =============================================================================
# wsm_function.R
# =============================================================================
# Purpose: Find the Weighted Sum Model for a dataset with given criteria and weights.
#
# Expects: 
#
# Outside this function: 
# 
# Required packages:
#
# Inline comments: plain-language walkthrough for explaining to a team what each
# step does and why.
# =============================================================================

#' Connectivity summary for a snapped dam network
#'
#' Given an `sfnetwork` where dams have already been snapped/blended onto a
#' directed river network, compute directional distances between current and
#' future dams for **future** dams, append **current** dams with `NA` in
#' cascade-only columns, and return one row per dam present on the network.
#'
#' @md
#'
#' @param dataset 
#' @param criteria 
#' @param criteria_weights
#' @param criteria_type
#'
#' @return 
#'
#' @details 
#'
#' @examples


wsm_function <- function(dataset, criteria, criteria_weights, criteria_type) {
  
  # Validate weights 
  if (abs(sum(criteria_weights) != 1)) {stop("Criteria weights must sum to 1.")}
  
  # --- Filter for columns inputed as criteria --- 
  dataset_filter <- dataset %>% select(any_of(criteria))
  
  # ---- Linear Normalize data ---- 
  
  # Benifical and Non benifical functions 
  benfical_norm_fun <- function(column) {
    norm = column / max(column, na.rm = TRUE)
    return(norm)
  }
  
  nonbenfical_norm_fun <- function(column) {
    norm = min(column, na.rm = TRUE) / column
    return(norm)
  }
  
  # Apply the benfical and nonbenfical functions to the correct criteria  
  for (i in names(criteria_type)) {
    if (criteria_type[i] == "beneficial") {
      dataset_filter[, i] <- benfical_norm_fun(dataset_filter[, i])
    } else {
      dataset_filter[, i] <- nonbenfical_norm_fun(dataset_filter[, i])
    }
  }
  
  #---- Calculate weighted sum scores --- 
  wsm_scores <- rowSums(critera_weights*dataset_filter)
  
  # Add scores to dataset 
  dataset$wsm_scores <- wsm_scores
  
  
  # Add letters 
  dataset <- dataset %>% 
    mutate(
      quartile = ntile(wsm_scores, n = 5), # Higher groups better 
      quartile_label = case_when( # Add letters based on quartiles 
        quartile == 1 ~ "F", 
        quartile == 2 ~ "D",
        quartile == 3 ~ "C",
        quartile == 4 ~ "B", 
        quartile == 5 ~ "A"
      )
    ) %>% 
    select(!quartile) # Remove this column 
  
  # return the dataset 
  return(dataset)
  
}
