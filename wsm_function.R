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
  
  # ---- Min-Max Normalize data ---- 
  
  # Beneficial and Non beneficial functions
  min_max_ben <- function(x) {
    norm = (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    return(norm)
  }
  
  min_max_NONben <- function(x) {
    norm = (max(x, na.rm = TRUE) - x) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    return(norm)
  }
  
  
  # Apply the benfical and nonbenfical functions to the correct criteria  
  for (i in names(criteria_type)) {
    if (criteria_type[[i]] == "beneficial") {
      dataset_filter[, i] <- min_max_ben(dataset_filter[, i])
    } else {
      dataset_filter[, i] <- min_max_NONben(dataset_filter[, i])
    }
  }
  
  #---- Calculate weighted sum scores --- 
  #wsm_scores <- rowSums(criteria_weights * dataset_filter) # dataset_filter might NOT align with criteria_weights
  wsm_scores <- rowSums(mapply(function(col, w) # iterating over each column (c) of dataset_filter & its corresponding weight
    col * w, dataset_filter, # multiply column & its matching weight 
    criteria_weights[names(dataset_filter)])) # reorders the weights vectors to match dataset_filter order, so everything get multi correctly
  
  # Add scores to dataset 
  dataset$wsm_scores <- wsm_scores
  
  
  # Add letters 
  dataset <- dataset %>% 
    mutate(score_letter = case_when(
      between(wsm_scores, 0, 0.25) ~ 'D',
      between(wsm_scores, 0.25, 0.50) ~ 'C',
      between(wsm_scores, 0.50, 0.75) ~ 'B',
      wsm_scores > 0.75 ~ 'A',
      TRUE ~ NA))
  
  # dataset <- dataset %>% 
  #   mutate(
  #     quartile = ntile(wsm_scores, n = 5), # Higher groups better 
  #     quartile_label = case_when( # Add letters based on quartiles 
  #       quartile == 1 ~ "F", 
  #       quartile == 2 ~ "D",
  #       quartile == 3 ~ "C",
  #       quartile == 4 ~ "B", 
  #       quartile == 5 ~ "A"
  #     )
  #   ) %>% 
  #   select(!quartile) # Remove this column 
  
  # return the dataset 
  
  return(dataset)
  
}
