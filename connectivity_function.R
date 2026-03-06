# =============================================================================
# connectivity_function.R
# =============================================================================
# Build connectivity matrices and reach-level dataframe from a river network
# with dams already snapped (e.g. output of st_network_blend). Use this for
# a single basin or sub-basin without re-running the full connectivity_df
# pipeline.
#
# Expects: sfnetwork with nodes that have dam_id and is_current_dam (from
# blending current and future dams via st_network_blend).
#
# Required packages: sfnetworks, igraph, dplyr, sf (load before sourcing or
# calling connectivity_from_network).
# =============================================================================

#' Connectivity summary for a snapped dam network
#'
#' Given an `sfnetwork` where dams have already been snapped/blended onto a
#' directed river network, compute directional distances between current and
#' future dams and return a one-row-per-future-dam summary.
#'
#' @md
#'
#' @param net_with_dams An `sfnetwork` created from rivers (`directed = TRUE`)
#'   and blended with dam points. Must have node columns `dam_id` and
#'   `is_current_dam` (TRUE = current, FALSE = future) and an edge column
#'   `weight` (typically metres from `edge_length()`).
#' @param cascade_threshold_km Optional numeric. If set, a future dam is
#'   `"cascade"` only if the nearest-current upstream *and* downstream distances
#'   are both `<= cascade_threshold_km` (km). If NULL, any upstream + any
#'   downstream counts as cascade.
#' @param same_trunk_only Logical. If TRUE, upstream/downstream connections and
#'   distances are restricted to current dams on the same trunk as each future
#'   dam. Trunks are identified by an edge column `bb_id` in `net_with_dams`
#'   (e.g. HydroROUT backbone IDs). Default TRUE = restrict connectivity and
#'   cascades to the same trunk; set to FALSE to use the full network regardless
#'   of trunk.
#' @param return_matrices Logical. If TRUE, return a list with `reach_df` plus
#'   the upstream/downstream distance matrices; otherwise return only `reach_df`.
#' @param subbasin Optional sf polygon. If provided, only future dams inside this
#'   polygon are included in the output. Default NULL = full network (no filter).
#' @param hybas_id Optional scalar (numeric or character). HydroBASINS sub-basin
#'   ID. When provided with `basins_sf`, restricts output to future dams inside
#'   that sub-basin (uses `basins_sf %>% filter(hybas_id == hybas_id)`).
#' @param basins_sf Optional sf object with HydroBASINS polygons and `hybas_id`
#'   column. Required when using `hybas_id`; ignored otherwise.
#'
#' @return If `return_matrices = FALSE`, a data frame with one row per future dam:
#'   `dam_id`, `dam_type`, `has_current_upstream`, `has_current_downstream`,
#'   `min_distance_upstream_km`, `min_distance_downstream_km`, `cascade_status`.
#'   If `return_matrices = TRUE`, a list with `reach_df`,
#'   `connectivity_matrix_downstream`, and `connectivity_matrix_upstream`.
#'
#' @details Distances are computed with `igraph::distances(..., mode = "out")`
#'   so travel follows edge direction (downstream). `Inf` means no path; `0`
#'   means the same node. When `same_trunk_only = TRUE`, only current dams that
#'   share the same trunk ID (`bb_id`) as each future dam are considered when
#'   evaluating upstream/downstream connectivity and minimum distances.
#'
#' @examples
#' \dontrun{
#' # Full network (no sub-basin filter)
#' reach_df <- connectivity_from_network(net_with_dams)
#'
#' out <- connectivity_from_network(net_with_dams, cascade_threshold_km = 15,
#'                                  return_matrices = TRUE)
#'
#' # Restrict to one sub-basin: get polygon from HydroBASINS (Level 4)
#' one_subbasin <- amazon_subbasins_lev4 %>%
#'   filter(hybas_id == 6040345000)
#'
#' reach_one <- connectivity_from_network(net_with_dams, subbasin = one_subbasin)
#'
#' # Or pass hybas_id + basins_sf directly (no need to create one_subbasin)
#' reach_one <- connectivity_from_network(net_with_dams, hybas_id = 6040345000,
#'                                        basins_sf = amazon_subbasins_lev4)
#'
#' # Same-trunk cascades only (requires edge column `bb_id`)
#' reach_trunk <- connectivity_from_network(net_with_dams,
#'                                          same_trunk_only = TRUE)
#' }
#'
#' @export
connectivity_from_network <- function(net_with_dams,
                                      cascade_threshold_km = NULL,
                                      same_trunk_only = TRUE,
                                      return_matrices = FALSE,
                                      subbasin = NULL,
                                      hybas_id = NULL,
                                      basins_sf = NULL) {
  
  # ---------------------------------------------------------------------------
  # Get node table and identify current vs future dam node indices
  # ---------------------------------------------------------------------------
  nodes <- net_with_dams %>%
    activate("nodes") %>%
    st_as_sf() %>%           # ensure sf
    st_drop_geometry() %>%   # drop geometry & sf class
    as_tibble() %>%
    mutate(node_id = row_number())
  
  current_dam_nodes <- nodes %>%
    filter(is_current_dam %in% TRUE) %>%
    pull(node_id)
  
  future_dam_nodes <- nodes %>%
    filter(is_current_dam %in% FALSE) %>%
    pull(node_id)
  
  # ---------------------------------------------------------------------------
  # Optional: restrict to future dams inside a sub-basin BEFORE computing
  # distance matrices (smaller matrices, faster for large basins)
  # ---------------------------------------------------------------------------
  subbasin_poly <- NULL
  if (!is.null(subbasin)) {
    subbasin_poly <- subbasin
  } else if (!is.null(hybas_id)) {
    if (is.null(basins_sf)) {
      stop("When using hybas_id, you must also provide basins_sf.", call. = FALSE)
    }
    subbasin_poly <- basins_sf %>% dplyr::filter(hybas_id == !!hybas_id)
    if (nrow(subbasin_poly) == 0L) {
      stop("No polygon found in basins_sf for hybas_id = ", hybas_id, ".", call. = FALSE)
    }
  }
  
  if (!is.null(subbasin_poly)) {
    # Future dam node geometries (order aligned with future_dam_nodes)
    future_dam_nodes_sf <- net_with_dams %>%
      activate("nodes") %>%
      dplyr::filter(is_current_dam %in% FALSE) %>%
      st_as_sf()
    
    # Use geometry (sfc) not full sf so st_filter works in all sf versions
    dam_ids_in_subbasin <- future_dam_nodes_sf %>%
      st_filter(st_geometry(subbasin_poly), .predicate = st_intersects) %>%
      dplyr::pull(dam_id)
    
    keep_idx <- which(future_dam_nodes_sf$dam_id %in% dam_ids_in_subbasin)
    future_dam_nodes <- future_dam_nodes[keep_idx]
    
    # If no future dams fall in the sub-basin, return empty result early
    if (length(future_dam_nodes) == 0L) {
      reach_df <- data.frame(
        dam_id = character(),
        dam_type = "future",
        has_current_upstream = logical(),
        has_current_downstream = logical(),
        min_distance_upstream_km = numeric(),
        min_distance_downstream_km = numeric(),
        cascade_status = character(),
        stringsAsFactors = FALSE
      )
      if (return_matrices) {
        return(list(
          reach_df = reach_df,
          connectivity_matrix_downstream = NULL,
          connectivity_matrix_upstream = NULL
        ))
      }
      return(reach_df)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Edge weights are usually produced by edge_length() / st_length().
  # We convert to numeric kilometres for igraph::distances().
  # ---------------------------------------------------------------------------
  weights <- net_with_dams %>%
    activate("edges") %>%
    pull(weight)
  
  weights_km <- tryCatch(
    {
      if (inherits(weights, "units")) {
        as.numeric(units::set_units(weights, "km", mode = "standard"))
      } else {
        as.numeric(weights) / 1000
      }
    },
    error = function(e) {
      stop(
        "Edge `weight` could not be converted to kilometres. ",
        "Make sure `weight` represents length (e.g. metres from edge_length()/st_length()) ",
        "before calling connectivity_from_network(). If your rivers are lon/lat, ",
        "either enable s2 geodesic lengths (sf::sf_use_s2(TRUE)) or project to a ",
        "metric CRS (st_transform) before computing `weight`.",
        call. = FALSE
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # Handle empty current or future dam set
  # ---------------------------------------------------------------------------
  if (length(current_dam_nodes) == 0L || length(future_dam_nodes) == 0L) {
    reach_df <- data.frame(
      dam_id = character(),
      dam_type = "future",
      has_current_upstream = logical(),
      has_current_downstream = logical(),
      min_distance_upstream_km = numeric(),
      min_distance_downstream_km = numeric(),
      cascade_status = character(),
      stringsAsFactors = FALSE
    )
    if (return_matrices) {
      return(list(
        reach_df = reach_df,
        connectivity_matrix_downstream = NULL,
        connectivity_matrix_upstream = NULL
      ))
    }
    return(reach_df)
  }
  
  # ---------------------------------------------------------------------------
  # Directional connectivity matrices (same logic as connectivity_df)
  # ---------------------------------------------------------------------------
  # Downstream: current -> future (columns = future dams)
  connectivity_matrix_downstream <- distances(
    net_with_dams,
    v = current_dam_nodes,
    to = future_dam_nodes,
    weights = weights_km,
    mode = "out"
  )
  
  # Upstream: future -> current (rows = future dams)
  connectivity_matrix_upstream <- distances(
    net_with_dams,
    v = future_dam_nodes,
    to = current_dam_nodes,
    weights = weights_km,
    mode = "out"
  )
  
  # ---------------------------------------------------------------------------
  # Future dam IDs for reach_df
  # ---------------------------------------------------------------------------
  future_dam_ids <- nodes %>%
    slice(future_dam_nodes) %>%
    pull(dam_id)
  
  # ---------------------------------------------------------------------------
  # Build reach_df: one row per future dam
  # ---------------------------------------------------------------------------
  reach_df <- data.frame(
    dam_id = future_dam_ids,
    dam_type = "future",
    stringsAsFactors = FALSE
  )
  
  # Has any current dam downstream of this future dam? (column = future dam)
  reach_df$has_current_downstream <- apply(
    connectivity_matrix_downstream, 2L,
    function(x) any(is.finite(x) & x > 0)
  )
  
  # Has any current dam upstream of this future dam? (row = future dam)
  reach_df$has_current_upstream <- apply(
    connectivity_matrix_upstream, 1L,
    function(x) any(is.finite(x) & x > 0)
  )
  
  # Min distance to nearest current dam downstream (column = future dam)
  reach_df$min_distance_downstream_km <- apply(
    connectivity_matrix_downstream, 2L,
    function(x) {
      z <- x[is.finite(x) & x > 0]
      if (length(z) > 0L) min(z) else NA_real_
    }
  )
  
  # Min distance to nearest current dam upstream (row = future dam)
  reach_df$min_distance_upstream_km <- apply(
    connectivity_matrix_upstream, 1L,
    function(x) {
      z <- x[is.finite(x) & x > 0]
      if (length(z) > 0L) min(z) else NA_real_
    }
  )
  
  # ---------------------------------------------------------------------------
  # Optional: restrict connectivity to same trunk (same bb_id)
  # ---------------------------------------------------------------------------
  if (isTRUE(same_trunk_only)) {
    edges_tbl <- net_with_dams %>%
      activate("edges") %>%
      st_as_sf() %>%
      st_drop_geometry() %>%
      as_tibble()
    
    if (!"bb_id" %in% names(edges_tbl)) {
      stop(
        "same_trunk_only = TRUE requires an edge attribute `bb_id` in net_with_dams.\n",
        "Ensure the rivers used to build the network contained `bb_id` and that\n",
        "this attribute was preserved in the edge table.",
        call. = FALSE
      )
    }
    
    # Derive a trunk id per node from incident edges
    incident_from <- edges_tbl %>%
      select(node_id = from, bb_id)
    incident_to <- edges_tbl %>%
      select(node_id = to, bb_id)
    
    incident <- bind_rows(incident_from, incident_to) %>%
      filter(!is.na(bb_id)) %>%
      group_by(node_id) %>%
      summarise(trunk_id = dplyr::first(bb_id), .groups = "drop")
    
    nodes_trunk <- nodes %>%
      left_join(incident, by = "node_id")
    
    current_trunk <- nodes_trunk %>%
      slice(current_dam_nodes) %>%
      pull(trunk_id)
    
    future_trunk <- nodes_trunk %>%
      slice(future_dam_nodes) %>%
      pull(trunk_id)
    
    reach_df$future_trunk <- future_trunk
    
    n_future <- length(future_dam_nodes)
    
    has_current_downstream_same <- logical(n_future)
    has_current_upstream_same <- logical(n_future)
    min_distance_downstream_same <- rep(NA_real_, n_future)
    min_distance_upstream_same <- rep(NA_real_, n_future)
    
    for (j in seq_len(n_future)) {
      if (is.na(future_trunk[j])) next
      
      same_idx <- which(current_trunk == future_trunk[j] & !is.na(current_trunk))
      if (length(same_idx) == 0L) next
      
      col_down <- connectivity_matrix_downstream[same_idx, j, drop = FALSE]
      row_up <- connectivity_matrix_upstream[j, same_idx, drop = FALSE]
      
      d_down <- col_down[is.finite(col_down) & col_down > 0]
      d_up <- row_up[is.finite(row_up) & row_up > 0]
      
      if (length(d_down) > 0L) {
        has_current_downstream_same[j] <- TRUE
        min_distance_downstream_same[j] <- min(d_down)
      }
      if (length(d_up) > 0L) {
        has_current_upstream_same[j] <- TRUE
        min_distance_upstream_same[j] <- min(d_up)
      }
    }
    
    # Override reach_df connectivity with same-trunk-only values
    reach_df$has_current_downstream <- has_current_downstream_same
    reach_df$has_current_upstream <- has_current_upstream_same
    reach_df$min_distance_downstream_km <- min_distance_downstream_same
    reach_df$min_distance_upstream_km <- min_distance_upstream_same
  }
  
  # ---------------------------------------------------------------------------
  # Cascade status (optional distance threshold)
  # ---------------------------------------------------------------------------
  threshold_km <- if (is.null(cascade_threshold_km)) NA_real_ else cascade_threshold_km
  
  if (is.na(threshold_km)) {
    reach_df$cascade_status <- ifelse(
      reach_df$has_current_upstream & reach_df$has_current_downstream,
      "cascade",
      "not in cascade"
    )
  } else {
    reach_df$cascade_status <- ifelse(
      reach_df$has_current_upstream & reach_df$has_current_downstream &
        !is.na(reach_df$min_distance_upstream_km) &
        !is.na(reach_df$min_distance_downstream_km) &
        reach_df$min_distance_upstream_km <= threshold_km &
        reach_df$min_distance_downstream_km <= threshold_km,
      "cascade",
      "not in cascade"
    )
  }
  
  # Sort by cascade status, then distances (match connectivity_df)
  reach_df <- reach_df %>%
    arrange(cascade_status, min_distance_upstream_km, min_distance_downstream_km)
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  if (return_matrices) {
    return(list(
      reach_df = reach_df,
      connectivity_matrix_downstream = connectivity_matrix_downstream,
      connectivity_matrix_upstream = connectivity_matrix_upstream
    ))
  }
  reach_df
}
