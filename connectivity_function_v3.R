# =============================================================================
# connectivity_function_v3.R
# =============================================================================
# Connectivity for a blended dam network (v3): builds `reach_df` (one row per
# dam node) plus a `debug` list of tibbles.
#
# Do not source connectivity_function.R in the same session; this file defines
# connectivity_from_network() with v3 behavior.
#
# Prerequisites: directed `sfnetwork` with edge `weight` (metres) and `bb_id`;
# blended nodes carry `dam_id` and `is_current_dam` (TRUE = current, FALSE =
# future). Load: sfnetworks, igraph, dplyr, sf, tidygraph before sourcing.
#
# Contract summary (v3):
# - Downstream “next” current dam: among future→current paths with finite
#   dist_km > 0 and within threshold_downstream_km, **fewest trunk hops** (bb_id
#   confluence graph), then **smallest river km**. `downstream_hops` = hop count.
# - Upstream “next” current dam: same tier rule on current→future paths within
#   threshold_upstream_km. `upstream_hops` = hop count.
# - cascade_level = max(downstream_hops, upstream_hops) when **both**
#   neighbors exist and both hop counts are finite; otherwise NA (“no cascade”).
# - connectivity_category: cascade_classic when both sides are same-trunk
#   (downstream_hops == 0 and upstream_hops == 0), then cascade1/2/3+ from
#   cascade_level, then FFR / downstream
#   / upstream arms from min distances vs the same threshold_*_km used for pools.
# - Thresholds are explicit km inputs: downstream/upstream/cascade.
# =============================================================================

#' Connectivity for a snapped dam network (`reach_df` + `debug`) — v3
#'
#' Builds directed river distances between **current** and **future** dam nodes
#' on a blended `sfnetwork`, applies tiered upstream/downstream rules with km
#' caps, trunk hop counts, `cascade_level`, and `connectivity_category`.
#'
#' **Dam nodes:** nodes with non-missing `dam_id`; `is_current_dam` is `TRUE` for
#' current dams and `FALSE` for future dams.
#'
#' **Edge weights:** column `weight` is treated as **metres** and converted to
#' **km** for [igraph::distances()] (`as.numeric(weight) / 1000`). Edges must
#' carry trunk id `bb_id` (confluence trunk graph).
#'
#' **Downstream neighbor:** finite future→current path, `dist_km > 0`, at most
#' `threshold_downstream_km`; pick **fewest trunk hops** then smallest `dist_km`.
#'
#' **Upstream neighbor:** finite current→future path, `dist_km > 0`, at most
#' `threshold_upstream_km`; same tier rule.
#'
#' **`cascade_level`:** `pmax(downstream_hops, upstream_hops)` when both
#' neighbors exist and both hop counts are non-NA; else `NA_integer_`.
#'
#' **`connectivity_category` cascade classes:** `cascade_classic` when both
#' neighbors are on the same trunk (`downstream_hops == 0` and
#' `upstream_hops == 0`), then `cascade1`, `cascade2`, `cascade3+` for max hops
#' 1, 2, and 3+ respectively.
#'
#' **`threshold_cascade_km`:** explicit km input (currently not used to pick a
#' second neighbor; kept as explicit cascade threshold parameter).
#'
#' @param net_with_dams Directed `sfnetwork` with edge `bb_id`, numeric edge
#'   `weight` (metres), and blended node columns `dam_id` / `is_current_dam`.
#' @param threshold_downstream_km Required km cap on downstream candidate paths.
#' @param threshold_upstream_km Required km cap on upstream candidate paths.
#' @param threshold_cascade_km Required cascade threshold in km (explicit input).
#'
#' @return Named list: `reach_df`, `debug` (intermediate tables), and
#'   `threshold_used` (named numeric km actually applied).
#'
#' @details
#' Distances use `igraph::distances(..., mode = "out")`. `Inf` = no path; `0` =
#' same vertex (excluded from candidate pools).
#'
#' @examples
#' \dontrun{
#' conn <- connectivity_from_network(
#'   out$net_with_dams,
#'   threshold_downstream_km = 200,
#'   threshold_upstream_km = 200,
#'   threshold_cascade_km = 200
#' )
#' reach_df <- conn$reach_df
#' }
#'
#' @md
#' @export
connectivity_from_network <- function(
    net_with_dams,
    threshold_downstream_km,
    threshold_upstream_km,
    threshold_cascade_km) {
  
  # Input checks + threshold capture
  #
  # Final return list always has: reach_df, debug, threshold_used.
  # cascade_level rule later in SECTION 9:
  #   - NA if either up/down dam link or hop value is missing
  #   - otherwise max(downstream_hops, upstream_hops)
  
  threshold_used <- c(
    threshold_downstream_km = threshold_downstream_km,
    threshold_upstream_km = threshold_upstream_km,
    threshold_cascade_km = threshold_cascade_km
  )
  
  # =============================================================================
  # SECTION 1 — flatten sfnetwork to node/edge tables
  #
  # INPUTS:
  #   - net_with_dams sfnetwork object
  # OUTPUTS:
  #   - nodes_tbl (node attributes, no geometry, node_id index)
  #   - edges_tbl (edge attributes, no geometry)
  #   - weights_km (numeric edge weights for igraph::distances())
  #   - node_to_bb_lookup (one representative bb_id per node)
  # DEBUG OBJECTS:
  #   - nodes_tbl, edges_tbl, weights_km, node_to_bb_lookup
  # =============================================================================
  
  nodes_tbl <- net_with_dams %>% # start from the sfnetwork object
    tidygraph::activate("nodes") %>% # switch active table to node attributes
    sf::st_as_sf() %>% # convert active node table to sf object
    sf::st_drop_geometry() %>% # remove POINT geometry; keep only attributes
    dplyr::as_tibble() %>% # standardize to tibble for dplyr pipelines
    dplyr::mutate(node_id = dplyr::row_number()) # create stable node index used in pair tables
  
  edges_tbl <- net_with_dams %>%
    tidygraph::activate("edges") %>%
    sf::st_as_sf() %>%
    sf::st_drop_geometry() %>%
    dplyr::as_tibble()
  
  weights_km <- as.numeric(edges_tbl$weight) / 1000
  
  # One representative bb_id per node (needed for early exits and pair tables).
  edge_from_node_bb <- edges_tbl %>% dplyr::select(node_id = .data$from, bb_id = .data$bb_id)
  edge_to_node_bb <- edges_tbl %>% dplyr::select(node_id = .data$to, bb_id = .data$bb_id)
  edge_node_bb_incidence <- dplyr::bind_rows(edge_from_node_bb, edge_to_node_bb)
  node_to_bb_lookup <- edge_node_bb_incidence %>%
    dplyr::filter(!is.na(.data$bb_id)) %>%
    dplyr::group_by(.data$node_id) %>%
    dplyr::summarise(bb_id = dplyr::first(.data$bb_id), .groups = "drop")
  
  # =============================================================================
  # SECTION 2 — split dam nodes into current vs future
  #
  # INPUTS:
  #   - nodes_tbl from SECTION 1
  # OUTPUTS:
  #   - dam_nodes_tbl, current_dams_tbl, future_dams_tbl
  #   - index vectors current_nodes / future_nodes + dam_id vectors
  # DEBUG OBJECTS:
  #   - dam_nodes_tbl, current_dams_tbl, future_dams_tbl
  # =============================================================================
  
  dam_nodes_tbl <- nodes_tbl %>% # begin with all nodes
    dplyr::filter(!is.na(.data$dam_id)) %>% # keep only snapped dam nodes (non-missing id)
    dplyr::mutate(dam_id = as.character(.data$dam_id)) # force character id for stable joins
  
  current_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% TRUE) %>%
    dplyr::select(node_id, dam_id, is_current_dam)
  
  future_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% FALSE) %>%
    dplyr::select(node_id, dam_id, is_current_dam)
  
  current_nodes <- current_dams_tbl$node_id
  future_nodes <- future_dams_tbl$node_id
  current_dam_ids <- current_dams_tbl$dam_id
  future_dam_ids <- future_dams_tbl$dam_id
  n_cur <- length(current_nodes)
  n_fut <- length(future_nodes)
  
  if (n_fut == 0L) {
    return(list(
      reach_df = empty_reach_connectivity_df(),
      debug = list(note = "no future dam nodes on network"),
      threshold_used = threshold_used
    ))
  }
  
  if (n_cur == 0L) {
    # No current dams: all connectivity NA; still return one row per future.
    future_trunks <- future_dams_tbl %>%
      dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
      dplyr::transmute(
        future_dam_id = as.character(.data$dam_id),
        future_node_id = .data$node_id,
        bb_id = .data$bb_id
      )
    
    reach_future <- data.frame(
      dam_id = future_dam_ids,
      dam_type = rep("future", n_fut),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(
        dplyr::transmute(future_trunks, dam_id = .data$future_dam_id, bb_id = .data$bb_id),
        by = "dam_id"
      ) %>%
      dplyr::mutate(
        has_current_downstream = FALSE,
        min_distance_downstream_km = NA_real_,
        dam_id_down = NA_character_,
        downstream_hops = NA_integer_,
        has_current_upstream = FALSE,
        min_distance_upstream_km = NA_real_,
        dam_id_up = NA_character_,
        upstream_hops = NA_integer_,
        cascade_level = NA_integer_,
        connectivity_category = NA_character_
      )
    
    reach_df <- reach_future
    return(list(
      reach_df = reach_df,
      debug = list(note = "no current dam nodes on network"),
      threshold_used = threshold_used
    ))
  }
  
  # =============================================================================
  # SECTION 3 — compute directed distance matrices (km)
  #
  # INPUTS:
  #   - net_with_dams graph, current_nodes, future_nodes, weights_km
  # OUTPUTS:
  #   - ds_mat (future -> current), us_mat (current -> future)
  # DEBUG OBJECTS:
  #   - ds_mat, us_mat
  # =============================================================================
  
  ds_mat <- igraph::distances( # directed shortest-path distances on graph
    net_with_dams, # graph object (sfnetwork/tidygraph backed igraph)
    v = future_nodes, # source vertices: future dams
    to = current_nodes, # target vertices: current dams
    weights = weights_km, # edge weights in km
    mode = "out" # follow edge direction from source downstream
  )
  
  us_mat <- igraph::distances(
    net_with_dams,
    v = current_nodes,
    to = future_nodes,
    weights = weights_km,
    mode = "out"
  )
  
  # =============================================================================
  # SECTION 4 — melt matrices into row-wise pair tables
  #
  # INPUTS:
  #   - ds_mat / us_mat and node/dam id index vectors
  # OUTPUTS:
  #   - ds_dam_all, us_dam_all (one row per future/current pair)
  # DEBUG OBJECTS:
  #   - ds_grid, us_grid, ds_dam_all, us_dam_all
  # =============================================================================
  
  ds_grid <- expand.grid(fut_idx = seq_len(n_fut), cur_idx = seq_len(n_cur)) # all (future, current) index pairs
  ds_grid$dist_km <- as.vector(ds_mat) # flatten matrix to align one distance per index-pair row
  
  ds_dam_all <- ds_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx],
      current_node_id = current_nodes[.data$cur_idx],
      future_dam_id = future_dam_ids[.data$fut_idx],
      current_dam_id = current_dam_ids[.data$cur_idx]
    ) %>%
    dplyr::select(
      future_dam_id = .data$future_dam_id,
      current_dam_id = .data$current_dam_id,
      future_node_id = .data$future_node_id,
      current_node_id = .data$current_node_id,
      dist_km = .data$dist_km
    )
  
  us_grid <- expand.grid(cur_idx = seq_len(n_cur), fut_idx = seq_len(n_fut))
  us_grid$dist_km <- as.vector(us_mat)
  
  us_dam_all <- us_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx],
      current_node_id = current_nodes[.data$cur_idx],
      future_dam_id = future_dam_ids[.data$fut_idx],
      current_dam_id = current_dam_ids[.data$cur_idx]
    ) %>%
    dplyr::select(
      future_dam_id = .data$future_dam_id,
      current_dam_id = .data$current_dam_id,
      future_node_id = .data$future_node_id,
      current_node_id = .data$current_node_id,
      dist_km = .data$dist_km
    )
  
  # =============================================================================
  # SECTION 5 — join bb_id (trunk) to every future/current pair row
  #
  # INPUTS:
  #   - ds_dam_all / us_dam_all + node_to_bb_lookup
  # OUTPUTS:
  #   - pair tables with future_bb_id/current_bb_id
  #   - future_trunks / current_trunks lookup tables
  # DEBUG OBJECTS:
  #   - ds_dam_all, us_dam_all, future_trunks, current_trunks
  # =============================================================================
  
  ds_dam_all <- ds_dam_all %>%
    dplyr::left_join(node_to_bb_lookup, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = .data$bb_id) %>%
    dplyr::left_join(node_to_bb_lookup, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = .data$bb_id)
  
  us_dam_all <- us_dam_all %>%
    dplyr::left_join(node_to_bb_lookup, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = .data$bb_id) %>%
    dplyr::left_join(node_to_bb_lookup, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = .data$bb_id)
  
  future_trunks <- future_dams_tbl %>%
    dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
    dplyr::rename(future_node_id = .data$node_id, future_dam_id = .data$dam_id, bb_id = .data$bb_id) %>%
    dplyr::select(.data$future_dam_id, .data$future_node_id, .data$bb_id)
  
  current_trunks <- current_dams_tbl %>%
    dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
    dplyr::rename(current_node_id = .data$node_id, current_dam_id = .data$dam_id, bb_id = .data$bb_id) %>%
    dplyr::select(.data$current_dam_id, .data$current_node_id, .data$bb_id)
  
  # =============================================================================
  # SECTION 6 — build undirected trunk graph and hop lookup matrix
  #
  # INPUTS:
  #   - edges_tbl endpoint incidence by node and bb_id
  # OUTPUTS:
  #   - bb_confluence_edges, bb_confluence_graph, bb_hop_distance_matrix
  #   - trunk_hops_vec() helper for vectorized hop lookup
  # DEBUG OBJECTS:
  #   - bb_confluence_edges, bb_confluence_graph, bb_hop_distance_matrix
  # =============================================================================
  
  junction_nodes <- unique(c(edges_tbl$from, edges_tbl$to)) # node ids that appear as any edge endpoint
  bb_confluence_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE) # accumulator of bb_id-to-bb_id links
  
  for (junction_node_id in junction_nodes) { # inspect each graph junction separately
    bb_ids_at_node <- edge_node_bb_incidence %>% dplyr::filter(.data$node_id == junction_node_id) %>% dplyr::pull(.data$bb_id) # collect incident trunk ids
    bb_ids_at_node <- unique(as.character(bb_ids_at_node)) # deduplicate and normalize type
    bb_ids_at_node <- bb_ids_at_node[!is.na(bb_ids_at_node) & nzchar(bb_ids_at_node)] # drop missing/blank bb ids
    if (length(bb_ids_at_node) >= 2) { # confluence only exists when 2+ trunks meet
      bb_id_pairs_at_node <- utils::combn(bb_ids_at_node, 2) # enumerate all unordered bb_id pairs
      bb_confluence_edges <- rbind( # append these pairs into the global edge list
        bb_confluence_edges,
        data.frame(
          from = as.character(bb_id_pairs_at_node[1, ]),
          to = as.character(bb_id_pairs_at_node[2, ]),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  bb_confluence_edges <- unique(bb_confluence_edges)
  bb_confluence_graph <- igraph::graph_from_data_frame(bb_confluence_edges, directed = FALSE)
  bb_hop_distance_matrix <- igraph::distances(bb_confluence_graph, mode = "all")
  if (igraph::vcount(bb_confluence_graph) > 0) {
    bb_vertex_names <- igraph::V(bb_confluence_graph)$name
    dimnames(bb_hop_distance_matrix) <- list(bb_vertex_names, bb_vertex_names)
  }
  
  # Helper: hop matrix lookup vectorized
  lookup_bb_hops_for_pairs <- function(future_bb_ids, current_bb_ids, hop_matrix) {
    future_bb_ids <- as.character(future_bb_ids)
    current_bb_ids <- as.character(current_bb_ids)
    matrix_row_bb_ids <- rownames(hop_matrix)
    if (length(matrix_row_bb_ids) == 0L) {
      return(rep(NA_integer_, length(future_bb_ids)))
    }
    future_bb_row_idx <- match(future_bb_ids, matrix_row_bb_ids)
    current_bb_col_idx <- match(current_bb_ids, colnames(hop_matrix))
    valid_matrix_indices <- !is.na(future_bb_row_idx) & !is.na(current_bb_col_idx)
    hop_counts_out <- rep(NA_integer_, length(future_bb_ids))
    hop_counts_out[valid_matrix_indices] <- as.integer(
      hop_matrix[cbind(future_bb_row_idx[valid_matrix_indices], current_bb_col_idx[valid_matrix_indices])]
    )
    hop_counts_out
  }
  
  # =============================================================================
  # SECTION 7 — downstream candidate pool + tiered winner
  #
  # IMPORTANT: this naturally does "same-trunk first" because same trunk is
  # trunk_step_dist == 0. If no same-trunk candidate is within threshold, the
  # next smallest hop tier (1, 2, ...) is used.
  #
  # INPUTS:
  #   - ds_dam_all + threshold_downstream_km + bb_hop_distance_matrix
  # OUTPUTS:
  #   - ds_reachable, ds_reachable_with_steps, ds_within_cap
  #   - chosen_downstream (one winner per future)
  #   - downstream_summary (pool diagnostics by future)
  # DEBUG OBJECTS:
  #   - ds_reachable, ds_same_trunk, ds_reachable_with_steps,
  #     ds_within_cap, chosen_downstream, downstream_summary
  # =============================================================================
  
  ds_reachable <- ds_dam_all %>% # start from all downstream direction pairs
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0) # keep only real reachable paths (exclude Inf/self 0)
  
  ds_same_trunk <- ds_reachable %>%
    dplyr::filter(.data$current_bb_id == .data$future_bb_id)
  
  future_bb_ids_ds <- as.character(ds_reachable$future_bb_id)
  current_bb_ids_ds <- as.character(ds_reachable$current_bb_id)
  ds_reachable_with_steps <- ds_reachable %>%
    dplyr::mutate(trunk_step_dist = lookup_bb_hops_for_pairs(future_bb_ids_ds, current_bb_ids_ds, bb_hop_distance_matrix))
  
  ds_within_cap <- ds_reachable_with_steps %>%
    dplyr::filter(.data$dist_km <= threshold_downstream_km)
  
  chosen_downstream <- ds_within_cap %>%
    dplyr::group_by(.data$future_dam_id) %>%
    # dplyr::coalesce puts unknown hop (NA) after all finite hop tiers.
    dplyr::arrange(
      dplyr::coalesce(.data$trunk_step_dist, .Machine$integer.max),
      .data$dist_km
    ) %>%
    # dplyr::slice(1L) keeps the best row in each future dam group.
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = .data$future_dam_id,
      dam_id_down = .data$current_dam_id,
      min_distance_downstream_km = .data$dist_km,
      downstream_hops = .data$trunk_step_dist
    )
  
  downstream_summary <- ds_within_cap %>%
    dplyr::group_by(.data$future_dam_id) %>%
    dplyr::summarise(
      has_current_downstream = dplyr::n() > 0,
      min_distance_downstream_km = ifelse(dplyr::n() > 0, min(.data$dist_km), NA_real_),
      .groups = "drop"
    )
  
  # =============================================================================
  # SECTION 8 — upstream candidate pool + tiered winner
  #
  # INPUTS:
  #   - us_dam_all + threshold_upstream_km + bb_hop_distance_matrix
  # OUTPUTS:
  #   - us_reachable, us_reachable_with_steps, us_within_cap
  #   - nearest_upstream (one winner per future)
  # DEBUG OBJECTS:
  #   - us_reachable, us_reachable_with_steps, us_within_cap, nearest_upstream
  # =============================================================================
  
  us_reachable <- us_dam_all %>% # start from all upstream direction pairs
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0) # keep only valid directed paths
  
  future_bb_ids_us <- as.character(us_reachable$future_bb_id)
  current_bb_ids_us <- as.character(us_reachable$current_bb_id)
  us_reachable_with_steps <- us_reachable %>%
    dplyr::mutate(trunk_step_dist = lookup_bb_hops_for_pairs(future_bb_ids_us, current_bb_ids_us, bb_hop_distance_matrix))
  
  us_within_cap <- us_reachable_with_steps %>%
    dplyr::filter(.data$dist_km <= threshold_upstream_km)
  
  nearest_upstream <- us_within_cap %>%
    dplyr::group_by(.data$future_dam_id) %>%
    # Same tier rule as downstream: hops first, then shortest river km.
    dplyr::arrange(
      dplyr::coalesce(.data$trunk_step_dist, .Machine$integer.max),
      .data$dist_km
    ) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = .data$future_dam_id,
      dam_id_up = .data$current_dam_id,
      min_distance_upstream_km = .data$dist_km,
      upstream_hops = .data$trunk_step_dist
    )
  
  # =============================================================================
  # SECTION 9 — assemble future rows + classify connectivity
  #
  # INPUTS:
  #   - future_dam_ids + future_trunks + chosen_downstream + nearest_upstream
  # OUTPUTS:
  #   - reach_future (future-only output rows with categories)
  # DEBUG OBJECTS:
  #   - reach_future is reflected in decision_table later
  # =============================================================================
  
  reach_future <- data.frame( # initialize one output row per future dam id
    dam_id = future_dam_ids, # output key
    dam_type = rep("future", n_fut), # explicit type label for downstream filtering
    stringsAsFactors = FALSE # keep character columns as character, not factor
  ) %>%
    dplyr::left_join(
      dplyr::select(future_trunks, future_dam_id, bb_id),
      by = c("dam_id" = "future_dam_id")
    ) %>%
    dplyr::left_join(chosen_downstream, by = c("dam_id" = "future_dam_id")) %>%
    dplyr::left_join(nearest_upstream, by = c("dam_id" = "future_dam_id")) %>%
    dplyr::mutate(
      has_current_downstream = !is.na(.data$dam_id_down),
      has_current_upstream = !is.na(.data$dam_id_up),
      cascade_level = dplyr::if_else(
        !is.na(.data$dam_id_down) & !is.na(.data$dam_id_up) &
          !is.na(.data$downstream_hops) & !is.na(.data$upstream_hops),
        as.integer(pmax(.data$downstream_hops, .data$upstream_hops)),
        NA_integer_
      ),
      connectivity_category = dplyr::case_when(
        !is.na(.data$cascade_level) &
          .data$downstream_hops == 0L & .data$upstream_hops == 0L ~ "cascade_classic",
        !is.na(.data$cascade_level) & .data$cascade_level == 1L ~ "cascade1",
        !is.na(.data$cascade_level) & .data$cascade_level == 2L ~ "cascade2",
        !is.na(.data$cascade_level) & .data$cascade_level >= 3L ~ "cascade3+",
        (is.na(.data$min_distance_upstream_km) | .data$min_distance_upstream_km > threshold_upstream_km) &
          (is.na(.data$min_distance_downstream_km) | .data$min_distance_downstream_km > threshold_downstream_km) ~ "FFR",
        !is.na(.data$min_distance_upstream_km) &
          .data$min_distance_upstream_km <= threshold_upstream_km &
          (is.na(.data$min_distance_downstream_km) | .data$min_distance_downstream_km > threshold_downstream_km) ~ "downstream",
        !is.na(.data$min_distance_downstream_km) &
          .data$min_distance_downstream_km <= threshold_downstream_km &
          (is.na(.data$min_distance_upstream_km) | .data$min_distance_upstream_km > threshold_upstream_km) ~ "upstream",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(
      dam_id = .data$dam_id,
      dam_type = .data$dam_type,
      bb_id = .data$bb_id,
      has_current_downstream = .data$has_current_downstream,
      min_distance_downstream_km = .data$min_distance_downstream_km,
      dam_id_down = .data$dam_id_down,
      downstream_hops = .data$downstream_hops,
      has_current_upstream = .data$has_current_upstream,
      min_distance_upstream_km = .data$min_distance_upstream_km,
      dam_id_up = .data$dam_id_up,
      upstream_hops = .data$upstream_hops,
      cascade_level = .data$cascade_level,
      connectivity_category = .data$connectivity_category
    )
  
  # =============================================================================
  # SECTION 10 — append current dam placeholder rows
  #
  # INPUTS:
  #   - current_dam_ids + current_trunks
  # OUTPUTS:
  #   - reach_current, then reach_df = bind_rows(reach_future, reach_current)
  # DEBUG OBJECTS:
  #   - reach_df returned; current placeholders are not added separately to debug
  # =============================================================================
  
  reach_current <- data.frame( # mirror output schema for current dams
    dam_id = current_dam_ids, # current dam ids from node table
    dam_type = rep("current", n_cur), # mark these rows as current
    bb_id = current_trunks$bb_id, # attach trunk id where available
    has_current_downstream = rep(NA, n_cur), # not defined for current reference rows
    min_distance_downstream_km = rep(NA_real_, n_cur), # placeholder numeric
    dam_id_down = rep(NA_character_, n_cur), # placeholder id
    downstream_hops = rep(NA_integer_, n_cur), # placeholder hops
    has_current_upstream = rep(NA, n_cur), # placeholder logical
    min_distance_upstream_km = rep(NA_real_, n_cur), # placeholder numeric
    dam_id_up = rep(NA_character_, n_cur), # placeholder id
    upstream_hops = rep(NA_integer_, n_cur), # placeholder hops
    cascade_level = rep(NA_integer_, n_cur), # placeholder cascade level
    connectivity_category = rep(NA_character_, n_cur), # placeholder category
    stringsAsFactors = FALSE # keep placeholders as character vectors
  )
  
  reach_df <- dplyr::bind_rows(reach_future, reach_current)
  
  # =============================================================================
  # SECTION 11 — build decision_table audit sheet
  #
  # INPUTS:
  #   - future ids + trunk lookups + downstream/upstream picks + final categories
  # OUTPUTS:
  #   - decision_table (one row per future dam, for QA and debugging)
  # DEBUG OBJECTS:
  #   - decision_table
  # =============================================================================
  
  decision_table <- data.frame(future_dam_id = future_dam_ids, stringsAsFactors = FALSE) %>% # start with one row per future dam id
    dplyr::left_join(
      dplyr::rename(future_trunks, future_dam_id = .data$future_dam_id),
      by = "future_dam_id"
    ) %>%
    dplyr::left_join(
      dplyr::rename(
        downstream_summary,
        has_ds_in_cap = .data$has_current_downstream,
        min_ds_km_in_cap_pool = .data$min_distance_downstream_km
      ),
      by = "future_dam_id"
    ) %>%
    dplyr::left_join(chosen_downstream, by = "future_dam_id") %>%
    dplyr::left_join(nearest_upstream, by = "future_dam_id") %>%
    dplyr::left_join(
      dplyr::select(reach_future, dam_id, cascade_level, connectivity_category),
      by = c("future_dam_id" = "dam_id")
    )
  
  # =============================================================================
  # SECTION 12 — collect debug list objects
  #
  # INPUTS:
  #   - all intermediate objects created above
  # OUTPUTS:
  #   - debug named list included in return value
  # DEBUG OBJECTS:
  #   - (this section defines the full debug payload)
  # =============================================================================
  
  debug <- list( # collect all intermediate objects for inspection
    nodes_tbl = nodes_tbl, # flattened node attributes with node_id
    edges_tbl = edges_tbl, # flattened edge attributes
    weights_km = weights_km, # numeric edge weight vector used in igraph
    dam_nodes_tbl = dam_nodes_tbl, # node subset where dam_id exists
    current_dams_tbl = current_dams_tbl, # current dam node subset
    future_dams_tbl = future_dams_tbl, # future dam node subset
    node_to_bb_lookup = node_to_bb_lookup, # node -> bb_id lookup table
    ds_mat = ds_mat, # future->current distance matrix (km)
    us_mat = us_mat, # current->future distance matrix (km)
    ds_grid = ds_grid, # index grid used to melt ds_mat
    us_grid = us_grid, # index grid used to melt us_mat
    ds_dam_all = ds_dam_all, # all downstream direction pair rows
    us_dam_all = us_dam_all, # all upstream direction pair rows
    ds_reachable = ds_reachable, # downstream pairs with finite positive distance
    ds_same_trunk = ds_same_trunk, # downstream reachable pairs where bb_id matches
    ds_reachable_with_steps = ds_reachable_with_steps, # downstream pairs with trunk hops
    ds_within_cap = ds_within_cap, # downstream pairs within km threshold
    chosen_downstream = chosen_downstream, # downstream winner per future dam
    downstream_summary = downstream_summary, # downstream pool summary per future dam
    bb_confluence_edges = bb_confluence_edges, # bb_id adjacency edges
    bb_confluence_graph = bb_confluence_graph, # igraph built from bb adjacency
    bb_hop_distance_matrix = bb_hop_distance_matrix, # all-pairs bb hop distances
    us_reachable = us_reachable, # upstream pairs with finite positive distance
    us_reachable_with_steps = us_reachable_with_steps, # upstream pairs with trunk hops
    us_within_cap = us_within_cap, # upstream pairs within km threshold
    nearest_upstream = nearest_upstream, # upstream winner per future dam
    decision_table = decision_table, # wide QA table for final decisions
    future_trunks = future_trunks, # future dam -> bb lookup
    current_trunks = current_trunks # current dam -> bb lookup
  )
  
  list(reach_df = reach_df, debug = debug, threshold_used = threshold_used)
}


# -----------------------------------------------------------------------------
# empty_reach_connectivity_df
# -----------------------------------------------------------------------------

empty_reach_connectivity_df <- function() {
  data.frame(
    dam_id = character(),
    dam_type = character(),
    bb_id = character(),
    has_current_downstream = logical(),
    min_distance_downstream_km = numeric(),
    dam_id_down = character(),
    downstream_hops = integer(),
    has_current_upstream = logical(),
    min_distance_upstream_km = numeric(),
    dam_id_up = character(),
    upstream_hops = integer(),
    cascade_level = integer(),
    connectivity_category = character(),
    stringsAsFactors = FALSE
  )
}
