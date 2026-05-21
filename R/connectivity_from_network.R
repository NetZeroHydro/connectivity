# =============================================================================
# connectivity_from_network.R
# =============================================================================
# Connectivity for a blended dam network: builds `reach_df` (one row per
# dam) plus a `debug` list of tibbles.
#
#
# Prerequisites: directed `sfnetwork` with edge `weight` (metres) and `bb_id`;
# blended nodes carry `dam_id` and `is_current_dam` (TRUE = current, FALSE =
# future). Load: sfnetworks, igraph, dplyr, sf, tidygraph before sourcing.
#
# =============================================================================
#'
#' Builds directed river distances between **current** and **future** dam nodes
#' on a blended `sfnetwork`, applies tiered upstream/downstream rules with km
#' caps, trunk hop counts, `cascade_level`, and `connectivity_category`.
#'
#' **Dam nodes:** nodes with non-missing `dam_id`; `is_current_dam` is `TRUE` for
#' current dams and `FALSE` for future dams.
#'
#' **Edge weights:** column `weight` is treated as **metres** and converted to
#' **km** for [igraph::distances()]. Edges must
#' carry trunk id `bb_id` for the trunk confluence graph.
#'
#' **Downstream neighbor:** finite future→current path where `dist_km > 0` and at most
#' `threshold_downstream_km`; pick **fewest trunk hops** then smallest `dist_km`.
#'
#' **Upstream neighbor:** finite current→future path, `dist_km > 0` and at most
#' `threshold_upstream_km`; same tier rule.
#'
#' **`cascade_level`:** `pmax(downstream_hops, upstream_hops)` when both
#' neighbors exist, both hop counts are non-NA, and both neighbor distances are
#' at most `threshold_cascade_km`; else `NA_integer_`.
#'
#' **`connectivity_category` cascade levels:** `cascade_classic` when both
#' neighbors are on the same trunk (`downstream_hops == 0` and
#' `upstream_hops == 0`), then `cascade1`, `cascade2+` for max hops
#' 1, 2+ respectively — only when both `min_distance_*_km` are within
#' `threshold_cascade_km` (independent of up/down search caps).
#'
#' **`threshold_cascade_km`:** km cap for cascade classification only; does not
#' change which up/down neighbors are chosen.
#'
#' @param net_with_dams Directed `sfnetwork` with edge `bb_id`, numeric edge
#'   `weight` (metres), and blended node columns `dam_id` / `is_current_dam`.
#' @param threshold_downstream_km Required km cap on downstream candidate paths.
#' @param threshold_upstream_km Required km cap on upstream candidate paths.
#' @param threshold_cascade_km Required cascade threshold in km.
#'
#' @return Named list: `reach_df`, `debug`, and `threshold_used`.
#'
#' @details
#' Distances use `igraph::distances(..., mode = "out")`. Excludes: `Inf` = no path; `0` =
#' same vertex (self)
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
  
  # Store thresholds for output
  
  threshold_used <- c(
    threshold_downstream_km = threshold_downstream_km,
    threshold_upstream_km = threshold_upstream_km,
    threshold_cascade_km = threshold_cascade_km
  )
  
  # =============================================================================
  # SECTION 1 — Convert sfnetwork to node/edge tables
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
    dplyr::mutate(node_id = dplyr::row_number()) # create node index used in pair tables
  
  edges_tbl <- net_with_dams %>% # same process as above but with edges
    tidygraph::activate("edges") %>%
    sf::st_as_sf() %>%
    sf::st_drop_geometry() %>%
    dplyr::as_tibble()
  
  weights_km <- as.numeric(edges_tbl$weight) / 1000 # convert m to km
  
  # node edge starts from
  edge_from_node_bb <- edges_tbl %>% dplyr::select(node_id = "from", bb_id = "bb_id")
  # node edge ends at
  edge_to_node_bb <- edges_tbl %>% dplyr::select(node_id = "to", bb_id = "bb_id")
  # all edge confluences
  edge_node_bb_incidence <- dplyr::bind_rows(edge_from_node_bb, edge_to_node_bb)
  
  # One bb_id per node used for early exits and pair tables.
  node_to_bb_lookup <- edge_node_bb_incidence %>%
    dplyr::filter(!is.na(bb_id)) %>%
    dplyr::group_by(node_id) %>%
    dplyr::summarise(bb_id = dplyr::first(bb_id), .groups = "drop")
  
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
    dplyr::filter(!is.na(dam_id)) %>% # keep only snapped dam nodes (non-missing id)
    dplyr::mutate(dam_id = as.character(dam_id)) # force character id for stable joins
  
  current_dams_tbl <- dam_nodes_tbl %>% # current dams table
    dplyr::filter(is_current_dam %in% TRUE) %>%
    dplyr::select(node_id, dam_id, is_current_dam)
  
  future_dams_tbl <- dam_nodes_tbl %>% # future dams table
    dplyr::filter(is_current_dam %in% FALSE) %>%
    dplyr::select(node_id, dam_id, is_current_dam)
  
  current_nodes <- current_dams_tbl$node_id # igraph indices
  future_nodes <- future_dams_tbl$node_id
  current_dam_ids <- current_dams_tbl$dam_id # dam labels
  future_dam_ids <- future_dams_tbl$dam_id
  n_cur <- length(current_nodes) # lengths
  n_fut <- length(future_nodes)
  
  # Empty dam data sets stop check. 
  # Stops the function call if missing either 
  # current, future or both dams. 
  # Returns empty reach_df instead of failing.
  
  if (n_fut == 0L) {
    return(list(
      reach_df = empty_reach_connectivity_df(),
      debug = list(note = "no future dam nodes on network"),
      threshold_used = threshold_used
    ))
  }
  
  if (n_cur == 0L) {
    future_trunks <- future_dams_tbl %>%
      dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
      dplyr::transmute(
        future_dam_id = as.character(dam_id),
        future_node_id = node_id,
        bb_id = bb_id
      )
    
    reach_future <- data.frame(
      dam_id = future_dam_ids,
      dam_type = rep("future", n_fut),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(
        dplyr::transmute(future_trunks, dam_id = future_dam_id, bb_id = bb_id),
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
  
  us_mat <- igraph::distances( # same process but for upstream
    net_with_dams,
    v = current_nodes,
    to = future_nodes,
    weights = weights_km,
    mode = "out"
  )
  
  # =============================================================================
  # SECTION 4 — convert matrices into pair tables
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
  
  ds_dam_all <- ds_grid %>% # Every future dam paired with every current dam, with downstream river distance in km
    dplyr::mutate(
      future_node_id = future_nodes[fut_idx],
      current_node_id = current_nodes[cur_idx],
      future_dam_id = future_dam_ids[fut_idx],
      current_dam_id = current_dam_ids[cur_idx]
    ) %>%
    dplyr::select(
      future_dam_id = future_dam_id,
      current_dam_id = current_dam_id,
      future_node_id = future_node_id,
      current_node_id = current_node_id,
      dist_km = dist_km
    )
  
  us_grid <- expand.grid(cur_idx = seq_len(n_cur), fut_idx = seq_len(n_fut)) # same process but for upstream
  us_grid$dist_km <- as.vector(us_mat)
  
  us_dam_all <- us_grid %>% # Every future dam paired with every current dam, with upstream river distance in km
    dplyr::mutate(
      future_node_id = future_nodes[fut_idx],
      current_node_id = current_nodes[cur_idx],
      future_dam_id = future_dam_ids[fut_idx],
      current_dam_id = current_dam_ids[cur_idx]
    ) %>%
    dplyr::select(
      future_dam_id = future_dam_id,
      current_dam_id = current_dam_id,
      future_node_id = future_node_id,
      current_node_id = current_node_id,
      dist_km = dist_km
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
  
  ds_dam_all <- ds_dam_all %>% # join bb_id to all downstream dam pairs
    dplyr::left_join(node_to_bb_lookup, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = bb_id) %>%
    dplyr::left_join(node_to_bb_lookup, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = bb_id)
  
  us_dam_all <- us_dam_all %>% # join bb_id to all upstream dam pairs
    dplyr::left_join(node_to_bb_lookup, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = bb_id) %>%
    dplyr::left_join(node_to_bb_lookup, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = bb_id)
  
  future_trunks <- future_dams_tbl %>% # join bb_id to all future dams
    dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
    dplyr::rename(future_node_id = node_id, future_dam_id = dam_id, bb_id = bb_id) %>%
    dplyr::select(future_dam_id, future_node_id, bb_id)
  
  current_trunks <- current_dams_tbl %>% # join bb_id to all current dams
    dplyr::left_join(node_to_bb_lookup, by = c("node_id" = "node_id")) %>%
    dplyr::rename(current_node_id = node_id, current_dam_id = dam_id, bb_id = bb_id) %>%
    dplyr::select(current_dam_id, current_node_id, bb_id)
  
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
  bb_confluence_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE) # all bb_id-to-bb_id links
  
  for (junction_node_id in junction_nodes) { # inspect each graph junction separately
    bb_ids_at_node <- edge_node_bb_incidence %>% dplyr::filter(node_id == junction_node_id) %>% dplyr::pull(bb_id) # collect incident trunk ids
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
  
  bb_confluence_edges <- unique(bb_confluence_edges)  # list of touching trunks
  bb_confluence_graph <- igraph::graph_from_data_frame(bb_confluence_edges, directed = FALSE) # graph of touching trunks
  bb_hop_distance_matrix <- igraph::distances(bb_confluence_graph, mode = "all") # hop distance matrix
  if (igraph::vcount(bb_confluence_graph) > 0) {
    bb_vertex_names <- igraph::V(bb_confluence_graph)$name
    dimnames(bb_hop_distance_matrix) <- list(bb_vertex_names, bb_vertex_names)
  }
  
  # Helper function looks up hops between all dams
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
    dplyr::filter(is.finite(dist_km) & dist_km > 0) # keep only real reachable paths (exclude Inf/self 0)
  
  ds_same_trunk <- ds_reachable %>% # get same trunk pairs
    dplyr::filter(current_bb_id == future_bb_id)
  
  future_bb_ids_ds <- as.character(ds_reachable$future_bb_id) # future downstream trunks
  current_bb_ids_ds <- as.character(ds_reachable$current_bb_id) # current downstream trunks
  ds_reachable_with_steps <- ds_reachable %>% # use helper function for hops downstream
    dplyr::mutate(trunk_step_dist = lookup_bb_hops_for_pairs(future_bb_ids_ds, current_bb_ids_ds, bb_hop_distance_matrix))
  
  ds_within_cap <- ds_reachable_with_steps %>% # filter using threshold
    dplyr::filter(dist_km <= threshold_downstream_km)
  
  chosen_downstream <- ds_within_cap %>% # chooses nearest neighbor
    dplyr::group_by(future_dam_id) %>%
    dplyr::arrange(
      dplyr::coalesce(trunk_step_dist, .Machine$integer.max), # puts unknown hop (NA) after all finite hop tiers
dist_km
    ) %>%
    dplyr::slice(1L) %>% # keep the best row in each future dam group
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = future_dam_id,
      dam_id_down = current_dam_id,
      min_distance_downstream_km = dist_km,
      downstream_hops = trunk_step_dist
    )
  
  downstream_summary <- ds_within_cap %>% # debug output
    dplyr::group_by(future_dam_id) %>%
    dplyr::summarise(
      has_current_downstream = dplyr::n() > 0,
      min_distance_downstream_km = ifelse(dplyr::n() > 0, min(dist_km), NA_real_),
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
    dplyr::filter(is.finite(dist_km) & dist_km > 0) # keep only valid directed paths
  
  future_bb_ids_us <- as.character(us_reachable$future_bb_id) # future upstream trunks
  current_bb_ids_us <- as.character(us_reachable$current_bb_id) # current upstream trunks
  us_reachable_with_steps <- us_reachable %>% # use helper function for hops upstream
    dplyr::mutate(trunk_step_dist = lookup_bb_hops_for_pairs(future_bb_ids_us, current_bb_ids_us, bb_hop_distance_matrix))
  
  us_within_cap <- us_reachable_with_steps %>% # filter for within threshold
    dplyr::filter(dist_km <= threshold_upstream_km)
  
  nearest_upstream <- us_within_cap %>%
    dplyr::group_by(future_dam_id) %>% # Same tier rule as downstream: hops first, then shortest river km.
    dplyr::arrange(
      dplyr::coalesce(trunk_step_dist, .Machine$integer.max),
dist_km
    ) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = future_dam_id,
      dam_id_up = current_dam_id,
      min_distance_upstream_km = dist_km,
      upstream_hops = trunk_step_dist
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
      has_current_downstream = !is.na(dam_id_down),
      has_current_upstream = !is.na(dam_id_up),
      within_cascade_km = !is.na(dam_id_down) &
        !is.na(dam_id_up) &
        !is.na(min_distance_upstream_km) &
        !is.na(min_distance_downstream_km) &
min_distance_upstream_km <= threshold_cascade_km &
min_distance_downstream_km <= threshold_cascade_km,
      cascade_level = dplyr::if_else(
within_cascade_km &
          !is.na(downstream_hops) & !is.na(upstream_hops),
        as.integer(pmax(downstream_hops, upstream_hops)),
        NA_integer_
      ),
      connectivity_category = dplyr::case_when( # define connectivity categories
within_cascade_km &
          !is.na(cascade_level) &
downstream_hops == 0L & upstream_hops == 0L ~ "cascade_classic",
within_cascade_km &
          !is.na(cascade_level) & cascade_level == 1L ~ "cascade1",
within_cascade_km &
          !is.na(cascade_level) & cascade_level >= 2L ~ "cascade2+",
        #!is.na(cascade_level) & cascade_level >= 3L ~ "cascade3+", edit here if more cascade levels are desired
        (is.na(min_distance_upstream_km) | min_distance_upstream_km > threshold_upstream_km) &
          (is.na(min_distance_downstream_km) | min_distance_downstream_km > threshold_downstream_km) ~ "undammed",
        !is.na(min_distance_upstream_km) &
min_distance_upstream_km <= threshold_upstream_km &
          (is.na(min_distance_downstream_km) | min_distance_downstream_km > threshold_downstream_km) ~ "downstream",
        !is.na(min_distance_downstream_km) &
min_distance_downstream_km <= threshold_downstream_km &
          (is.na(min_distance_upstream_km) | min_distance_upstream_km > threshold_upstream_km) ~ "upstream",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select( # finalize output columns
      dam_id = dam_id,
      dam_type = dam_type,
      bb_id = bb_id,
      has_current_downstream = has_current_downstream,
      min_distance_downstream_km = min_distance_downstream_km,
      dam_id_down = dam_id_down,
      downstream_hops = downstream_hops,
      has_current_upstream = has_current_upstream,
      min_distance_upstream_km = min_distance_upstream_km,
      dam_id_up = dam_id_up,
      upstream_hops = upstream_hops,
      cascade_level = cascade_level,
      connectivity_category = connectivity_category
    )
  
  # =============================================================================
  # SECTION 10 — add current dam rows
  #
  # INPUTS:
  #   - current_dam_ids + current_trunks
  # OUTPUTS:
  #   - reach_current, then reach_df = bind_rows(reach_future, reach_current)
  # DEBUG OBJECTS:
  #   - reach_df returned
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
  
  reach_df <- dplyr::bind_rows(reach_future, reach_current) # final dataframe
  
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
      dplyr::rename(future_trunks, future_dam_id = future_dam_id),
      by = "future_dam_id"
    ) %>%
    dplyr::left_join(
      dplyr::rename(
        downstream_summary,
        has_ds_in_cap = has_current_downstream,
        min_ds_km_in_cap_pool = min_distance_downstream_km
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
