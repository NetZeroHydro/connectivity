# =============================================================================
# connectivity_function.R
# =============================================================================
# Build connectivity matrices and reach-level dataframe from a river network
# with dams already snapped. Use this for
# a single basin or sub-basin.
#
# Expects: sfnetwork with nodes that have dam_id and is_current_dam when dams
# were blended via st_network_blend; otherwise returns an empty template.
#
# Choose your study area outside this function: crop rivers and dams to your
# basin or sub-basin polygon, then as_sfnetwork() + st_network_blend(), then
# pass that net_with_dams here. Or use net_with_dams_from_hydrobasin() in
# basin_network.R (in-memory loads; rivers: continent + Strahler, then basin crop) to build
# net_with_dams, then add_ffr_attr() after this function for FFR columns.
#
# Reach-level connectivity uses trunk (bb_id) scope: cascade_trunk_steps = 0 is
# same-trunk only; steps > 0 adds trunks within that many edges on a graph built
# from confluences (network nodes where two or more distinct bb_id meet).
#
# Required packages: sfnetworks, igraph, dplyr, sf (load before sourcing or calling
# connectivity_from_network). Edge `weight` should be length in metres; values are
# divided by 1000 for km after as.numeric().
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
#' @param net_with_dams An `sfnetwork` created from rivers (`directed = TRUE`)
#'   and optionally blended with dam points. When blended, node columns include
#'   `dam_id` and `is_current_dam` (TRUE = current, FALSE = future). Edge column
#'   `weight` is typically metres from `edge_length()` (converted with
#'   `as.numeric()` then divided by 1000 for km). Edge column `bb_id` (trunk ID)
#'   is required. If there are no dam nodes, returns a zero-row data frame with
#'   the standard columns.
#' @param cascade_threshold_km Optional numeric. If set, a future dam is
#'   `"cascade"` only if the nearest-current upstream *and* downstream distances
#'   (among allowed trunks) are both `<= cascade_threshold_km` (km). If NULL,
#'   any upstream + any downstream among those trunks counts as cascade.
#' @param cascade_trunk_steps Non-negative integer. `0` = only current dams on
#'   the **same** `bb_id` as the future dam are considered. `1, 2, ...` = also
#'   current dams on trunks within that many steps on an **undirected trunk graph**
#'   inferred from the network: at each graph node, if two or more distinct
#'   non-missing `bb_id` appear on incident edges, those trunks are linked.
#'   Requires confluences to appear as single nodes with correct `bb_id` on edges.
#'
#' @return A data frame with one row per dam node: `dam_id`, `dam_type`
#'   (`"future"` or `"current"`), `has_current_downstream`, `has_current_upstream`,
#'   `min_distance_downstream_km`, `min_distance_upstream_km`, `cascade_status`,
#'   `bb_id`. Current-dam rows have `NA` in connectivity and cascade columns;
#'   future-dam rows are sorted by cascade status and distances, then current
#'   rows are appended.
#'
#' @details Distances are computed with `igraph::distances(..., mode = "out")`
#'   on the full directed river graph. `Inf` means no path; `0` means the same
#'   node. Allowed current dams for each future dam are those whose trunk lies
#'   within `cascade_trunk_steps` steps of the future dam’s trunk on the inferred
#'   trunk graph (or exactly the same trunk when `cascade_trunk_steps = 0`).
#'
#' @examples
#' \dontrun{
#' reach_df <- connectivity_from_network(net_with_dams)
#'
#' out <- connectivity_from_network(net_with_dams, cascade_threshold_km = 15)
#'
#' reach_wider <- connectivity_from_network(net_with_dams, cascade_trunk_steps = 1,
#'                                          cascade_threshold_km = 500)
#' }
#'
#' @export
connectivity_from_network <- function(net_with_dams,
                                      cascade_threshold_km = 100,
                                      cascade_trunk_steps = NULL) {
  
  # =============================================================================
  # SECTION 0 — what this function returns
  # =============================================================================
  
  # This function returns a list with two items.
  # - reach_df: the final per-dam output table
  # - debug: a list of tables showing every step of the calculation
  
  # =============================================================================
  # SECTION 1 — nodes and edges as plain tables
  # =============================================================================
  
  # Pull nodes into a tibble so we can filter and join easily.
  nodes_tbl <- net_with_dams %>%
    tidygraph::activate("nodes") %>% # activate node table
    sf::st_as_sf() %>% # convert to sf so geometry can be dropped safely
    sf::st_drop_geometry() %>% # remove point geometry
    dplyr::as_tibble() %>% # convert to tibble
    dplyr::mutate(node_id = dplyr::row_number()) # store igraph node id
  
  # Pull edges into a tibble so we can read bb_id and weight.
  edges_tbl <- net_with_dams %>%
    tidygraph::activate("edges") %>% # activate edge table
    sf::st_as_sf() %>% # convert to sf
    sf::st_drop_geometry() %>% # remove line geometry
    dplyr::as_tibble() # convert to tibble
  
  # Convert edge weights from metres to kilometres (vector aligns with edges).
  weights_km <- as.numeric(edges_tbl$weight) / 1000 # weights in km
  
  # =============================================================================
  # SECTION 2 — identify dam nodes (current vs future)
  # =============================================================================
  
  # A dam node is a node where dam_id is NOT NA (these are created by blending).
  dam_nodes_tbl <- nodes_tbl %>%
    dplyr::filter(!is.na(.data$dam_id)) %>% # keep only dam nodes
    dplyr::mutate(dam_id = as.character(.data$dam_id)) # force dam_id to character
  
  # Current dams are dam nodes flagged TRUE.
  current_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% TRUE) %>% # keep current dams only
    dplyr::select(node_id, dam_id, is_current_dam) # keep only key columns
  
  # Future dams are dam nodes flagged FALSE.
  future_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% FALSE) %>% # keep future dams only
    dplyr::select(node_id, dam_id, is_current_dam) # keep only key columns
  
  # Extract node id vectors (these are the igraph vertex ids).
  current_nodes <- current_dams_tbl$node_id # current dam node ids
  
  # Extract node id vectors (these are the igraph vertex ids).
  future_nodes <- future_dams_tbl$node_id # future dam node ids
  
  # Extract dam id vectors (these are the dam identifiers).
  current_dam_ids <- current_dams_tbl$dam_id # current dam ids
  
  # Extract dam id vectors (these are the dam identifiers).
  future_dam_ids <- future_dams_tbl$dam_id # future dam ids
  
  # Count current dams (used to build grids).
  n_cur <- length(current_nodes) # number of current dams
  
  # Count future dams (used to build grids).
  n_fut <- length(future_nodes) # number of future dams
  
  # =============================================================================
  # SECTION 3 — distance matrices (km) for every dam pair
  # =============================================================================
  
  # Downstream matrix: distance from each future dam to each current dam (future -> current).
  ds_mat <- igraph::distances(
    net_with_dams, # the river graph
    v = future_nodes, # start at future dams
    to = current_nodes, # end at current dams
    weights = weights_km, # use edge length in km
    mode = "out" # follow edge direction downstream
  )
  
  # Upstream matrix: distance from each current dam to each future dam (current -> future).
  us_mat <- igraph::distances(
    net_with_dams, # the river graph
    v = current_nodes, # start at current dams
    to = future_nodes, # end at future dams
    weights = weights_km, # use edge length in km
    mode = "out" # follow edge direction downstream (current to future means current is upstream)
  )
  
  # =============================================================================
  # SECTION 4 — build ds_dam and us_dam tables (ALL dam pairs)
  # =============================================================================
  
  # Create a full grid of indices for downstream pairs (future index, current index).
  ds_grid <- expand.grid( # create all future/current combinations
    fut_idx = seq_len(n_fut), # 1..n_fut
    cur_idx = seq_len(n_cur) # 1..n_cur
  )
  
  # Flatten the downstream matrix into a vector that matches expand.grid row order.
  ds_grid$dist_km <- as.vector(ds_mat) # downstream distance values
  
  # Map indices to ids/node_ids to create a readable downstream dam-pair table.
  ds_dam_all <- ds_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx], # future node id
      current_node_id = current_nodes[.data$cur_idx], # current node id
      future_dam_id = future_dam_ids[.data$fut_idx], # future dam id
      current_dam_id = current_dam_ids[.data$cur_idx] # current dam id
    ) %>%
    dplyr::select(
      .data$future_dam_id, # future dam id
      .data$current_dam_id, # current dam id
      .data$future_node_id, # future node id
      .data$current_node_id, # current node id
      dist_km = .data$dist_km # downstream distance (future -> current)
    )
  
  # Create a full grid of indices for upstream pairs (future index, current index).
  us_grid <- expand.grid( # create all future/current combinations
    cur_idx = seq_len(n_cur), # 1..n_cur
    fut_idx = seq_len(n_fut) # 1..n_fut
  )
  
  # Flatten the upstream matrix into a vector that matches expand.grid row order.
  us_grid$dist_km <- as.vector(us_mat) # upstream distance values
  
  # Map indices to ids/node_ids to create a readable upstream dam-pair table.
  us_dam_all <- us_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx], # future node id
      current_node_id = current_nodes[.data$cur_idx], # current node id
      future_dam_id = future_dam_ids[.data$fut_idx], # future dam id
      current_dam_id = current_dam_ids[.data$cur_idx] # current dam id
    ) %>%
    dplyr::select(
      .data$future_dam_id, # future dam id
      .data$current_dam_id, # current dam id
      .data$future_node_id, # future node id
      .data$current_node_id, # current node id
      dist_km = .data$dist_km # upstream distance (current -> future)
    )
  
  # =============================================================================
  # SECTION 5 — bb_finder: node -> bb_id mapping and add bb_id to pair tables
  # =============================================================================
  
  # Build an incident table for edge endpoints (from side).
  inc_from <- edges_tbl %>% dplyr::select(node_id = .data$from, bb_id = .data$bb_id) # from-endpoint bb_id
  
  # Build an incident table for edge endpoints (to side).
  inc_to <- edges_tbl %>% dplyr::select(node_id = .data$to, bb_id = .data$bb_id) # to-endpoint bb_id
  
  # Combine incident rows so each node sees all bb_id values touching it.
  inc_all <- dplyr::bind_rows(inc_from, inc_to) # all incident node/bb_id rows
  
  # Collapse incident bb_id values into one bb_id per node (first one is enough for dams).
  node_bb <- inc_all %>%
    dplyr::filter(!is.na(.data$bb_id)) %>% # drop missing trunk ids
    dplyr::group_by(.data$node_id) %>% # group by node
    dplyr::summarise(bb_id = dplyr::first(.data$bb_id), .groups = "drop") # choose the first bb_id
  
  # Add bb_id for the future and current node to every downstream pair row.
  ds_dam_all <- ds_dam_all %>%
    dplyr::left_join(node_bb, by = c("future_node_id" = "node_id")) %>% # join future node bb_id
    dplyr::rename(future_bb_id = .data$bb_id) %>% # rename future bb_id
    dplyr::left_join(node_bb, by = c("current_node_id" = "node_id")) %>% # join current node bb_id
    dplyr::rename(current_bb_id = .data$bb_id) # rename current bb_id
  
  # Add bb_id for the future and current node to every upstream pair row.
  us_dam_all <- us_dam_all %>%
    dplyr::left_join(node_bb, by = c("future_node_id" = "node_id")) %>% # join future node bb_id
    dplyr::rename(future_bb_id = .data$bb_id) %>% # rename future bb_id
    dplyr::left_join(node_bb, by = c("current_node_id" = "node_id")) %>% # join current node bb_id
    dplyr::rename(current_bb_id = .data$bb_id) # rename current bb_id
  
  # Build a future dam trunk table (used for final output).
  future_trunks <- future_dams_tbl %>%
    dplyr::left_join(node_bb, by = c("node_id" = "node_id")) %>% # attach bb_id
    dplyr::rename(future_node_id = .data$node_id, future_dam_id = .data$dam_id, bb_id = .data$bb_id) %>% # rename columns
    dplyr::select(.data$future_dam_id, .data$future_node_id, .data$bb_id) # keep columns
  
  # Build a current dam trunk table (used for appending current rows).
  current_trunks <- current_dams_tbl %>%
    dplyr::left_join(node_bb, by = c("node_id" = "node_id")) %>% # attach bb_id
    dplyr::rename(current_node_id = .data$node_id, current_dam_id = .data$dam_id, bb_id = .data$bb_id) %>% # rename columns
    dplyr::select(.data$current_dam_id, .data$current_node_id, .data$bb_id) # keep columns
  
  # =============================================================================
  # SECTION 6 — ds_filter: strict downstream SAME-TRUNK candidates + chosen downstream
  # =============================================================================
  
  # Keep only real downstream connections (finite, > 0).
  ds_reachable <- ds_dam_all %>%
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0) # downstream reachable pairs
  
  # Keep only downstream connections where the current dam is on the SAME trunk as the future dam.
  ds_same_trunk <- ds_reachable %>%
    dplyr::filter(.data$current_bb_id == .data$future_bb_id) # strict same-trunk downstream
  
  # For each future dam, choose the nearest downstream current dam on the same trunk.
  chosen_downstream <- ds_same_trunk %>%
    dplyr::group_by(.data$future_dam_id) %>% # group by future dam
    dplyr::arrange(.data$dist_km) %>% # sort by smallest distance
    dplyr::slice(1L) %>% # keep the closest row
    dplyr::ungroup() %>% # remove grouping
    dplyr::rename( # rename columns to show these are the chosen downstream values
      chosen_ds_current_dam_id = .data$current_dam_id, # chosen downstream current dam
      chosen_ds_dist_km = .data$dist_km # chosen downstream distance
    ) %>%
    dplyr::select(.data$future_dam_id, .data$chosen_ds_current_dam_id, .data$chosen_ds_dist_km) # keep only key columns
  
  # Build a per-future downstream summary table (flag + minimum distance).
  downstream_summary <- ds_same_trunk %>%
    dplyr::group_by(.data$future_dam_id) %>% # group by future dam
    dplyr::summarise(
      has_current_downstream = dplyr::n() > 0, # TRUE if any same-trunk downstream current exists
      min_distance_downstream_km = ifelse(dplyr::n() > 0, min(.data$dist_km), NA_real_), # minimum downstream distance
      .groups = "drop" # drop grouping
    )
  
  # =============================================================================
  # SECTION 7 — trunk graph: compute trunk-step distance for upstream search
  # =============================================================================
  
  # Collect all node ids in the network (potential confluences).
  junction_nodes <- unique(c(edges_tbl$from, edges_tbl$to)) # all junction nodes
  
  # Prepare an empty edge list for trunk adjacency (filled below).
  trunk_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE) # trunk adjacency list
  
  # Loop through junction nodes and connect trunks that meet at the node.
  for (nid in junction_nodes) { # for each node
    bb_here <- inc_all %>% dplyr::filter(.data$node_id == nid) %>% dplyr::pull(.data$bb_id) # all bb_ids at node
    bb_here <- unique(as.character(bb_here)) # unique trunk ids at node
    bb_here <- bb_here[!is.na(bb_here) & nzchar(bb_here)] # drop NA/empty
    if (length(bb_here) >= 2) { # confluence if 2+ trunks meet
      pairs <- utils::combn(bb_here, 2) # all trunk pairs
      trunk_edges <- rbind( # add to edge list
        trunk_edges, # existing rows
        data.frame(from = as.character(pairs[1, ]), to = as.character(pairs[2, ]), stringsAsFactors = FALSE) # new rows
      )
    }
  }
  
  # Remove duplicate trunk edges.
  trunk_edges <- unique(trunk_edges) # unique edges
  
  # Build an undirected graph where vertices are trunks and edges connect neighboring trunks.
  g_trunk <- igraph::graph_from_data_frame(
    trunk_edges, # trunk edge list
    directed = FALSE # undirected adjacency
  )
  
  # =============================================================================
  # SECTION 8 — us_filter: upstream candidates within threshold + trunk-step distance
  # =============================================================================
  
  # Keep only real upstream connections (finite, > 0).
  us_reachable <- us_dam_all %>%
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0) # upstream reachable pairs
  
  # Apply the distance threshold (default 100 km) to upstream connections.
  us_within_threshold <- us_reachable %>%
    dplyr::filter(.data$dist_km <= cascade_threshold_km) # upstream within threshold
  
  # Prepare a trunk distance matrix between all trunks in the trunk graph.
  trunk_dist_mat <- igraph::distances(g_trunk, mode = "all") # trunk-to-trunk hop distances
  
  # Give names to the trunk distance matrix so we can index by bb_id values.
  if (igraph::vcount(g_trunk) > 0) { # only if graph has vertices
    vnames <- igraph::V(g_trunk)$name # trunk ids used as vertex names
    dimnames(trunk_dist_mat) <- list(vnames, vnames) # label rows/cols by trunk id
  }
  
  # Add the trunk-step distance for each upstream pair (future trunk -> current trunk).
  # Extract the two trunk-id columns as character vectors (one value per row).
  fbb <- as.character(us_within_threshold$future_bb_id) # future trunk id for each upstream pair
  
  # Extract the two trunk-id columns as character vectors (one value per row).
  cbb <- as.character(us_within_threshold$current_bb_id) # current trunk id for each upstream pair
  
  # Convert trunk ids into row indices on the trunk distance matrix.
  ri <- match(fbb, rownames(trunk_dist_mat)) # row index per pair (future trunk)
  
  # Convert trunk ids into column indices on the trunk distance matrix.
  ci <- match(cbb, colnames(trunk_dist_mat)) # column index per pair (current trunk)
  
  # Identify which rows have a valid (row, col) lookup.
  ok <- !is.na(ri) & !is.na(ci) # TRUE where both indices exist
  
  # Initialize the output with NA (same length as the upstream table).
  trunk_step_dist <- rep(NA_integer_, length(fbb)) # hop distance (trunk steps)
  
  # Fill hop distances for valid pairs using one vectorized matrix lookup.
  trunk_step_dist[ok] <- as.integer(trunk_dist_mat[cbind(ri[ok], ci[ok])]) # hop count per valid row
  
  # Add the trunk-step distance column back onto the upstream candidate table.
  us_with_steps <- us_within_threshold %>%
    dplyr::mutate(trunk_step_dist = trunk_step_dist) # attach hop-distance column
  
  # =============================================================================
  # SECTION 9 — choose the closest upstream dam and compute cascade_level
  # =============================================================================
  
  # For each future dam, choose the closest upstream current dam within the threshold.
  chosen_upstream <- us_with_steps %>%
    dplyr::group_by(.data$future_dam_id) %>% # group by future dam
    dplyr::arrange(.data$dist_km) %>% # choose closest by distance
    dplyr::slice(1L) %>% # keep the closest upstream row
    dplyr::ungroup() %>% # drop grouping
    dplyr::rename(
      chosen_us_current_dam_id = .data$current_dam_id, # chosen upstream current dam
      chosen_us_dist_km = .data$dist_km, # chosen upstream distance
      chosen_us_trunk_steps = .data$trunk_step_dist # chosen upstream trunk-step distance
    ) %>%
    dplyr::select(.data$future_dam_id, .data$chosen_us_current_dam_id, .data$chosen_us_dist_km, .data$chosen_us_trunk_steps) # keep key columns
  
  # Build the future output table (one row per future dam).
  reach_future <- data.frame(
    dam_id = future_dam_ids, # future dam id
    dam_type = rep("future", n_fut), # label future dams
    stringsAsFactors = FALSE # keep strings
  ) %>%
    dplyr::left_join(dplyr::select(future_trunks, future_dam_id, bb_id), by = c("dam_id" = "future_dam_id")) %>% # add bb_id
    dplyr::left_join(downstream_summary, by = c("dam_id" = "future_dam_id")) %>% # add downstream flags/distances
    dplyr::mutate(has_current_downstream = ifelse(is.na(.data$has_current_downstream), FALSE, .data$has_current_downstream)) %>% # NA -> FALSE
    dplyr::left_join(chosen_upstream, by = c("dam_id" = "future_dam_id")) %>% # add chosen upstream row
    dplyr::mutate( # fill upstream outputs + cascade level
      has_current_upstream = !is.na(.data$chosen_us_current_dam_id), # upstream exists if we chose one
      min_distance_upstream_km = .data$chosen_us_dist_km, # chosen upstream distance is the min within threshold
      cascade_level = dplyr::if_else( # apply strict downstream requirement first
        .data$has_current_downstream & .data$has_current_upstream, # both must exist
        as.integer(.data$chosen_us_trunk_steps + 1L), # cascade level = trunk steps + 1
        0L # otherwise not a cascade
      )
    ) %>%
    dplyr::select( # keep only the standard reach_df columns
      .data$dam_id,
      .data$dam_type,
      bb_id = .data$bb_id,
      .data$has_current_downstream,
      .data$min_distance_downstream_km,
      .data$has_current_upstream,
      .data$min_distance_upstream_km,
      .data$cascade_level
    )
  
  # =============================================================================
  # SECTION 10 — append current dams (kept for a complete dam list)
  # =============================================================================
  
  # Build the current output rows (same columns, but cascade fields are NA).
  reach_current <- data.frame(
    dam_id = current_dam_ids, # current dam id
    dam_type = rep("current", n_cur), # label current dams
    bb_id = current_trunks$bb_id, # trunk id
    has_current_downstream = rep(NA, n_cur), # not computed for current dams
    min_distance_downstream_km = rep(NA_real_, n_cur), # not computed for current dams
    has_current_upstream = rep(NA, n_cur), # not computed for current dams
    min_distance_upstream_km = rep(NA_real_, n_cur), # not computed for current dams
    cascade_level = rep(NA_integer_, n_cur), # NA for current dams
    stringsAsFactors = FALSE # keep strings
  )
  
  # Combine future + current rows.
  reach_df <- dplyr::bind_rows(reach_future, reach_current) # final output
  
  # =============================================================================
  # SECTION 11 — decision table (one row per future dam)
  # =============================================================================
  
  # Build a one-row-per-future table that shows the decisions used to compute cascade_level.
  decision_table <- data.frame(
    future_dam_id = future_dam_ids, # future dam id
    stringsAsFactors = FALSE # keep strings
  ) %>%
    dplyr::left_join(dplyr::rename(future_trunks, future_dam_id = .data$future_dam_id), by = "future_dam_id") %>% # add future bb_id
    dplyr::left_join(downstream_summary, by = c("future_dam_id" = "future_dam_id")) %>% # add downstream summary
    dplyr::left_join(chosen_downstream, by = c("future_dam_id" = "future_dam_id")) %>% # add chosen downstream
    dplyr::left_join(chosen_upstream, by = c("future_dam_id" = "future_dam_id")) %>% # add chosen upstream
    dplyr::left_join(dplyr::select(reach_future, dam_id, cascade_level), by = c("future_dam_id" = "dam_id")) # add final cascade_level
  
  # =============================================================================
  # SECTION 12 — debug output: save EVERY step as a table
  # =============================================================================
  
  debug <- list(
    nodes_tbl = nodes_tbl, # all nodes (including non-dam nodes)
    edges_tbl = edges_tbl, # all edges
    weights_km = weights_km, # weight vector used in distances()
    dam_nodes_tbl = dam_nodes_tbl, # all dam nodes
    current_dams_tbl = current_dams_tbl, # current dam nodes only
    future_dams_tbl = future_dams_tbl, # future dam nodes only
    node_bb = node_bb, # node -> bb_id mapping
    ds_mat = ds_mat, # downstream distance matrix (future -> current)
    us_mat = us_mat, # upstream distance matrix (current -> future)
    ds_grid = ds_grid, # downstream index grid (fut_idx/cur_idx + dist_km)
    us_grid = us_grid, # upstream index grid (cur_idx/fut_idx + dist_km)
    ds_dam_all = ds_dam_all, # all downstream dam pairs with ids + bb_id
    us_dam_all = us_dam_all, # all upstream dam pairs with ids + bb_id
    ds_reachable = ds_reachable, # downstream pairs where a path exists
    ds_same_trunk = ds_same_trunk, # strict same-trunk downstream pairs
    chosen_downstream = chosen_downstream, # chosen downstream dam per future dam
    downstream_summary = downstream_summary, # per-future downstream flags/distances
    trunk_edges = trunk_edges, # trunk adjacency edges
    g_trunk = g_trunk, # trunk adjacency graph
    trunk_dist_mat = trunk_dist_mat, # trunk-step distance matrix
    us_reachable = us_reachable, # upstream pairs where a path exists
    us_within_threshold = us_within_threshold, # upstream pairs within km threshold
    us_with_steps = us_with_steps, # upstream pairs with trunk-step distance
    chosen_upstream = chosen_upstream, # chosen upstream dam per future dam
    decision_table = decision_table, # per-future decision summary table
    future_trunks = future_trunks, # future dam -> bb_id
    current_trunks = current_trunks # current dam -> bb_id
  )
  
  # Return final output and debug tables.
  list(reach_df = reach_df, debug = debug)
}


# -----------------------------------------------------------------------------
# empty_reach_connectivity_df
#
# Zero rows but correct column names and types, so callers always get a familiar
# table shape when there are no blended dams on the network.
# -----------------------------------------------------------------------------
empty_reach_connectivity_df <- function() {
  data.frame(
    dam_id = character(),
    dam_type = character(),
    has_current_downstream = logical(),
    has_current_upstream = logical(),
    min_distance_downstream_km = numeric(),
    min_distance_upstream_km = numeric(),
    cascade_level = integer(),
    bb_id = character(),
    stringsAsFactors = FALSE
  )
}
