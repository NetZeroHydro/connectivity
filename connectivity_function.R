# =============================================================================
# connectivity_function.R
# =============================================================================
# Build connectivity matrices and reach-level dataframe from a river network
# with dams already snapped (e.g. output of st_network_blend). Use this for
# a single basin or sub-basin without re-running the full connectivity_df
# pipeline.
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
                                      cascade_threshold_km = NULL,
                                      cascade_trunk_steps = 0L) {
  
  # ---------------------------------------------------------------------------
  # connectivity_from_network — what this function does (walkthrough)
  #
  # Directed river graph (downstream along edges). For each future dam, among
  # current dams on the same trunk (or nearby trunks if steps > 0), find whether
  # some current dam is upstream and some downstream in river km, then label
  # cascade (optionally requiring both distances under a km cap). Append current
  # dam rows with NA in those cascade columns.
  # ---------------------------------------------------------------------------
  
  # How far we walk in the “trunk graph” when steps > 0 (must be whole, non-negative)
  steps <- suppressWarnings(as.integer(cascade_trunk_steps))
  if (length(steps) != 1L || is.na(steps) || steps < 0L) {
    stop("cascade_trunk_steps must be a non-negative integer.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Nodes as a plain table
  #
  # We need dam ids and which nodes are current vs future. node_id is just the
  # row number in the node table — that matches how igraph indexes vertices in
  # this sfnetwork. If the network was never blended, those columns are missing
  # and we return an empty result with the right column names.
  # ---------------------------------------------------------------------------
  nodes <- net_with_dams %>%
    tidygraph::activate("nodes") %>%
    sf::st_as_sf() %>%
    sf::st_drop_geometry() %>% # avoid sf quirks in joins later
    dplyr::as_tibble() %>%
    dplyr::mutate(node_id = dplyr::row_number())
  
  if (!"dam_id" %in% names(nodes) || !"is_current_dam" %in% names(nodes)) {
    return(empty_reach_connectivity_df())
  }
  
  # Node indices (not dam_id strings yet) for igraph distance calls
  current_dam_nodes <- nodes %>%
    dplyr::filter(.data$is_current_dam %in% TRUE) %>%
    dplyr::pull(.data$node_id)
  
  future_dam_nodes <- nodes %>%
    dplyr::filter(.data$is_current_dam %in% FALSE) %>%
    dplyr::pull(.data$node_id)
  
  n_cur <- length(current_dam_nodes)
  n_fut <- length(future_dam_nodes)
  
  # ---------------------------------------------------------------------------
  # Edge lengths for shortest paths
  #
  # igraph wants a numeric weight per edge in km here. We divide metres by 1000.
  # If weight is missing or wrong type, we stop with a hint about projection.
  # ---------------------------------------------------------------------------
  edges_tbl <- net_with_dams %>%
    tidygraph::activate("edges") %>%
    sf::st_as_sf() %>%
    sf::st_drop_geometry() %>%
    dplyr::as_tibble()
  
  weights_km <- tryCatch(
    {
      as.numeric(edges_tbl$weight) / 1000
    },
    error = function(e) {
      stop(
        "Edge `weight` could not be converted to kilometres. ",
        "Make sure `weight` represents length in metres (e.g. edge_length()/st_length()) ",
        "before calling connectivity_from_network(). If your rivers are lon/lat, ",
        "either enable s2 geodesic lengths (sf::sf_use_s2(TRUE)) or project to a ",
        "metric CRS (st_transform) before computing `weight`.",
        call. = FALSE
      )
    }
  )
  
  # ---------------------------------------------------------------------------
  # All-pairs distances used later per future dam
  #
  # Downstream matrix: rows = current dam nodes, cols = future dam nodes, values =
  # km along directed edges from that current to that future. Upstream matrix:
  # from future to current (still mode "out", but start at future). If either
  # side has no dams, we skip igraph (empty v/to breaks) and use placeholder
  # matrices so the rest of the code can still run.
  # ---------------------------------------------------------------------------
  if (n_cur > 0L && n_fut > 0L) {
    connectivity_matrix_downstream <- igraph::distances(
      net_with_dams,
      v = current_dam_nodes,
      to = future_dam_nodes,
      weights = weights_km,
      mode = "out"
    )
    
    connectivity_matrix_upstream <- igraph::distances(
      net_with_dams,
      v = future_dam_nodes,
      to = current_dam_nodes,
      weights = weights_km,
      mode = "out"
    )
  } else {
    connectivity_matrix_downstream <- matrix(NA_real_, nrow = n_cur, ncol = n_fut)
    connectivity_matrix_upstream <- matrix(NA_real_, nrow = n_fut, ncol = n_cur)
  }
  
  # ---------------------------------------------------------------------------
  # Which river trunk (bb_id) is each dam on?
  #
  # Each edge carries bb_id. At each node, look at all touching edges and pick
  # one bb_id as that node’s trunk (first in group — good enough when dams sit
  # on a single reach). We need this to restrict which current dams count for
  # each future dam when steps = 0 (same trunk only).
  # ---------------------------------------------------------------------------
  if (!"bb_id" %in% names(edges_tbl)) {
    stop(
      "Edges must have attribute `bb_id` (backbone / trunk ID).",
      call. = FALSE
    )
  }
  
  # Endpoints of edges: both “from” and “to” get a row with that edge’s bb_id
  incident_from <- edges_tbl %>%
    dplyr::select(dplyr::all_of(c("from", "bb_id"))) %>%
    dplyr::rename(node_id = .data$from)
  incident_to <- edges_tbl %>%
    dplyr::select(dplyr::all_of(c("to", "bb_id"))) %>%
    dplyr::rename(node_id = .data$to)
  
  incident <- dplyr::bind_rows(incident_from, incident_to) %>%
    dplyr::filter(!is.na(.data$bb_id)) %>%
    dplyr::group_by(.data$node_id) %>%
    dplyr::summarise(trunk_id = dplyr::first(.data$bb_id), .groups = "drop")
  
  nodes_trunk <- nodes %>%
    dplyr::left_join(incident, by = "node_id")
  
  current_trunk <- if (n_cur > 0L) {
    nodes_trunk %>%
      dplyr::slice(current_dam_nodes) %>%
      dplyr::pull(.data$trunk_id)
  } else {
    character(0)
  }
  
  future_trunk <- if (n_fut > 0L) {
    nodes_trunk %>%
      dplyr::slice(future_dam_nodes) %>%
      dplyr::pull(.data$trunk_id)
  } else {
    character(0)
  }
  
  future_dam_ids <- if (n_fut > 0L) {
    nodes %>%
      dplyr::slice(future_dam_nodes) %>%
      dplyr::pull(.data$dam_id)
  } else {
    character(0)
  }
  
  # ---------------------------------------------------------------------------
  # Optional: treat nearby trunks as “neighbors”
  #
  # At a confluence, several bb_id values meet. We build an undirected graph
  # whose vertices are trunk ids and edges connect trunks that share a node.
  # Then “allowed” trunks for a future dam are those within `steps` hops on that
  # graph. If steps = 0, we only allow the future dam’s own trunk id.
  # ---------------------------------------------------------------------------
  bb_from_edges <- as.character(stats::na.omit(unique(as.character(edges_tbl$bb_id))))
  bb_from_dams <- as.character(stats::na.omit(unique(c(
    as.character(current_trunk),
    as.character(future_trunk)
  ))))
  bb_vertices <- unique(c(bb_from_edges, bb_from_dams))
  
  if (steps == 0L) {
    allowed_trunks_for <- function(tf) {
      tf_chr <- as.character(tf)
      if (length(tf_chr) != 1L || is.na(tf_chr) || !nzchar(tf_chr)) {
        return(character(0))
      }
      tf_chr
    }
  } else {
    node_ids <- unique(c(edges_tbl$from, edges_tbl$to))
    
    trunk_pairs_from <- character()
    trunk_pairs_to <- character()
    
    for (n in node_ids) {
      inc <- edges_tbl[edges_tbl$from == n | edges_tbl$to == n, , drop = FALSE]
      bb <- unique(as.character(inc$bb_id))
      bb <- bb[!is.na(bb) & nzchar(bb)]
      if (length(bb) >= 2L) {
        cm <- utils::combn(bb, 2L) # all pairs of trunks at this junction
        trunk_pairs_from <- c(trunk_pairs_from, as.character(cm[1L, ]))
        trunk_pairs_to <- c(trunk_pairs_to, as.character(cm[2L, ]))
      }
    }
    
    trunk_edges_df <- data.frame(
      from = trunk_pairs_from,
      to = trunk_pairs_to,
      stringsAsFactors = FALSE
    )
    if (nrow(trunk_edges_df) > 0L) trunk_edges_df <- unique(trunk_edges_df)
    
    g_trunk <- igraph::graph_from_data_frame(
      trunk_edges_df,
      directed = FALSE,
      vertices = data.frame(name = bb_vertices, stringsAsFactors = FALSE)
    )
    trunk_dist <- igraph::distances(g_trunk, mode = "all")
    vnames <- igraph::V(g_trunk)$name
    dimnames(trunk_dist) <- list(vnames, vnames)
    
    allowed_trunks_for <- function(tf) {
      tf_chr <- as.character(tf)
      if (length(tf_chr) != 1L || is.na(tf_chr) || !nzchar(tf_chr)) {
        return(character(0))
      }
      if (!tf_chr %in% rownames(trunk_dist)) {
        return(tf_chr)
      }
      di <- trunk_dist[tf_chr, , drop = TRUE]
      names(di)[is.finite(di) & di <= steps]
    }
  }
  
  # ---------------------------------------------------------------------------
  # Fill in connectivity columns for each future dam
  #
  # For future dam j, only look at current dams whose trunk is in allowed_bb.
  # Downstream distance: positive finite entry in downstream matrix means there
  # is a path from current to future along the river; we want the smallest such
  # distance. Upstream: from future to current, same idea. Zero is same node —
  # we ignore for “nearest along river” by requiring > 0.
  # ---------------------------------------------------------------------------
  reach_future <- data.frame(
    dam_id = as.character(future_dam_ids),
    dam_type = rep("future", n_fut),
    stringsAsFactors = FALSE
  )
  
  has_current_downstream_tr <- logical(n_fut)
  has_current_upstream_tr <- logical(n_fut)
  min_distance_downstream_tr <- rep(NA_real_, n_fut)
  min_distance_upstream_tr <- rep(NA_real_, n_fut)
  
  current_trunk_chr <- as.character(current_trunk)
  
  for (j in seq_len(n_fut)) {
    allowed_bb <- allowed_trunks_for(future_trunk[j])
    if (length(allowed_bb) == 0L) next
    
    same_idx <- which(!is.na(current_trunk) & current_trunk_chr %in% allowed_bb)
    if (length(same_idx) == 0L) next
    
    if (n_cur > 0L) {
      col_down <- connectivity_matrix_downstream[same_idx, j, drop = FALSE]
      row_up <- connectivity_matrix_upstream[j, same_idx, drop = FALSE]
      
      d_down <- col_down[is.finite(col_down) & col_down > 0]
      d_up <- row_up[is.finite(row_up) & row_up > 0]
      
      if (length(d_down) > 0L) {
        has_current_downstream_tr[j] <- TRUE
        min_distance_downstream_tr[j] <- min(d_down)
      }
      if (length(d_up) > 0L) {
        has_current_upstream_tr[j] <- TRUE
        min_distance_upstream_tr[j] <- min(d_up)
      }
    }
  }
  
  reach_future$has_current_downstream <- has_current_downstream_tr
  reach_future$has_current_upstream <- has_current_upstream_tr
  reach_future$min_distance_downstream_km <- min_distance_downstream_tr
  reach_future$min_distance_upstream_km <- min_distance_upstream_tr
  reach_future$bb_id <- as.character(future_trunk)
  
  # ---------------------------------------------------------------------------
  # Cascade label
  #
  # Simple rule: need both an upstream and a downstream current dam in the
  # allowed trunk set. If cascade_threshold_km is set, those two shortest
  # distances must also both be under that many km.
  # ---------------------------------------------------------------------------
  threshold_km <- if (is.null(cascade_threshold_km)) NA_real_ else cascade_threshold_km
  
  if (is.na(threshold_km)) {
    reach_future$cascade_status <- ifelse(
      reach_future$has_current_upstream & reach_future$has_current_downstream,
      "cascade",
      "not in cascade"
    )
  } else {
    reach_future$cascade_status <- ifelse(
      reach_future$has_current_upstream & reach_future$has_current_downstream &
        !is.na(reach_future$min_distance_upstream_km) &
        !is.na(reach_future$min_distance_downstream_km) &
        reach_future$min_distance_upstream_km <= threshold_km &
        reach_future$min_distance_downstream_km <= threshold_km,
      "cascade",
      "not in cascade"
    )
  }
  
  reach_future <- reach_future %>%
    dplyr::arrange(
      .data$cascade_status,
      .data$min_distance_upstream_km,
      .data$min_distance_downstream_km
    )
  
  # ---------------------------------------------------------------------------
  # Current dams: keep them in the output for a complete dam list
  #
  # We do not recompute cascade for them; distances and cascade columns stay NA.
  # bb_id still shows which trunk they sit on (useful for maps and joins).
  # ---------------------------------------------------------------------------
  reach_current <- NULL
  if (n_cur > 0L) {
    current_dam_ids <- nodes %>%
      dplyr::slice(current_dam_nodes) %>%
      dplyr::pull(.data$dam_id)
    
    reach_current <- data.frame(
      dam_id = as.character(current_dam_ids),
      dam_type = rep("current", n_cur),
      has_current_downstream = rep(NA, n_cur),
      has_current_upstream = rep(NA, n_cur),
      min_distance_downstream_km = rep(NA_real_, n_cur),
      min_distance_upstream_km = rep(NA_real_, n_cur),
      cascade_status = rep(NA_character_, n_cur),
      bb_id = as.character(current_trunk),
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(reach_current)) {
    reach_future
  } else {
    dplyr::bind_rows(reach_future, reach_current)
  }
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
    cascade_status = character(),
    bb_id = character(),
    stringsAsFactors = FALSE
  )
}
