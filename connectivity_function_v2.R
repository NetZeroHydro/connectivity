# =============================================================================
# connectivity_function.R
# =============================================================================
# Connectivity for a blended dam network: builds `reach_df` (one row per dam
# node) plus a `debug` list of tibbles.
#
# Prerequisites: directed `sfnetwork` with edge `weight` (metres) and `bb_id`;
# blended nodes carry `dam_id` and `is_current_dam` (TRUE = current, FALSE =
# future). Load: sfnetworks, igraph, dplyr, sf, tidygraph before sourcing.
#
# - Downstream “next” current dam: closest same-trunk downstream (river km).
# - Upstream “next” current dam: **tier by trunk hops first** (0 = same trunk,
#   then 1 hop on the confluence trunk graph, then 2, …); **within** a tier,
#   pick smallest river km. No km cap on connectivity. Rows with unknown hop
#   count (NA) sort after all finite tiers.
# - Cascade upstream: same tier rule, but only pairs with
#   `dist_km <= cascade_threshold_km`; then same-trunk downstream must exist for
#   `cascade_level > 0`.
# =============================================================================

#' Connectivity for a snapped dam network (`reach_df` + `debug`)
#'
#' Builds directed river distances between **current** and **future** dam nodes
#' on a blended `sfnetwork`, applies cascade rules, and returns a **list** with
#' `reach_df` (one row per dam) and `debug` (intermediate tables for inspection).
#'
#' **Dam nodes:** nodes with non-missing `dam_id`; `is_current_dam` is `TRUE` for
#' current dams and `FALSE` for future dams.
#'
#' **Edge weights:** column `weight` is treated as **metres** and converted to
#' **km** for [igraph::distances()] (`as.numeric(weight) / 1000`). Edges must
#' carry trunk id `bb_id` (used for same-trunk checks and the trunk graph).
#'
#' **Downstream (strict same trunk):** nearest downstream **current** dam on the
#' same `bb_id` as the future dam (finite path, `dist_km > 0`). Exposed as
#' `min_distance_downstream_km` and `dam_id_down`.
#'
#' **Upstream (no km threshold; tiered by trunk):** among finite current-to-future
#' paths with `dist_km > 0`, choose the upstream current with **smallest trunk-graph
#' hop distance** from the future dam’s trunk (0 = same trunk, then 1, 2, …). If
#' several ties on hop count, pick the **smallest river km** in that tier. That
#' choice defines `has_current_upstream`, `min_distance_upstream_km`, `dam_id_up`,
#' and `us_trunks_away`. Pairs with unknown hop count (`NA`) are considered only
#' after all finite tiers.
#'
#' **Cascade (`cascade_threshold_km` only):** same tier rule, but only pairs with
#' `dist_km <= cascade_threshold_km`. If there is also a qualifying same-trunk
#' downstream current dam, `cascade_level = trunk_hops + 1` for that cascade
#' upstream pair; otherwise `0`. Current-dam rows use `NA` for connectivity/cascade
#' fields.
#'
#' @param net_with_dams Directed `sfnetwork` with edge `bb_id`, numeric edge
#'   `weight` (metres), and blended node columns `dam_id` / `is_current_dam` when
#'   dams were snapped on.
#' @param cascade_threshold_km Maximum river-network distance (km) for the
#'   **cascade** upstream candidate only (not for `min_distance_upstream_km`).
#'   Default `100`.
#'
#' @return A named list:
#'
#'  {`reach_df`} Tibble/data frame with `dam_id`, `dam_type`, `bb_id`,
#'     downstream/upstream flags, min distances (km), `dam_id_down`, `dam_id_up`,
#'     `us_trunks_away`, and integer `cascade_level`.
#'   {`debug`} List of intermediate objects (pairwise tables, trunk graph,
#'     nearest vs cascade upstream, `decision_table`, etc.).
#' }
#'
#' @details
#' Distances use `igraph::distances(..., mode = "out")` on the directed graph.
#' `Inf` means no path; `0` means same vertex (pairs with distance 0 are not used
#' as reachable in the filters below).
#'
#' @examples
#' \dontrun{
#' out <- connectivity_from_network(out$net_with_dams, cascade_threshold_km = 100)
#' reach_df <- out$reach_df
#' dbg <- out$debug
#' }
#'
#' @md
#' @export
connectivity_from_network <- function(net_with_dams,
                                      cascade_threshold_km = 100) {
  
  # =============================================================================
  # SECTION 0 — what this function returns
  #
  #   - reach_df: one row per dam node (future rows filled in; current rows
  #     carry NA in connectivity/cascade columns because those notions are
  #     defined from the future dam’s perspective).
  #   - debug: named data frames so you can inspect each filtering / join step
  #     when debugging a basin (print head(), join to maps, etc.).
  # =============================================================================
  
  # =============================================================================
  # SECTION 1 — Strip nodes and edges to plain tables
  #
  # sfnetwork stores geometry on nodes (points) and edges (lines). For table logic
  # we drop geometry and work in dplyr: that keeps joins fast and avoids
  # accidentally carrying big sf columns through every mutate.
  # We add node_id = row number because that matches igraph vertex indices when
  # the network was built with vertices in row order (typical tidygraph path).
  # =============================================================================
  
  # Activate the node table inside the sfnetwork object, then flatten to tibble.
  nodes_tbl <- net_with_dams %>%
    tidygraph::activate("nodes") %>% # switch context from edges to nodes
    sf::st_as_sf() %>% # need sf first so st_drop_geometry() is defined
    sf::st_drop_geometry() %>% # drop POINT geometry; keep attribute columns only
    dplyr::as_tibble() %>% # standard dplyr table
    dplyr::mutate(node_id = dplyr::row_number()) # 1..N == igraph vertex id
  
  # Same pattern for edges: we need from, to, weight, bb_id as plain columns.
  edges_tbl <- net_with_dams %>%
    tidygraph::activate("edges") %>% # edge table holds line geometry + attrs
    sf::st_as_sf() %>%
    sf::st_drop_geometry() %>% # drop LINESTRING; lengths already in `weight`
    dplyr::as_tibble()
  
  # igraph distances() expects one positive weight per edge; we use km not m.
  weights_km <- as.numeric(edges_tbl$weight) / 1000 # coerce then metres → km
  
  # =============================================================================
  # SECTION 2 — identify dam nodes (current vs future)
  #
  # After st_network_blend(), dam sites become graph nodes with dam_id set.
  # Non-dam river junctions have NA dam_id — we ignore them for dam–dam logic.
  # is_current_dam distinguishes operational (TRUE) vs planned (FALSE) dams.
  # =============================================================================
  
  # Keep only vertices that represent a dam; coerce dam_id to character for joins.
  dam_nodes_tbl <- nodes_tbl %>%
    dplyr::filter(!is.na(.data$dam_id)) %>% # NA dam_id => not a dam vertex
    dplyr::mutate(dam_id = as.character(.data$dam_id)) # avoid factor surprises
  
  # Current dams: the ones already built / in the “current” inventory.
  current_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% TRUE) %>% # strict TRUE (not NA)
    dplyr::select(node_id, dam_id, is_current_dam) # minimal columns for indices
  
  # Future dams: planned sites (FALSE); same column set for symmetry.
  future_dams_tbl <- dam_nodes_tbl %>%
    dplyr::filter(.data$is_current_dam %in% FALSE) %>%
    dplyr::select(node_id, dam_id, is_current_dam)
  
  # Vertex indices passed to igraph::distances (must align with node_id order).
  current_nodes <- current_dams_tbl$node_id # length n_cur
  future_nodes <- future_dams_tbl$node_id # length n_fut
  
  # Readable dam identifiers (parallel to the node index vectors above).
  current_dam_ids <- current_dams_tbl$dam_id
  future_dam_ids <- future_dams_tbl$dam_id
  
  # Counts drive expand.grid sizes and matrix dimensions.
  n_cur <- length(current_nodes) # number of current dam vertices
  n_fut <- length(future_nodes) # number of future dam vertices
  
  # =============================================================================
  # SECTION 3 — distance matrices (km) for every future–current pair
  #
  # mode = "out" follows directed edges downstream. So:
  #   - ds_mat[i,j]: start at future i, can we reach current j downstream? how far?
  #   - us_mat[j,i]: start at current j, can we reach future i downstream? if yes,
  #     current j lies UPSTREAM of future i along the network (water flows from
  #     current toward future in the directed sense).
  # We exclude self-pairs with dist 0 later — those are same-node, not “along river”.
  # =============================================================================
  
  # Future → current: “is this current dam downstream of this future dam?”
  ds_mat <- igraph::distances(
    net_with_dams, # whole directed river graph
    v = future_nodes, # rows: each future dam vertex
    to = current_nodes, # cols: each current dam vertex
    weights = weights_km, # edge lengths in km (Inf if no edge weight)
    mode = "out" # follow arrows downstream from future
  )
  
  # Current → future: positive entry means a directed path current → future.
  us_mat <- igraph::distances(
    net_with_dams,
    v = current_nodes, # start at each current dam
    to = future_nodes, # end at each future dam
    weights = weights_km,
    mode = "out" # current upstream of future if such a path exists
  )
  
  # =============================================================================
  # SECTION 4 — build long tables ds_dam_all and us_dam_all (ALL pairs)
  #
  # Matrices are hard to read and join; we expand to one row per pair.
  # CRITICAL: as.vector(ds_mat) goes column-wise; expand.grid must use the same
  # ordering (fut_idx runs fast, cur_idx slow) so row i of the grid matches
  # element i of the vectorized matrix.
  # =============================================================================
  
  # All combinations of (future index, current index) for downstream matrix.
  ds_grid <- expand.grid(
    fut_idx = seq_len(n_fut), # 1 .. n_fut
    cur_idx = seq_len(n_cur) # 1 .. n_cur
  )
  # Flatten ds_mat in column-major order to align with expand.grid rows.
  ds_grid$dist_km <- as.vector(ds_mat) # future→current km (may be Inf/0)
  
  # Attach real node ids and dam ids so the table is self-describing.
  ds_dam_all <- ds_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx], # map index → vertex id
      current_node_id = current_nodes[.data$cur_idx],
      future_dam_id = future_dam_ids[.data$fut_idx],
      current_dam_id = current_dam_ids[.data$cur_idx]
    ) %>%
    dplyr::select(
      .data$future_dam_id,
      .data$current_dam_id,
      .data$future_node_id,
      .data$current_node_id,
      dist_km = .data$dist_km # keep one distance column with clear name
    )
  
  # Upstream matrix rows = current, cols = future; grid order matches as.vector(us_mat).
  us_grid <- expand.grid(
    cur_idx = seq_len(n_cur), # must match us_mat row order
    fut_idx = seq_len(n_fut) # must match us_mat col order
  )
  us_grid$dist_km <- as.vector(us_mat) # current→future km
  
  us_dam_all <- us_grid %>%
    dplyr::mutate(
      future_node_id = future_nodes[.data$fut_idx],
      current_node_id = current_nodes[.data$cur_idx],
      future_dam_id = future_dam_ids[.data$fut_idx],
      current_dam_id = current_dam_ids[.data$cur_idx]
    ) %>%
    dplyr::select(
      .data$future_dam_id,
      .data$current_dam_id,
      .data$future_node_id,
      .data$current_node_id,
      dist_km = .data$dist_km # river km along directed path current → future
    )
  
  # =============================================================================
  # SECTION 5 — node → bb_id (“which trunk is this vertex on?”)
  #
  # Each edge carries bb_id (hydrobasin trunk segment id). A node can touch
  # multiple edges; we take dplyr::first() after group_by as a deterministic
  # representative trunk for that vertex (good enough for dam vertices on one reach).
  # We then join future_bb_id and current_bb_id onto every pair row for same-trunk tests.
  # =============================================================================
  
  # Every directed edge has a “from” and “to” vertex; both need a trunk label.
  inc_from <- edges_tbl %>% dplyr::select(node_id = .data$from, bb_id = .data$bb_id)
  inc_to <- edges_tbl %>% dplyr::select(node_id = .data$to, bb_id = .data$bb_id)
  
  # Stack so each node appears once per incident edge (may duplicate node_id).
  inc_all <- dplyr::bind_rows(inc_from, inc_to)
  
  # One bb_id per node_id: first non-group column order from incident edges.
  node_bb <- inc_all %>%
    dplyr::filter(!is.na(.data$bb_id)) %>% # cannot assign trunk without bb_id
    dplyr::group_by(.data$node_id) %>%
    dplyr::summarise(bb_id = dplyr::first(.data$bb_id), .groups = "drop")
  
  # Downstream pairs: need both endpoints’ trunks to enforce same-trunk rule.
  ds_dam_all <- ds_dam_all %>%
    dplyr::left_join(node_bb, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = .data$bb_id) %>% # trunk at future dam vertex
    dplyr::left_join(node_bb, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = .data$bb_id) # trunk at current dam vertex
  
  # Upstream pairs: same joins — needed for trunk hop distance between trunks.
  us_dam_all <- us_dam_all %>%
    dplyr::left_join(node_bb, by = c("future_node_id" = "node_id")) %>%
    dplyr::rename(future_bb_id = .data$bb_id) %>%
    dplyr::left_join(node_bb, by = c("current_node_id" = "node_id")) %>%
    dplyr::rename(current_bb_id = .data$bb_id)
  
  # Small lookup tables for output bb_id on reach_df rows.
  future_trunks <- future_dams_tbl %>%
    dplyr::left_join(node_bb, by = c("node_id" = "node_id")) %>%
    dplyr::rename(future_node_id = .data$node_id, future_dam_id = .data$dam_id, bb_id = .data$bb_id) %>%
    dplyr::select(.data$future_dam_id, .data$future_node_id, .data$bb_id)
  
  current_trunks <- current_dams_tbl %>%
    dplyr::left_join(node_bb, by = c("node_id" = "node_id")) %>%
    dplyr::rename(current_node_id = .data$node_id, current_dam_id = .data$dam_id, bb_id = .data$bb_id) %>%
    dplyr::select(.data$current_dam_id, .data$current_node_id, .data$bb_id)
  
  # =============================================================================
  # SECTION 6 — downstream: strict same trunk, nearest current = dam_id_down
  #
  # Policy: “downstream current neighbor” must sit on the SAME bb_id as the future
  # dam. That avoids counting currents on a parallel branch that met at a confluence.
  # We still require dist_km > 0 so we never treat “same graph vertex” as a link.
  # =============================================================================
  
  # Any finite strictly positive path future → current counts as “downstream reachable”.
  ds_reachable <- ds_dam_all %>%
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0) # drop Inf and 0
  
  # Subset to pairs where both ends share one trunk id (string compare).
  ds_same_trunk <- ds_reachable %>%
    dplyr::filter(.data$current_bb_id == .data$future_bb_id)
  
  # Per future dam: pick the CLOSEST same-trunk downstream current (km along network).
  chosen_downstream <- ds_same_trunk %>%
    dplyr::group_by(.data$future_dam_id) %>% # one winner per future
    dplyr::arrange(.data$dist_km) %>% # smallest river distance first
    dplyr::slice(1L) %>% # take the top row after sort
    dplyr::ungroup() %>%
    dplyr::rename(
      chosen_ds_current_dam_id = .data$current_dam_id, # becomes dam_id_down later
      chosen_ds_dist_km = .data$dist_km
    ) %>%
    dplyr::select(.data$future_dam_id, .data$chosen_ds_current_dam_id, .data$chosen_ds_dist_km)
  
  # Aggregate: boolean “any downstream?” + min km (NA if none on same trunk).
  downstream_summary <- ds_same_trunk %>%
    dplyr::group_by(.data$future_dam_id) %>%
    dplyr::summarise(
      has_current_downstream = dplyr::n() > 0, # TRUE if ds_same_trunk non-empty
      min_distance_downstream_km = ifelse(dplyr::n() > 0, min(.data$dist_km), NA_real_),
      .groups = "drop"
    )
  
  # =============================================================================
  # SECTION 7 — trunk graph at confluences (undirected bb_id adjacency)
  #
  # Where two different bb_id values meet at one graph node, we add an undirected
  # “trunk–trunk” edge. That graph’s shortest-path LENGTH (in hops) is how many
  # trunks apart two dams are for labeling (us_trunks_away, cascade trunk_hops).
  # This is NOT the same as river km — it counts confluence steps between segments.
  # =============================================================================
  
  # Every edge endpoint is a candidate confluence.
  junction_nodes <- unique(c(edges_tbl$from, edges_tbl$to))
  
  # Accumulate unique (bb_id, bb_id) pairs seen at shared nodes.
  trunk_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  
  for (nid in junction_nodes) { # scan each graph vertex
    # All trunk labels on edges touching this node.
    bb_here <- inc_all %>% dplyr::filter(.data$node_id == nid) %>% dplyr::pull(.data$bb_id)
    bb_here <- unique(as.character(bb_here)) # dedupe labels at this node
    bb_here <- bb_here[!is.na(bb_here) & nzchar(bb_here)] # drop bad labels
    if (length(bb_here) >= 2) { # confluence: two or more distinct trunks meet
      pairs <- utils::combn(bb_here, 2) # all unordered trunk pairs at this node
      trunk_edges <- rbind(
        trunk_edges,
        data.frame(from = as.character(pairs[1, ]), to = as.character(pairs[2, ]), stringsAsFactors = FALSE)
      )
    }
  }
  
  trunk_edges <- unique(trunk_edges) # drop duplicate trunk–trunk edges from multiple nodes
  g_trunk <- igraph::graph_from_data_frame(trunk_edges, directed = FALSE) # undirected
  
  # All-pairs shortest PATH LENGTHS in the trunk graph (integer hop counts).
  trunk_dist_mat <- igraph::distances(g_trunk, mode = "all")
  # Name rows/cols by bb_id string so we can match() later (igraph stores names on V()).
  if (igraph::vcount(g_trunk) > 0) {
    vnames <- igraph::V(g_trunk)$name # vertex name = bb_id from data frame
    dimnames(trunk_dist_mat) <- list(vnames, vnames)
  }
  
  # =============================================================================
  # SECTION 8 — upstream reachable pairs + trunk hop column (threshold split later)
  #
  # us_reachable: every current that is upstream of a future (finite km, > 0).
  # trunk_step_dist = hop count on the undirected trunk graph (0 = same bb_id).
  #
  # Selection rule (both connectivity and cascade): **trunk tier first**, then
  # river km. So same-trunk upstream candidates always beat cross-trunk ones,
  # even if a cross-trunk dam is closer in km. Within one tier, smallest dist_km wins.
  #
  #   - nearest_upstream: all us_reachable rows, tiered sort, no km cap
  #   - us_within_threshold: rows with dist_km <= cascade_threshold_km (cascade pool)
  #   - cascade_upstream: same tiered sort on that pool → cascade_level input
  # =============================================================================
  
  # Upstream in network sense: directed path current → future, positive length.
  us_reachable <- us_dam_all %>%
    dplyr::filter(is.finite(.data$dist_km) & .data$dist_km > 0)
  
  # Vectorized trunk hop lookup: map bb_id strings to matrix row/col indices.
  fbb_all <- as.character(us_reachable$future_bb_id) # future’s trunk per row
  cbb_all <- as.character(us_reachable$current_bb_id) # current’s trunk per row
  ri_all <- match(fbb_all, rownames(trunk_dist_mat)) # NA if trunk not in graph
  ci_all <- match(cbb_all, colnames(trunk_dist_mat))
  ok_all <- !is.na(ri_all) & !is.na(ci_all) # both ends must exist in trunk graph
  trunk_step_all <- rep(NA_integer_, length(fbb_all)) # default NA if lookup fails
  # Single matrix index per row: cbind(row_i, col_i) picks one cell each.
  trunk_step_all[ok_all] <- as.integer(trunk_dist_mat[cbind(ri_all[ok_all], ci_all[ok_all])])
  
  # Carry hops beside river km for every upstream candidate row.
  us_reachable_with_steps <- us_reachable %>%
    dplyr::mutate(trunk_step_dist = trunk_step_all) # 0 = same trunk, 1+ = hops away
  
  # Cascade-only slice: still full rows, just filtered by river km policy.
  us_within_threshold <- us_reachable_with_steps %>%
    dplyr::filter(.data$dist_km <= cascade_threshold_km) # parameter applies HERE only
  
  # CONNECTIVITY upstream neighbor: fewest trunk hops first, then min river km.
  nearest_upstream <- us_reachable_with_steps %>%
    dplyr::group_by(.data$future_dam_id) %>%
    dplyr::arrange(
      dplyr::coalesce(.data$trunk_step_dist, .Machine$integer.max), # 0,1,2,… then NA last
      .data$dist_km # within same hop tier, closest along the network
    ) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = .data$future_dam_id,
      dam_id_up = .data$current_dam_id,
      min_distance_upstream_km = .data$dist_km,
      us_trunks_away = .data$trunk_step_dist
    )
  
  # CASCADE upstream neighbor: same tier rule, pool = within km threshold only.
  cascade_upstream <- us_within_threshold %>%
    dplyr::group_by(.data$future_dam_id) %>%
    dplyr::arrange(
      dplyr::coalesce(.data$trunk_step_dist, .Machine$integer.max),
      .data$dist_km
    ) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      future_dam_id = .data$future_dam_id,
      cascade_us_dam_id = .data$current_dam_id,
      cascade_us_dist_km = .data$dist_km,
      cascade_us_trunk_steps = .data$trunk_step_dist
    )
  
  # =============================================================================
  # SECTION 9 — assemble reach_future (one row per future dam)
  #
  # Join order: skeleton dam_id + dam_type, bb_id, downstream summary, dam_id_down,
  # nearest_upstream (tiered upstream pick), cascade_upstream (tiered + threshold),
  # then has_current_upstream and cascade_level, then select() for stable columns.
  # =============================================================================
  
  # Skeleton: every future dam id appears exactly once.
  reach_future <- data.frame(
    dam_id = future_dam_ids,
    dam_type = rep("future", n_fut),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(dplyr::select(future_trunks, future_dam_id, bb_id), by = c("dam_id" = "future_dam_id")) %>%
    dplyr::left_join(downstream_summary, by = c("dam_id" = "future_dam_id")) %>%
    # left_join gives NA for has_* when no downstream row; coerce to FALSE for clarity
    dplyr::mutate(has_current_downstream = ifelse(is.na(.data$has_current_downstream), FALSE, .data$has_current_downstream)) %>%
    dplyr::left_join(
      chosen_downstream %>% dplyr::transmute(
        future_dam_id = .data$future_dam_id,
        dam_id_down = .data$chosen_ds_current_dam_id # explicit neighbor id column
      ),
      by = c("dam_id" = "future_dam_id")
    ) %>%
    dplyr::left_join(nearest_upstream, by = c("dam_id" = "future_dam_id")) %>%
    dplyr::left_join(cascade_upstream, by = c("dam_id" = "future_dam_id")) %>%
    dplyr::mutate(
      # If nearest_upstream missing, dam_id_up is NA → no upstream current found.
      has_current_upstream = !is.na(.data$dam_id_up),
      # Cascade needs BOTH same-trunk downstream AND a threshold-qualified upstream.
      cascade_level = dplyr::if_else(
        .data$has_current_downstream & !is.na(.data$cascade_us_dam_id),
        # trunk hops + 1 is your severity index; coalesce handles NA hops as 0
        as.integer(dplyr::coalesce(.data$cascade_us_trunk_steps, 0L) + 1L),
        0L # not a cascade if either leg fails
      )
    ) %>%
    dplyr::select(
      .data$dam_id,
      .data$dam_type,
      bb_id = .data$bb_id,
      .data$has_current_downstream,
      .data$min_distance_downstream_km,
      .data$dam_id_down,
      .data$has_current_upstream,
      .data$min_distance_upstream_km,
      .data$dam_id_up,
      .data$us_trunks_away,
      .data$cascade_level
    )
  
  # =============================================================================
  # SECTION 10 — append current dam rows (placeholders for a complete dam inventory)
  #
  # Current dams are included so reach_df lines up with “all dam_id on network”.
  # We do not compute up/down connectivity relative to them here — all NA/NA_real_.
  # =============================================================================
  
  reach_current <- data.frame(
    dam_id = current_dam_ids,
    dam_type = rep("current", n_cur),
    bb_id = current_trunks$bb_id,
    has_current_downstream = rep(NA, n_cur),
    min_distance_downstream_km = rep(NA_real_, n_cur),
    dam_id_down = rep(NA_character_, n_cur),
    has_current_upstream = rep(NA, n_cur),
    min_distance_upstream_km = rep(NA_real_, n_cur),
    dam_id_up = rep(NA_character_, n_cur),
    us_trunks_away = rep(NA_integer_, n_cur),
    cascade_level = rep(NA_integer_, n_cur),
    stringsAsFactors = FALSE
  )
  
  # Future rows first (analysis focus), then current rows.
  reach_df <- dplyr::bind_rows(reach_future, reach_current)
  
  # =============================================================================
  # SECTION 11 — decision_table (one row per future dam, wide audit sheet)
  #
  # Joins every intermediate “choice” so you can see downstream pick, nearest up,
  # cascade up, and final cascade_level side by side in one table.
  # =============================================================================
  
  decision_table <- data.frame(
    future_dam_id = future_dam_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(dplyr::rename(future_trunks, future_dam_id = .data$future_dam_id), by = "future_dam_id") %>%
    dplyr::left_join(downstream_summary, by = c("future_dam_id" = "future_dam_id")) %>%
    dplyr::left_join(chosen_downstream, by = c("future_dam_id" = "future_dam_id")) %>%
    dplyr::left_join(nearest_upstream, by = c("future_dam_id" = "future_dam_id")) %>%
    dplyr::left_join(cascade_upstream, by = c("future_dam_id" = "future_dam_id")) %>%
    dplyr::left_join(dplyr::select(reach_future, dam_id, cascade_level), by = c("future_dam_id" = "dam_id"))
  
  # =============================================================================
  # SECTION 12 — debug list: expose every dataframe for step-by-step inspection
  #
  # Naming matches section flow so you can `names(dbg)` and jump to the right step.
  # =============================================================================
  
  debug <- list(
    # --- Raw network tables (no geometry) ---
    nodes_tbl = nodes_tbl,           # All graph vertices + attributes; node_id = igraph index
    edges_tbl = edges_tbl,           # All directed edges: from, to, weight, bb_id, etc.
    weights_km = weights_km,         # Edge lengths in km (same order as edges_tbl rows); passed to igraph::distances
    
    # --- Dam vertices only ---
    dam_nodes_tbl = dam_nodes_tbl,   # Nodes with non-NA dam_id (blended dams)
    current_dams_tbl = current_dams_tbl,  # Subset: is_current_dam == TRUE
    future_dams_tbl = future_dams_tbl,    # Subset: is_current_dam == FALSE
    
    # --- Trunk (bb_id) per node ---
    node_bb = node_bb,               # One representative bb_id per node_id (from incident edges)
    
    # --- Distance matrices (before melting to long tables) ---
    ds_mat = ds_mat,                 # Future → current, km along directed river (downstream from future)
    us_mat = us_mat,                 # Current → future, km (positive ⇒ current is upstream of future)
    
    # --- Index grids matching vectorized matrix order ---
    ds_grid = ds_grid,               # fut_idx × cur_idx + dist_km for ds_mat
    us_grid = us_grid,               # cur_idx × fut_idx + dist_km for us_mat
    
    # --- All future–current pairs as data frames (+ bb_id on both ends) ---
    ds_dam_all = ds_dam_all,         # Every pair, downstream direction (future → current), with dist_km
    us_dam_all = us_dam_all,         # Every pair, upstream test direction (current → future), with dist_km
    
    # --- Downstream filtering chain ---
    ds_reachable = ds_reachable,     # Pairs with finite dist_km > 0 (real downstream path)
    ds_same_trunk = ds_same_trunk,   # Subset: current_bb_id == future_bb_id (policy: downstream same trunk only)
    chosen_downstream = chosen_downstream,  # One row per future: closest same-trunk downstream current + km
    downstream_summary = downstream_summary, # Per future: has_current_downstream, min_distance_downstream_km
    
    # --- Trunk adjacency (confluence graph), not river km ---
    trunk_edges = trunk_edges,       # Undirected edges between bb_ids that meet at a node
    g_trunk = g_trunk,               # igraph object for that trunk graph
    trunk_dist_mat = trunk_dist_mat, # All-pairs hop counts between trunks (dimnames = bb_id)
    
    # --- Upstream filtering and picks ---
    us_reachable = us_reachable,     # Upstream pairs: finite dist_km > 0
    us_reachable_with_steps = us_reachable_with_steps,  # Same + trunk_step_dist (hops between trunks)
    us_within_threshold = us_within_threshold,          # Subset for cascade: dist_km <= cascade_threshold_km
    nearest_upstream = nearest_upstream,   # Per future: tiered pick (hops then km) for connectivity / dam_id_up
    cascade_upstream = cascade_upstream,     # Per future: same tier rule on threshold pool only; feeds cascade_level
    
    # --- Audit / output helpers ---
    decision_table = decision_table,       # One row per future: joins summaries + chosen up/down + cascade_level
    future_trunks = future_trunks,         # future_dam_id, node_id, bb_id
    current_trunks = current_trunks        # current_dam_id, node_id, bb_id
  )
  
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
    bb_id = character(),
    has_current_downstream = logical(),
    min_distance_downstream_km = numeric(),
    dam_id_down = character(),
    has_current_upstream = logical(),
    min_distance_upstream_km = numeric(),
    dam_id_up = character(),
    us_trunks_away = integer(),
    cascade_level = integer(),
    stringsAsFactors = FALSE
  )
}
