# =============================================================================
# add_ffr_attr.R
# =============================================================================
# PURPOSE
#   Enrich (append columns onto) `reach_df` from connectivity with:
#   - river-edge attributes from `net_with_dams` (join by bb_id)
#   - optional future-dam attributes from a dam table (join by dam_id)
#
# IMPORTANT CONTRACT
#   `reach_enriched` should keep EVERYTHING from `reach_df` and only ADD columns.
#   This means we do NOT drop columns via `output_cols` inside this function.
# =============================================================================

# -----------------------------------------------------------------------------
# REACH_DF_OUTPUT_COLS — core connectivity columns from `reach_df`
#
# Use this as a base “keep list” when you want a slim export after enrichment:
#   keep <- c(REACH_DF_OUTPUT_COLS, "csi", "ord_stra", "bas_name")
#   reach_enriched %>% dplyr::select(dplyr::all_of(keep))
# -----------------------------------------------------------------------------
REACH_DF_OUTPUT_COLS <- c( # define the default core columns
  "dam_id", # dam identifier
  "dam_type", # future/current label
  "has_current_upstream", # upstream presence flag
  "has_current_downstream", # downstream presence flag
  "min_distance_upstream_km", # nearest upstream distance
  "min_distance_downstream_km", # nearest downstream distance
  "dam_id_up", # nearest upstream current dam id
  "dam_id_down", # nearest same-trunk downstream current dam id
  "us_trunks_away", # trunk-graph hops to nearest upstream current dam
  "cascade_level", # cascade level (0..)
  "bb_id" # trunk id
) # end vector

#' Add river-edge (FFR) and optional future-dam columns to reach_df
#'
#' Takes the per-dam output from `connectivity_from_network()` and **only adds**
#' columns by joining:
#'
#' - **Edge (river / FFR) attributes** from the `net_with_dams` edge table by `bb_id`
#' - **Dam-level attributes** from an optional table by `dam_id`
#'
#' This function does **not** recompute connectivity and does **not** drop any
#' existing columns from `reach_df` (keep-all contract).
#'
#' ## How to select more FFR columns
#'
#' `edge_attr_cols` is the knob that controls which **edge columns** get copied
#' onto `reach_df`. Typical FFR fields include `csi`, `ord_stra`, `bas_name`, etc.
#'
#' Workflow:
#'
#' 1. Inspect available columns on cropped edges:
#'
#' \preformatted{
#' edges_tbl <- net_with_dams %>%
#'   tidygraph::activate("edges") %>%
#'   sf::st_as_sf() %>%
#'   sf::st_drop_geometry() %>%
#'   dplyr::as_tibble()
#' names(edges_tbl)
#' }
#'
#' 2. Pick the columns you want and pass them in:
#'
#' \preformatted{
#' reach_enriched <- add_ffr_attr(
#'   net_with_dams,
#'   reach_df,
#'   edge_attr_cols = c("csi", "ord_stra", "bas_name")
#' )
#' }
#'
#' If a requested name doesn't exist on the edges table, it is skipped with a warning.
#'
#' ## Keeping / dropping columns
#'
#' `output_cols` is intentionally ignored to preserve the keep-all contract. If
#' you want a slim export, do it explicitly after enrichment with `select()`.
#'
#' @md
#'
#' @param net_with_dams `sfnetwork` with edge attribute `bb_id` and any edge
#'   columns requested in `edge_attr_cols`.
#' @param reach_df Data frame from `connectivity_from_network()` (or equivalent).
#'   Must contain `dam_id`. Must contain `bb_id` when joining edge attributes.
#' @param edge_attr_cols Character vector of edge-column names to copy from
#'   `net_with_dams` edges onto `reach_df` via `bb_id`.
#' @param future_extra Optional `data.frame` or `sf` with a `dam_id` column.
#'   If `sf`, geometry is dropped before joining. Left-joined by `dam_id`.
#' @param output_cols Ignored (kept for compatibility with earlier versions).
#'
#' @return `reach_df` with additional columns appended.
#'
#' @export
add_ffr_attr <- function(
    net_with_dams,
    reach_df,
    edge_attr_cols = NULL,
    future_extra = NULL,
    output_cols = NULL) {
  
  # ---------------------------------------------------------------------------
  # SECTION 0 — start from reach_df and only add columns
  # ---------------------------------------------------------------------------
  
  out <- reach_df # working output table (we will only add columns to this)
  
  # ---------------------------------------------------------------------------
  # SECTION 1 — join edge attributes by bb_id
  #
  # This is how you add more FFR fields: list their edge-column names in
  # edge_attr_cols, and they'll be joined onto reach_df by bb_id.
  # ---------------------------------------------------------------------------
  
  if (length(edge_attr_cols) > 0L) { # only do the join if user requested columns
    
    # Pull edge attributes into a plain table (drop geometry so joins are cheap).
    edges_tbl <- net_with_dams %>% # start with sfnetwork
      tidygraph::activate("edges") %>% # activate edges
      sf::st_as_sf() %>% # convert to sf to drop geometry
      sf::st_drop_geometry() %>% # remove geometry so we only have attributes
      dplyr::as_tibble() # convert to tibble
    
    # Warn on requested columns that do not exist (typos, different schema, etc.).
    miss_edge <- setdiff(edge_attr_cols, names(edges_tbl)) # columns requested but missing
    
    if (length(miss_edge) > 0L) { # if any requested columns are missing
      warning( # warn but do not stop
        "edge_attr_cols not found on edges and skipped: ", # warning prefix
        paste(miss_edge, collapse = ", "), # list missing
        call. = FALSE # no call stack
      ) # end warning
    } # end missing edge columns warning
    
    # Only keep columns that exist on edges_tbl.
    keep_edge <- intersect(edge_attr_cols, names(edges_tbl)) # keep only columns that exist
    
    if (length(keep_edge) > 0L) { # if we have any columns to actually join
      
      # Reduce edges to one row per bb_id so the join is stable/deterministic.
      bb_attr <- edges_tbl %>% # begin with edge attributes
        dplyr::select(bb_id, dplyr::all_of(keep_edge)) %>% # keep bb_id + requested columns
        dplyr::filter(!is.na(.data$bb_id)) %>% # drop missing bb_id
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>% # ensure bb_id is character for join
        dplyr::group_by(.data$bb_id) %>% # group by bb_id
        dplyr::slice(1L) %>% # keep one representative edge row per bb_id
        dplyr::ungroup() # drop grouping
      
      # Join edge attributes onto reach_df by bb_id (type-coerce for safety).
      out <- out %>% # modify output
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>% # ensure join key type matches
        dplyr::left_join(bb_attr, by = "bb_id") # add edge attributes
      
    } # end keep_edge join
    
  } # end edge attribute join
  
  # ---------------------------------------------------------------------------
  # SECTION 2 — join future-dam attributes by dam_id (optional)
  #
  # Use future_extra when you want FHReD (or other) dam-level fields appended to
  # the reach table. This is separate from edge_attr_cols (edge join is by bb_id).
  # ---------------------------------------------------------------------------
  
  if (!is.null(future_extra)) { # only do this if user passed a table
    
    fe <- future_extra # local variable for clarity
    
    if (inherits(fe, "sf")) { # if the dam table is sf
      fe <- sf::st_drop_geometry(fe) # drop geometry before joining
    } # end sf geometry drop
    
    if (!"dam_id" %in% names(fe)) { # dam_id is required for join
      stop("future_extra must contain column `dam_id`.", call. = FALSE) # stop if missing
    } # end dam_id check
    
    fe <- fe %>% dplyr::mutate(dam_id = as.character(.data$dam_id)) # force dam_id to character
    
    out <- out %>% # modify output
      dplyr::mutate(dam_id = as.character(.data$dam_id)) %>% # ensure join key type matches
      dplyr::left_join(fe, by = "dam_id", suffix = c("", ".future_extra")) # join dam-level attributes
    
    # If capacity_mw exists, set it to NA for current dams (future-only).
    if ("capacity_mw" %in% names(out)) { # only if column exists after join
      out <- out %>% # modify output
        dplyr::mutate( # mutate columns
          capacity_mw = dplyr::if_else( # keep for future, NA for current
            .data$dam_type == "future", # condition
            .data$capacity_mw, # keep value
            NA_real_ # set to NA for current
          ) # end if_else
        ) # end mutate
    } # end capacity_mw future-only behavior
    
  } # end future_extra join
  
  # ---------------------------------------------------------------------------
  # SECTION 3 — output_cols (ignored for keep-all contract)
  # ---------------------------------------------------------------------------
  
  # output_cols is intentionally ignored so `reach_enriched` always keeps all `reach_df` columns.
  # If you want to drop columns, do it explicitly in your Qmd with dplyr::select().
  
  out # return enriched table
}

# =============================================================================
# add_ffr_attr.R
# =============================================================================
# PURPOSE
#   Enrich (append columns onto) `reach_df` from connectivity with:
#   - river-edge attributes from `net_with_dams` (join by bb_id)
#   - optional future-dam attributes from a dam table (join by dam_id)
#
# IMPORTANT CONTRACT (per your request)
#   `reach_enriched` should keep EVERYTHING from `reach_df` and only ADD columns.
#   This means we do NOT drop columns via `output_cols` inside this function.
# =============================================================================

# -----------------------------------------------------------------------------
# REACH_DF_OUTPUT_COLS — core connectivity columns from `reach_df`
# -----------------------------------------------------------------------------
REACH_DF_OUTPUT_COLS <- c( # define the default core columns
  "dam_id", # dam identifier
  "dam_type", # future/current label
  "has_current_upstream", # upstream presence flag
  "has_current_downstream", # downstream presence flag
  "min_distance_upstream_km", # nearest upstream distance
  "min_distance_downstream_km", # nearest downstream distance
  "cascade_level", # cascade level (0..)
  "bb_id" # trunk id
) # end vector

#' Add river-edge (FFR) and optional future-dam columns to reach_df
#'
#' @md
#' @export
add_ffr_attr <- function(
    net_with_dams,
    reach_df,
    edge_attr_cols = NULL,
    future_extra = NULL,
    output_cols = NULL) {
  
  # ---------------------------------------------------------------------------
  # SECTION 0 — start from reach_df and only add columns
  # ---------------------------------------------------------------------------
  
  out <- reach_df # working output table (we will only add columns to this)
  
  # ---------------------------------------------------------------------------
  # SECTION 1 — join edge attributes by bb_id
  # ---------------------------------------------------------------------------
  
  if (length(edge_attr_cols) > 0L) { # only do the join if user requested columns
    
    edges_tbl <- net_with_dams %>% # start with sfnetwork
      tidygraph::activate("edges") %>% # activate edges
      sf::st_as_sf() %>% # convert to sf to drop geometry
      sf::st_drop_geometry() %>% # remove geometry so we only have attributes
      dplyr::as_tibble() # convert to tibble
    
    miss_edge <- setdiff(edge_attr_cols, names(edges_tbl)) # columns requested but missing
    
    if (length(miss_edge) > 0L) { # if any requested columns are missing
      warning( # warn but do not stop
        "edge_attr_cols not found on edges and skipped: ", # warning prefix
        paste(miss_edge, collapse = ", "), # list missing
        call. = FALSE # no call stack
      ) # end warning
    } # end missing edge columns warning
    
    keep_edge <- intersect(edge_attr_cols, names(edges_tbl)) # keep only columns that exist
    
    if (length(keep_edge) > 0L) { # if we have any columns to actually join
      
      bb_attr <- edges_tbl %>% # begin with edge attributes
        dplyr::select(bb_id, dplyr::all_of(keep_edge)) %>% # keep bb_id + requested columns
        dplyr::filter(!is.na(.data$bb_id)) %>% # drop missing bb_id
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>% # ensure bb_id is character for join
        dplyr::group_by(.data$bb_id) %>% # group by bb_id
        dplyr::slice(1L) %>% # keep one representative edge row per bb_id
        dplyr::ungroup() # drop grouping
      
      out <- out %>% # modify output
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>% # ensure join key type matches
        dplyr::left_join(bb_attr, by = "bb_id") # add edge attributes
      
    } # end keep_edge join
    
  } # end edge attribute join
  
  # ---------------------------------------------------------------------------
  # SECTION 2 — join future-dam attributes by dam_id (optional)
  # ---------------------------------------------------------------------------
  
  if (!is.null(future_extra)) { # only do this if user passed a table
    
    fe <- future_extra # local variable for clarity
    
    if (inherits(fe, "sf")) { # if the dam table is sf
      fe <- sf::st_drop_geometry(fe) # drop geometry before joining
    } # end sf geometry drop
    
    if (!"dam_id" %in% names(fe)) { # dam_id is required for join
      stop("future_extra must contain column `dam_id`.", call. = FALSE) # stop if missing
    } # end dam_id check
    
    fe <- fe %>% dplyr::mutate(dam_id = as.character(.data$dam_id)) # force dam_id to character
    
    out <- out %>% # modify output
      dplyr::mutate(dam_id = as.character(.data$dam_id)) %>% # ensure join key type matches
      dplyr::left_join(fe, by = "dam_id", suffix = c("", ".future_extra")) # join dam-level attributes
    
    # If capacity_mw exists, set it to NA for current dams (future-only).
    if ("capacity_mw" %in% names(out)) { # only if column exists after join
      out <- out %>% # modify output
        dplyr::mutate( # mutate columns
          capacity_mw = dplyr::if_else( # keep for future, NA for current
            .data$dam_type == "future", # condition
            .data$capacity_mw, # keep value
            NA_real_ # set to NA for current
          ) # end if_else
        ) # end mutate
    } # end capacity_mw future-only behavior
    
  } # end future_extra join
  
  # ---------------------------------------------------------------------------
  # SECTION 3 — output_cols (ignored for keep-all contract)
  # ---------------------------------------------------------------------------
  
  # output_cols is intentionally ignored so `reach_enriched` always keeps all `reach_df` columns.
  # If you want to drop columns, do it explicitly in your Qmd with dplyr::select().
  
  out # return enriched table
}

