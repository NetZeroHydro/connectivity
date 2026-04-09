# =============================================================================
# basin_network.R
# =============================================================================
# Build net_with_dams from HydroBASINS + GRanD + FHReD + world_rivers_ffr (all
# in-memory; no read_* inside). Filter rivers by continent + Strahler on
# world_rivers_ffr (row filters shrink the table before geometry), then crop to
# the study basin polygon, then as_sfnetwork and optionally st_network_blend
# when there is at least one dam point.
#
# After connectivity_from_network(), use add_ffr_attr() to attach FFR edge
# fields (by bb_id) and optional future-dam attributes (by dam_id).
#
# Required packages: sf, dplyr, sfnetworks, tidygraph, igraph (igraph pulled by
# sfnetworks for distances if you call connectivity_from_network next).
#
# Comments below are written for walking a teammate through the flow: what each
# chunk is for, and short notes on non-obvious lines.
# =============================================================================

# -----------------------------------------------------------------------------
# REACH_DF_OUTPUT_COLS — “minimal reach table” column names (no magic).
#
# connectivity_from_network() only returns these connectivity columns (plus
# dam_type). add_ffr_attr() can *add* more columns via edge_attr_cols and
# future_extra; output_cols then *selects* which names to keep at the end.
# Edit this vector to match the slim export you want (add joined names you
# actually need, e.g. "csi", "ord_stra").
# -----------------------------------------------------------------------------
REACH_DF_OUTPUT_COLS <- c(
  "dam_id",
  "dam_type",
  "has_current_upstream",
  "has_current_downstream",
  "min_distance_upstream_km",
  "min_distance_downstream_km",
  "cascade_status",
  "bb_id"
)

#' Build snapped dam network for one HydroBASINS unit
#'
#' Takes the same four objects you would normally load in a Qmd: current GRanD
#' dams, FHReD future dams table, a regional HydroBASINS `sf`, and
#' `world_rivers_ffr`. Filters rivers by **`continent`** (column on
#' `world_rivers_ffr`) and **`ord_stra`**, then [sf::st_filter()] to the **study**
#' basin. Passing the right `continent` string (e.g. `"South America"` with
#' `hybas_sa`) avoids running the spatial join against the full global river
#' layer. Then builds a directed `sfnetwork` with edge `weight` (m) and blends
#' dam points when any dams lie in the basin.
#'
#' **Filter order (speed):** row filters `continent` + `ord_stra` first, then
#' spatial crop to `basin_poly`.
#'
#' HydroBASINS `hybas_id` is typically 10 digits: digit 1 = region; digits 2–3 =
#' Pfafstetter level (`01`–`12`, or `00` for level 0); digits 4–9 = unit id;
#' digit 10 = side in some products. Load the **shapefile whose level matches**
#' the basin you want (e.g. lev04 file for a lev04 code).
#'
#' @section Dam input columns (swapping data sources):
#'
#' **Current dams (`current_dams`, `sf`):**
#'
#' - Required: `grand_id` (renamed internally to `dam_id`); valid `sf` geometry
#'   (POINT or else coerced with [sf::st_point_on_surface()] before blending).
#' - CRS: any; transformed to the river layer CRS.
#'
#' **Future dams (`future_dams`, `data.frame`, not `sf`):**
#'
#' - Required: `lon_cleaned`, `lat_cleaned` (WGS84-style coordinates matching
#'   your cleaning pipeline); same CRS reference as `current_dams` for
#'   [sf::st_as_sf(coords = ...)].
#' - ID: `dam_id` if present; else `grand_id` or `id`; else synthetic row-number id.
#'
#' **Rivers (`world_rivers_ffr`, `sf`):**
#'
#' - Required columns: `continent`, `ord_stra`, `bb_id`, plus geometry.
#'
#' **Basins (`basins_sf`, `sf`):**
#'
#' - Required: `hybas_id` (character-coercible; must match the scalar `hybas_id`
#'   argument after [as.character()]).
#'
#' @md
#'
#' @param current_dams `sf` of current dams (e.g. GRanD after `clean_names()`),
#'   with column `grand_id` (used as `dam_id`).
#' @param future_dams `data.frame` of future dams (FHReD after `clean_names()`),
#'   **not** `sf`. Must include `lon_cleaned` and `lat_cleaned`; converted with
#'   [sf::st_as_sf()] using CRS from `current_dams` (before that object is
#'   reprojected), then transformed to the river CRS.
#' @param basins_sf Regional HydroBASINS `sf` (e.g. `hybas_sa`, `hybas_af`). The
#'   official regional products use the **same attribute schema** across continents;
#'   only which file and Pfafstetter **level** you load changes. Must contain
#'   column `hybas_id`.
#' @param world_rivers_ffr Combined rivers + FFR attributes (typically `sf`).
#'   Must include a `continent` column (see argument `continent`).
#' @param hybas_id **Single** target basin ID on `basins_sf`. Must match
#'   `hybas_id` after [as.character()].
#' @param continent Single string exactly matching values in
#'   `world_rivers_ffr$continent` for the region you are studying (e.g.
#'   `"South America"` with South American Hybas and rivers).
#' @param ord_stra_min Minimum Strahler order to keep (default `3L`).
#' @param blend_tolerance_m Passed to [sfnetworks::st_network_blend()] (metres);
#'   ignored when there are no dams in the basin (network is left unblended).
#'
#' @details Current dams must include GRanD column `grand_id` (renamed to
#'   `dam_id`). Rivers are cast to `LINESTRING`; non-POINT current dam geometry
#'   uses [sf::st_point_on_surface()] before blending. If there are **no** dams
#'   in the basin after intersection, the returned network has no blended dam
#'   nodes (no `dam_id` / `is_current_dam` on nodes from blending).
#'
#' @return Named list: `net_with_dams`, `basin_poly`, `rivers_cropped`,
#'   `dams_current_cropped`, `dams_future_cropped`, `all_dams` (`dam_id`,
#'   `is_current_dam`, geometry — may have 0 rows).
#'
#' @export
net_with_dams_from_hydrobasin <- function(
    current_dams,
    future_dams,
    basins_sf,
    world_rivers_ffr,
    hybas_id,
    continent,
    ord_stra_min = 3L,
    blend_tolerance_m = 30000) {
  
  # ---------------------------------------------------------------------------
  # net_with_dams_from_hydrobasin — what this function does (walkthrough)
  #
  # One basin id + data already in memory: (1) basin polygon in the river CRS,
  # (2) river lines filtered by continent and Strahler, then clipped to polygon,
  # (3) dams clipped to polygon, (4) directed river network with edge lengths,
  # (5) snap dams onto the network if any exist, else leave rivers-only graph.
  # ---------------------------------------------------------------------------
  
  # --- Check continent string: must be one label that exists on the river layer ---
  if (missing(continent) || length(continent) != 1L || !nzchar(as.character(continent)[1])) {
    stop(
      "continent must be a single non-empty string matching world_rivers_ffr$continent.",
      call. = FALSE
    )
  }
  continent_val <- as.character(continent)[1] # single string for dplyr filter
  
  # --- One basin only: keeps the function simple (no multi-basin union here) ---
  if (length(hybas_id) != 1L) {
    stop("hybas_id must be length 1 (single basin).", call. = FALSE)
  }
  study_id <- as.character(hybas_id)[1] # match Hybas id column after coercion
  
  if (!"hybas_id" %in% names(basins_sf)) {
    stop("basins_sf must contain column `hybas_id`.", call. = FALSE)
  }
  
  # --- Strahler cutoff: coerce so callers can pass "3" from a form if needed ---
  ord_stra_min <- suppressWarnings(as.numeric(ord_stra_min))
  if (length(ord_stra_min) != 1L || is.na(ord_stra_min)) {
    stop("ord_stra_min must be a single numeric value.", call. = FALSE)
  }
  
  # --- River layer must look like the object we expect (sf + key columns) ---
  if (!inherits(world_rivers_ffr, "sf")) {
    stop("world_rivers_ffr must be an sf object.", call. = FALSE)
  }
  if (!"continent" %in% names(world_rivers_ffr)) {
    stop("world_rivers_ffr must have a column `continent`.", call. = FALSE)
  }
  if (!"ord_stra" %in% names(world_rivers_ffr)) {
    stop("world_rivers_ffr must have a column `ord_stra`.", call. = FALSE)
  }
  if (!"bb_id" %in% names(world_rivers_ffr)) {
    stop("world_rivers_ffr must have edge attribute `bb_id` on river features.", call. = FALSE)
  }
  
  crs_r <- sf::st_crs(world_rivers_ffr) # everything else gets reprojected here
  
  # ---------------------------------------------------------------------------
  # Study basin polygon
  #
  # Pull the polygon for this hybas_id from the Hybas shapefile, fix invalid
  # shapes if any, and move it to the river CRS so intersects and lengths work.
  # If the table has duplicate rows for the same id, union them into one polygon.
  # ---------------------------------------------------------------------------
  basin_poly <- basins_sf %>%
    dplyr::filter(as.character(.data$hybas_id) == study_id) %>%
    sf::st_make_valid() %>%
    sf::st_transform(crs_r)
  
  if (nrow(basin_poly) == 0L) {
    stop(
      "No polygon in basins_sf for hybas_id: ", study_id,
      call. = FALSE
    )
  }
  if (nrow(basin_poly) > 1L) {
    u_study <- sf::st_union(sf::st_geometry(basin_poly)) # merge pieces
    basin_poly <- sf::st_sf(geometry = u_study, crs = sf::st_crs(basin_poly))
  }
  
  # ---------------------------------------------------------------------------
  # Rivers inside the study area
  #
  # First drop whole continents and tiny streams using columns only (fast).
  # Then keep lines that actually touch the basin polygon. That order avoids
  # feeding the spatial filter the entire global river set.
  # ---------------------------------------------------------------------------
  ord_min <- ord_stra_min
  rivers_crop <- world_rivers_ffr %>%
    dplyr::filter(.data$continent == continent_val, .data$ord_stra >= ord_min) %>%
    sf::st_filter(y = basin_poly, .predicate = sf::st_intersects)
  
  if (nrow(rivers_crop) == 0L) {
    stop(
      "No river reaches intersect the basin after continent = \"", continent_val,
      "\" and ord_stra >= ", ord_stra_min,
      ". Check continent spelling against unique(world_rivers_ffr$continent), hybas_id, and CRS.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Clean geometries for the network builder
  #
  # sfnetwork wants simple line features. Cast to LINESTRING and drop empties
  # so blend and shortest paths do not hit broken geometries.
  # ---------------------------------------------------------------------------
  rivers_crop <- suppressWarnings(sf::st_cast(rivers_crop, "LINESTRING"))
  rivers_crop <- rivers_crop[!sf::st_is_empty(rivers_crop), , drop = FALSE]
  
  if (nrow(rivers_crop) == 0L) {
    stop("No LINESTRING reaches left after casting; check geometries.", call. = FALSE)
  }
  
  # ---------------------------------------------------------------------------
  # Current dams (GRanD) in the basin
  #
  # Require grand_id (GRanD’s id). Reproject to river CRS, then intersect with
  # the basin polygon so only dams inside the study area are used.
  # ---------------------------------------------------------------------------
  if (!"grand_id" %in% names(current_dams)) {
    stop("current_dams must contain column `grand_id` (GRanD).", call. = FALSE)
  }
  
  crs_dam_points <- sf::st_crs(current_dams) # keep for building future sf from lon/lat
  current_dams <- sf::st_make_valid(current_dams) %>% sf::st_transform(crs_r)
  
  # ---------------------------------------------------------------------------
  # Future dams (table) in the basin
  #
  # FHReD-style tables are not sf yet: we build points from lon/lat, using the
  # same CRS as the current dams file, then reproject to rivers. dam_id is
  # created or renamed from common column names so downstream code has one id.
  # ---------------------------------------------------------------------------
  if (inherits(future_dams, "sf")) {
    stop("future_dams must be a data.frame (FHReD table), not sf.", call. = FALSE)
  }
  if (!is.data.frame(future_dams)) {
    stop("future_dams must be a data.frame.", call. = FALSE)
  }
  if (!"lon_cleaned" %in% names(future_dams) || !"lat_cleaned" %in% names(future_dams)) {
    stop("future_dams must contain columns `lon_cleaned` and `lat_cleaned`.", call. = FALSE)
  }
  
  if (!"dam_id" %in% names(future_dams)) {
    if ("grand_id" %in% names(future_dams)) {
      future_dams <- future_dams %>% dplyr::rename(dam_id = grand_id)
    } else if ("id" %in% names(future_dams)) {
      future_dams <- future_dams %>% dplyr::rename(dam_id = id)
    } else {
      future_dams <- future_dams %>%
        dplyr::mutate(dam_id = as.character(dplyr::row_number())) # fallback id
    }
  }
  
  future_dams <- sf::st_as_sf(
    future_dams,
    coords = c("lon_cleaned", "lat_cleaned"),
    crs = crs_dam_points,
    remove = FALSE # keep lon/lat columns in the table if you need them later
  ) %>%
    sf::st_make_valid() %>%
    sf::st_transform(crs_r)
  
  dams_current_cropped <- sf::st_intersection(basin_poly, current_dams) %>%
    dplyr::rename(dam_id = grand_id) %>%
    dplyr::mutate(dam_id = as.character(.data$dam_id))
  
  dams_future_cropped <- sf::st_intersection(basin_poly, future_dams) %>%
    dplyr::mutate(dam_id = as.character(.data$dam_id))
  
  # ---------------------------------------------------------------------------
  # Tag current vs future for the network blend
  #
  # Blending expects a point per dam. If intersection returned polygons or
  # lines, take a representative point. Empty dam sets still get the right
  # column types so bind_rows does not break.
  # ---------------------------------------------------------------------------
  if (nrow(dams_current_cropped) > 0L) {
    gt <- unique(as.character(sf::st_geometry_type(dams_current_cropped)))
    if (!identical(gt, "POINT")) {
      dams_current_cropped <- sf::st_point_on_surface(dams_current_cropped)
    }
  }
  
  if (nrow(dams_current_cropped) > 0L) {
    dams_current_cropped <- dams_current_cropped %>% dplyr::mutate(is_current_dam = TRUE)
  } else {
    dams_current_cropped <- dams_current_cropped %>% dplyr::mutate(is_current_dam = logical())
  }
  
  if (nrow(dams_future_cropped) > 0L) {
    dams_future_cropped <- dams_future_cropped %>% dplyr::mutate(is_current_dam = FALSE)
  } else {
    dams_future_cropped <- dams_future_cropped %>% dplyr::mutate(is_current_dam = logical())
  }
  
  # --- Minimal columns for blend: id + flag (geometry still attached on each sf) ---
  all_dams <- dplyr::bind_rows(
    dams_current_cropped %>% dplyr::select(dam_id, is_current_dam),
    dams_future_cropped %>% dplyr::select(dam_id, is_current_dam)
  )
  
  # ---------------------------------------------------------------------------
  # River network and optional dam snap
  #
  # Directed graph: water flows along edge direction. weight is length in metres
  # for downstream distance in connectivity_from_network. If there are no dam
  # points, skip blend — the graph still works for rivers only.
  # ---------------------------------------------------------------------------
  net <- sfnetworks::as_sfnetwork(rivers_crop, directed = TRUE) %>%
    tidygraph::activate("edges") %>%
    dplyr::mutate(weight = sfnetworks::edge_length())
  
  if (nrow(all_dams) > 0L) {
    net_with_dams <- net %>% sfnetworks::st_network_blend(all_dams, tolerance = blend_tolerance_m)
  } else {
    net_with_dams <- net
  }
  
  list(
    net_with_dams = net_with_dams,
    basin_poly = basin_poly,
    rivers_cropped = rivers_crop,
    dams_current_cropped = dams_current_cropped,
    dams_future_cropped = dams_future_cropped,
    all_dams = all_dams
  )
}


#' Add river-edge (FFR) and optional dam-table columns to reach_df
#'
#' @description
#' **What this function does:** `reach_df` from [connectivity_from_network()] is
#' one row per dam; each row already has a trunk id `bb_id` (river segment) and
#' `dam_id`. This helper **adds extra columns** by joining to (1) the **network
#' edge** table and/or (2) an optional **per-dam** table. It does **not** change
#' connectivity logic.
#'
#' **Three arguments, three roles:**
#'
#' \describe{
#'   \item{\code{edge_attr_cols}}{Adds columns from the \strong{river edges} in
#'     \code{net_with_dams} (cropped FFR reaches): e.g. \code{csi}, \code{ord_stra},
#'     \code{bas_name}. Join key is \code{bb_id}. This is \strong{not} the same as
#'     \code{output_cols}.}
#'   \item{\code{future_extra}}{Optional. Adds columns from a \strong{dam-level}
#'     table (e.g. FHReD) joined by \code{dam_id}. Leave \code{NULL} if you do not
#'     need extra dam attributes beyond \code{reach_df} and edge fields.}
#'   \item{\code{output_cols}}{Does not add data: after joins, \strong{selects}
#'     only these names. Often \code{c(REACH_DF_OUTPUT_COLS, "csi", "ord_stra")}.
#'     \code{NULL} keeps all columns.}
#' }
#'
#' **What is `REACH_DF_OUTPUT_COLS`?** A named vector at the top of
#' `basin_network.R` listing the **core connectivity columns** from
#' [connectivity_from_network()] (dam_id, dam_type, flags, distances,
#' cascade_status, bb_id). Use it as the base list when building `output_cols`,
#' then append any edge or `future_extra` names you want to keep.
#'
#' @md
#'
#' @param net_with_dams `sfnetwork` from [net_with_dams_from_hydrobasin()] (or
#'   equivalent) with edge column `bb_id`.
#' @param reach_df Data frame from [connectivity_from_network()] (must contain
#'   `bb_id` and `dam_id` when edge attributes are requested).
#' @param edge_attr_cols Names of columns on **network edges** to copy onto
#'   `reach_df` via `bb_id` (e.g. FFR fields). Not the same as `output_cols`.
#' @param future_extra Optional `data.frame` or `sf` (geometry dropped) with a
#'   `dam_id` column; left-joined for extra **dam-level** attributes. Often
#'   `NULL`; only needed if you want FHReD (or other) table columns on the reach
#'   table.
#' @param output_cols Optional character vector of **final** column names to
#'   **select** after joins. `NULL` keeps everything.
#'
#' @return `reach_df` with additional columns (possibly subset with `output_cols`).
#'
#' @export
add_ffr_attr <- function(
    net_with_dams,
    reach_df,
    edge_attr_cols = NULL,
    future_extra = NULL,
    output_cols = NULL) {
  
  # ---------------------------------------------------------------------------
  # add_ffr_attr — what this function does (walkthrough)
  #
  # reach_df already has one row per dam (connectivity + bb_id). Here we enrich:
  # join river/FFR columns by bb_id, optionally join extra dam table by dam_id,
  # then optionally keep only the column names you list in output_cols.
  # ---------------------------------------------------------------------------
  
  if (!inherits(net_with_dams, "sfnetwork")) {
    stop("net_with_dams must be an sfnetwork.", call. = FALSE)
  }
  if (!is.data.frame(reach_df)) {
    stop("reach_df must be a data frame.", call. = FALSE)
  }
  
  out <- reach_df # working copy we keep joining into
  
  # ---------------------------------------------------------------------------
  # Join river-segment attributes by bb_id
  #
  # Edges in the network are the same river reaches (with bb_id) after crop.
  # We collapse to one row per bb_id so the join is stable, then left_join onto
  # reach_df so every dam row picks up the FFR fields for its snapped reach.
  # ---------------------------------------------------------------------------
  if (length(edge_attr_cols) > 0L) {
    edges <- net_with_dams %>%
      tidygraph::activate("edges") %>%
      sf::st_as_sf() %>%
      sf::st_drop_geometry() %>% # only need attributes, not line shapes
      dplyr::as_tibble()
    
    if (!"bb_id" %in% names(edges)) {
      stop("Network edges must include `bb_id`.", call. = FALSE)
    }
    
    miss <- setdiff(edge_attr_cols, names(edges)) # user typos
    if (length(miss) > 0L) {
      warning(
        "edge_attr_cols not found on edges and skipped: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    cols <- intersect(edge_attr_cols, names(edges)) # only existing columns
    if (length(cols) > 0L) {
      if (!"bb_id" %in% names(out)) {
        stop("reach_df must contain column `bb_id` for edge attribute join.", call. = FALSE)
      }
      
      bb_attr <- edges %>%
        dplyr::select(bb_id, dplyr::all_of(cols)) %>%
        dplyr::filter(!is.na(.data$bb_id)) %>%
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>% # match types on join
        dplyr::group_by(.data$bb_id) %>%
        dplyr::slice(1L) %>% # one representative row per trunk id
        dplyr::ungroup()
      
      out <- out %>%
        dplyr::mutate(bb_id = as.character(.data$bb_id)) %>%
        dplyr::left_join(bb_attr, by = "bb_id")
    }
  }
  
  # ---------------------------------------------------------------------------
  # Optional dam-level table
  #
  # Same dam_id as in reach_df; left join adds columns (e.g. from FHReD) without
  # dropping dams that have no match in the extra table.
  # ---------------------------------------------------------------------------
  if (!is.null(future_extra)) {
    fe <- future_extra
    if (inherits(fe, "sf")) {
      fe <- sf::st_drop_geometry(fe) # join keys are ids, not shapes
    }
    if (!"dam_id" %in% names(fe)) {
      stop("future_extra must contain column `dam_id`.", call. = FALSE)
    }
    fe <- fe %>% dplyr::mutate(dam_id = as.character(.data$dam_id))
    out <- out %>%
      dplyr::mutate(dam_id = as.character(.data$dam_id)) %>%
      dplyr::left_join(fe, by = "dam_id", suffix = c("", ".future_extra"))
  }
  
  # ---------------------------------------------------------------------------
  # Final column subset
  #
  # If you listed output_cols, keep only those that exist; warn on typos.
  # Empty keep would wipe the table, so we only select when something matches.
  # ---------------------------------------------------------------------------
  if (length(output_cols) > 0L) {
    miss_out <- setdiff(output_cols, names(out))
    if (length(miss_out) > 0L) {
      warning(
        "output_cols not found in result and skipped: ",
        paste(miss_out, collapse = ", "),
        call. = FALSE
      )
    }
    keep <- intersect(output_cols, names(out))
    if (length(keep) > 0L) {
      out <- out %>% dplyr::select(dplyr::all_of(keep))
    }
  }
  
  out
}

