# =============================================================================
# net_with_dams_from_hydrobasin.R
# =============================================================================
# Build a directed `sfnetwork` for a single HydroBASINS polygon, optionally
# blending (snapping) current + future dam points onto the network.
# =============================================================================

#' Build snapped dam network for one HydroBASINS unit
#'
#' Takes in-memory objects (no file I/O) and returns a directed `sfnetwork`
#' representing river reaches inside a single HydroBASINS polygon, with optional
#' dam nodes blended into the network.
#'
#' **Filter order (performance):**
#' \describe{
#'   \item{Attribute filter}{Filter `world_rivers_ffr` by `ord_stra >= ord_stra_min`
#'   (and `continent` if provided).}
#'   \item{Spatial crop}{Then crop to the basin polygon via `st_filter(..., st_intersects)`.}
#' }
#'
#' **Required columns / types**
#' \describe{
#'   \item{`basins_sf` (sf)}{Must contain `hybas_id` and polygon geometry.}
#'   \item{`world_rivers_ffr` (sf)}{Must contain `ord_stra`, `bb_id`, geometry (castable to LINESTRING).
#'     If `continent` filtering is used, must also contain `continent`.}
#'   \item{`current_dams` (sf)}{Must contain `grand_id` and geometry. Non-POINT results after intersection are
#'     coerced to points via `st_point_on_surface()` before blending.}
#'   \item{`future_dams` (data.frame)}{Must contain `lon_cleaned`, `lat_cleaned`. If missing `dam_id`, the function
#'     renames `grand_id`/`id` or creates a synthetic id.}
#' }
#'
#' **CRS behavior**
#' \describe{
#'   \item{Working CRS}{All spatial operations are done in the CRS of `world_rivers_ffr`.}
#'   \item{Future dam point creation}{Future dams are converted from lon/lat using the CRS of `current_dams`,
#'   then transformed into the working CRS.}
#' }
#'
#' @param current_dams `sf` of current dams (GRanD-like) with column `grand_id`.
#' @param future_dams `data.frame` of future dams (FHReD-like) with `lon_cleaned` and `lat_cleaned`.
#' @param basins_sf Regional HydroBASINS `sf` with column `hybas_id`.
#' @param world_rivers_ffr `sf` of rivers with FFR attributes; must include `ord_stra` and `bb_id`.
#' @param hybas_id Single basin id matching `basins_sf$hybas_id` (after `as.character()`).
#' @param continent Optional single string used to pre-filter rivers by `world_rivers_ffr$continent`.
#'   If missing/NULL/NA/empty, no continent filter is applied.
#' @param ord_stra_min Minimum Strahler order to keep.
#' @param blend_tolerance_m Tolerance (metres) passed to `sfnetworks::st_network_blend()`.
#'
#' @return Named list with:
#' \describe{
#'   \item{`net_with_dams`}{Directed `sfnetwork` of rivers with optional blended dam nodes.}
#'   \item{`basin_poly`}{Single-row `sf` polygon used for cropping (unioned if needed).}
#'   \item{`rivers_cropped`}{Cropped `sf` of river reaches (LINESTRING).}
#'   \item{`dams_current_cropped`}{Cropped `sf` of current dams, with `dam_id` and `is_current_dam`.}
#'   \item{`dams_future_cropped`}{Cropped `sf` of future dams, with `dam_id` and `is_current_dam`.}
#'   \item{`all_dams`}{Combined `sf` of all dams (`dam_id`, `is_current_dam`) used for blending.}
#' }
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
  
  # --- Continent is OPTIONAL (fast pre-filter; skip if missing/blank). ---     # continent rule: optional for flexibility
  use_continent <- !(missing(continent) || is.null(continent) ||               # TRUE only if continent argument exists
                       length(continent) != 1L ||                              # require scalar input when present
                       is.na(as.character(continent)[1]) ||                    # reject NA
                       !nzchar(as.character(continent)[1]))                    # reject empty string
  
  continent_val <- if (use_continent) {                                        # compute the actual filter value
    as.character(continent)[1]                                                 # coerce to scalar character
  } else {                                                                     # if not filtering by continent
    NA_character_                                                              # sentinel meaning “no continent filter”
  }                                                                            # end continent_val
  
  # --- Ensure `hybas_id` is a single basin id. ---                             # this function builds ONE basin
  if (length(hybas_id) != 1L) {                                                # enforce scalar basin id
    stop("hybas_id must be length 1 (single basin).", call. = FALSE)           # fail fast with readable error
  }                                                                            # end hybas_id check
  study_id <- as.character(hybas_id)[1]                                        # coerce basin id to character for matching
  
  # --- Ensure basin table has the id column. ---                               # basins must have a hybas id field
  if (!"hybas_id" %in% names(basins_sf)) {                                     # require hybas_id column
    stop("basins_sf must contain column `hybas_id`.", call. = FALSE)           # stop if missing
  }                                                                            # end basins_sf check
  
  # --- Coerce Strahler cutoff to numeric. ---                                  # callers may pass "3" from UI/forms
  ord_stra_min <- suppressWarnings(as.numeric(ord_stra_min))                   # coerce (quietly)
  if (length(ord_stra_min) != 1L || is.na(ord_stra_min)) {                     # ensure scalar numeric
    stop("ord_stra_min must be a single numeric value.", call. = FALSE)        # stop if invalid
  }                                                                            # end ord_stra_min validation
  
  # --- Validate river layer shape and required columns. ---                    # protect downstream dplyr/sf calls
  if (!inherits(world_rivers_ffr, "sf")) {                                     # rivers must be an sf object
    stop("world_rivers_ffr must be an sf object.", call. = FALSE)              # stop if not
  }                                                                            # end sf check
  if (use_continent && !"continent" %in% names(world_rivers_ffr)) {            # only required when we actually filter
    stop("world_rivers_ffr must have a column `continent` when continent filtering is used.", call. = FALSE) # stop
  }                                                                            # end conditional continent col check
  if (!"ord_stra" %in% names(world_rivers_ffr)) {                              # need Strahler column for filtering
    stop("world_rivers_ffr must have a column `ord_stra`.", call. = FALSE)     # stop if missing
  }                                                                            # end ord_stra col check
  if (!"bb_id" %in% names(world_rivers_ffr)) {                                 # need trunk id for downstream joins
    stop("world_rivers_ffr must have edge attribute `bb_id` on river features.", call. = FALSE) # stop if missing
  }                                                                            # end bb_id col check
  
  # --- Choose the CRS of the river layer as the working CRS. ---               # all layers will be transformed here
  crs_r <- sf::st_crs(world_rivers_ffr)                                        # capture rivers CRS (authoritative)
  
  # --------------------------------------------------------------------------- # begin basin polygon build
  # Basin polygon (study area)                                                  # pull the basin polygon for one hybas id
  # --------------------------------------------------------------------------- # end basin polygon header
  basin_poly <- basins_sf %>%                                                  # start from basins layer
    dplyr::filter(as.character(.data$hybas_id) == study_id) %>%                # keep the one matching id
    sf::st_make_valid() %>%                                                    # fix any invalid polygon geometry
    sf::st_transform(crs_r)                                                    # transform to rivers CRS
  
  if (nrow(basin_poly) == 0L) {                                                # if basin id not found
    stop("No polygon in basins_sf for hybas_id: ", study_id, call. = FALSE)    # stop with useful message
  }                                                                            # end missing basin stop
  
  if (nrow(basin_poly) > 1L) {                                                 # if multiple parts/rows share same id
    u_study <- sf::st_union(sf::st_geometry(basin_poly))                       # union geometries into one multipart polygon
    basin_poly <- sf::st_sf(geometry = u_study, crs = sf::st_crs(basin_poly))  # rebuild single-row sf with CRS
  }                                                                            # end union for multi-piece basin
  
  # --------------------------------------------------------------------------- # begin river filtering/cropping
  # Rivers: filter fast by attributes, then spatially crop to basin polygon      # attribute filters reduce cost of st_filter
  # --------------------------------------------------------------------------- # end river header
  rivers_base <- world_rivers_ffr %>%                                          # start from all rivers
    dplyr::filter(.data$ord_stra >= ord_stra_min)                              # apply Strahler threshold first (fast)
  
  if (use_continent) {                                                         # if the caller provided a continent
    rivers_base <- rivers_base %>%                                             # continue filtering
      dplyr::filter(.data$continent == continent_val)                          # apply continent row filter (fast)
  }                                                                            # end optional continent filter
  
  rivers_crop <- rivers_base %>%                                               # take attribute-filtered rivers
    sf::st_filter(y = basin_poly, .predicate = sf::st_intersects)              # keep only those touching the basin polygon
  
  if (nrow(rivers_crop) == 0L) {                                               # if nothing intersects after filters
    cont_msg <- if (use_continent) paste0("continent = \"", continent_val, "\" and ") else "" # build readable message
    stop(                                                                      # stop with actionable guidance
      "No river reaches intersect the basin after ", cont_msg,                 # show continent only when used
      "ord_stra >= ", ord_stra_min,                                            # show threshold used
      ". Check hybas_id and CRS (and continent spelling if provided).",        # hints
      call. = FALSE                                                            # suppress call stack
    )                                                                          # end stop call
  }                                                                            # end no-rivers stop
  
  rivers_crop <- suppressWarnings(sf::st_cast(rivers_crop, "LINESTRING"))      # ensure LINESTRING geometries for sfnetwork
  rivers_crop <- rivers_crop[!sf::st_is_empty(rivers_crop), , drop = FALSE]    # drop empty geometries after casting
  
  if (nrow(rivers_crop) == 0L) {                                               # if all became empty after casting
    stop("No LINESTRING reaches left after casting; check geometries.", call. = FALSE) # stop with geometry hint
  }                                                                            # end empty-after-cast stop
  
  # --------------------------------------------------------------------------- # begin current dams handling
  # Current dams: validate and crop to basin polygon                             # crop GRanD points to study basin
  # --------------------------------------------------------------------------- # end current dams header
  if (!"grand_id" %in% names(current_dams)) {                                  # require GRanD id column
    stop("current_dams must contain column `grand_id` (GRanD).", call. = FALSE) # stop if missing
  }                                                                            # end current_dams id check
  
  crs_dam_points <- sf::st_crs(current_dams)                                   # store input dam CRS for future lon/lat conversion
  current_dams <- sf::st_make_valid(current_dams) %>%                          # fix invalid dam geometry (if any)
    sf::st_transform(crs_r)                                                    # transform current dams to rivers CRS
  
  # --------------------------------------------------------------------------- # begin future dams handling
  # Future dams: ensure table, build sf points from lon/lat, then crop           # convert FHReD-style table to sf points
  # --------------------------------------------------------------------------- # end future dams header
  if (inherits(future_dams, "sf")) {                                           # reject if already sf (we expect data.frame)
    stop("future_dams must be a data.frame (FHReD table), not sf.", call. = FALSE) # stop if wrong type
  }                                                                            # end sf rejection
  if (!is.data.frame(future_dams)) {                                           # ensure table-like input
    stop("future_dams must be a data.frame.", call. = FALSE)                   # stop if wrong type
  }                                                                            # end data.frame check
  if (!"lon_cleaned" %in% names(future_dams) || !"lat_cleaned" %in% names(future_dams)) { # require cleaned coords
    stop("future_dams must contain columns `lon_cleaned` and `lat_cleaned`.", call. = FALSE) # stop if missing
  }                                                                            # end coords check
  
  if (!"dam_id" %in% names(future_dams)) {                                     # if no standard id column
    if ("grand_id" %in% names(future_dams)) {                                  # prefer grand_id if present
      future_dams <- future_dams %>%                                           # transform the table
        dplyr::rename(dam_id = grand_id)                                       # rename to dam_id
    } else if ("id" %in% names(future_dams)) {                                 # else accept generic id
      future_dams <- future_dams %>%                                           # transform the table
        dplyr::rename(dam_id = id)                                             # rename to dam_id
    } else {                                                                   # else create a synthetic id
      future_dams <- future_dams %>%                                           # transform the table
        dplyr::mutate(dam_id = as.character(dplyr::row_number()))              # assign row-number ids as character
    }                                                                          # end id fallback
  }                                                                            # end dam_id creation/rename
  
  future_dams <- sf::st_as_sf(                                                 # convert table into sf points
    future_dams,                                                               # input FHReD table
    coords = c("lon_cleaned", "lat_cleaned"),                                  # coordinate columns
    crs = crs_dam_points,                                                      # interpret coords in dam-file CRS
    remove = FALSE                                                             # keep lon/lat columns in the table
  ) %>%                                                                        # continue pipeline
    sf::st_make_valid() %>%                                                    # fix any invalid point geometry
    sf::st_transform(crs_r)                                                    # transform to rivers CRS
  
  dams_current_cropped <- sf::st_intersection(basin_poly, current_dams) %>%    # clip current dams to basin polygon
    dplyr::rename(dam_id = grand_id) %>%                                       # standardize id name to dam_id
    dplyr::mutate(dam_id = as.character(.data$dam_id))                         # ensure join-safe character ids
  
  dams_future_cropped <- sf::st_intersection(basin_poly, future_dams) %>%      # clip future dams to basin polygon
    dplyr::mutate(dam_id = as.character(.data$dam_id))                         # ensure join-safe character ids
  
  # --------------------------------------------------------------------------- # begin blending prep
  # Tag current vs future for blending                                           # blending keeps attributes on inserted nodes
  # --------------------------------------------------------------------------- # end blending prep header
  if (nrow(dams_current_cropped) > 0L) {                                       # only if there are current dams
    gt <- unique(as.character(sf::st_geometry_type(dams_current_cropped)))     # check geometry type(s)
    if (!identical(gt, "POINT")) {                                             # if intersection returned non-point geometry
      dams_current_cropped <- sf::st_point_on_surface(dams_current_cropped)    # coerce to representative points
    }                                                                          # end point coercion
  }                                                                            # end current dam geometry fix
  
  if (nrow(dams_current_cropped) > 0L) {                                       # if any current dams remain
    dams_current_cropped <- dams_current_cropped %>%                           # mutate the sf object
      dplyr::mutate(is_current_dam = TRUE)                                     # tag as current
  } else {                                                                     # if none remain
    dams_current_cropped <- dams_current_cropped %>%                           # mutate the empty sf object
      dplyr::mutate(is_current_dam = logical())                                # preserve column type on empty
  }                                                                            # end current dam tagging
  
  if (nrow(dams_future_cropped) > 0L) {                                        # if any future dams remain
    dams_future_cropped <- dams_future_cropped %>%                             # mutate the sf object
      dplyr::mutate(is_current_dam = FALSE)                                    # tag as future
  } else {                                                                     # if none remain
    dams_future_cropped <- dams_future_cropped %>%                             # mutate the empty sf object
      dplyr::mutate(is_current_dam = logical())                                # preserve column type on empty
  }                                                                            # end future dam tagging
  
  all_dams <- dplyr::bind_rows(                                                # combine current + future into one sf object
    dams_current_cropped %>% dplyr::select(dam_id, is_current_dam),             # keep minimal columns for blending
    dams_future_cropped %>% dplyr::select(dam_id, is_current_dam)               # keep minimal columns for blending
  )                                                                            # end bind_rows
  
  # --------------------------------------------------------------------------- # begin network build
  # Build directed sfnetwork and blend dams (if any)                             # build graph and optionally snap dams
  # --------------------------------------------------------------------------- # end network header
  net <- sfnetworks::as_sfnetwork(rivers_crop, directed = TRUE) %>%            # convert cropped rivers into directed sfnetwork
    tidygraph::activate("edges") %>%                                           # switch to edge table
    dplyr::mutate(weight = sfnetworks::edge_length())                          # compute edge length (metres) as weight
  
  if (nrow(all_dams) > 0L) {                                                   # only blend when there are any dam points
    net_with_dams <- net %>%                                                   # start from rivers-only network
      sfnetworks::st_network_blend(all_dams, tolerance = blend_tolerance_m)    # snap dams into the network as nodes
  } else {                                                                     # if there are no dams
    net_with_dams <- net                                                       # keep rivers-only network unchanged
  }                                                                            # end blend switch
  
  list(                                                                        # return a named list for downstream steps
    net_with_dams = net_with_dams,                                             # sfnetwork with optional blended dam nodes
    basin_poly = basin_poly,                                                   # basin polygon used for cropping
    rivers_cropped = rivers_crop,                                              # cropped river sf used to build network
    dams_current_cropped = dams_current_cropped,                               # current dams within basin
    dams_future_cropped = dams_future_cropped,                                 # future dams within basin
    all_dams = all_dams                                                        # combined dam sf used for blending
  )                                                                            # end return list
}                                                                              # end function
