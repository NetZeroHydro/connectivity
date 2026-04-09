# =============================================================================
# net_with_dams_from_hydrobasin.R
# =============================================================================
# PURPOSE
#   Build `net_with_dams` for ONE HydroBASINS unit by:
#   - filtering rivers by continent + Strahler
#   - cropping rivers to the basin polygon
#   - cropping current + future dams to the basin polygon
#   - building a directed sfnetwork and adding edge lengths as `weight`
#   - blending dam points onto the network (if any dams exist)
# =============================================================================

#' Build snapped dam network for one HydroBASINS unit
#'
#' @md
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
  
  # Treat continent as optional: NULL/NA/"" means “don’t filter by continent”
  use_continent <- !(missing(continent) || is.null(continent) || is.na(continent) || !nzchar(as.character(continent)[1]))
  
  if (use_continent) {
    continent_val <- as.character(continent)[1]
    rivers_crop <- world_rivers_ffr %>%
      dplyr::filter(.data$continent == continent_val, .data$ord_stra >= ord_stra_min) %>%
      sf::st_filter(y = basin_poly, .predicate = sf::st_intersects)
  } else {
    rivers_crop <- world_rivers_ffr %>%
      dplyr::filter(.data$ord_stra >= ord_stra_min) %>%
      sf::st_filter(y = basin_poly, .predicate = sf::st_intersects)
  }
  
  # --- Ensure `hybas_id` is a single basin id. ---
  if (length(hybas_id) != 1L) { # enforce scalar basin id
    stop("hybas_id must be length 1 (single basin).", call. = FALSE) # error if not
  } # end hybas_id length check
  study_id <- as.character(hybas_id)[1] # convert basin id to character for matching
  
  # --- Ensure basin table has the id column. ---
  if (!"hybas_id" %in% names(basins_sf)) { # check required column exists
    stop("basins_sf must contain column `hybas_id`.", call. = FALSE) # stop if missing
  } # end basins_sf column check
  
  # --- Coerce Strahler cutoff to numeric. ---
  ord_stra_min <- suppressWarnings(as.numeric(ord_stra_min)) # allow callers to pass "3"
  if (length(ord_stra_min) != 1L || is.na(ord_stra_min)) { # validate scalar numeric
    stop("ord_stra_min must be a single numeric value.", call. = FALSE) # stop if bad
  } # end ord_stra_min validation
  
  # --- Validate river layer shape and columns. ---
  if (!inherits(world_rivers_ffr, "sf")) { # rivers must be sf
    stop("world_rivers_ffr must be an sf object.", call. = FALSE) # stop if not
  } # end sf check
  if (!"continent" %in% names(world_rivers_ffr)) { # continent column needed for fast filter
    stop("world_rivers_ffr must have a column `continent`.", call. = FALSE) # stop if missing
  } # end continent col check
  if (!"ord_stra" %in% names(world_rivers_ffr)) { # Strahler column needed
    stop("world_rivers_ffr must have a column `ord_stra`.", call. = FALSE) # stop if missing
  } # end ord_stra col check
  if (!"bb_id" %in% names(world_rivers_ffr)) { # trunk id column needed
    stop("world_rivers_ffr must have edge attribute `bb_id` on river features.", call. = FALSE) # stop if missing
  } # end bb_id col check
  
  # --- Choose the CRS of the river layer as the working CRS. ---
  crs_r <- sf::st_crs(world_rivers_ffr) # read CRS from rivers
  
  # ---------------------------------------------------------------------------
  # Basin polygon (study area)
  # ---------------------------------------------------------------------------
  
  # Select the basin polygon for the requested basin id.
  basin_poly <- basins_sf %>% # start with all basins
    dplyr::filter(as.character(.data$hybas_id) == study_id) %>% # keep the one matching id
    sf::st_make_valid() %>% # fix invalid polygon geometries
    sf::st_transform(crs_r) # transform basin polygon to river CRS
  
  # Stop if basin id was not found.
  if (nrow(basin_poly) == 0L) { # no polygon rows
    stop("No polygon in basins_sf for hybas_id: ", study_id, call. = FALSE) # stop with message
  } # end no-basin check
  
  # If basin id has multiple pieces, union them into one geometry.
  if (nrow(basin_poly) > 1L) { # multiple rows for same id
    u_study <- sf::st_union(sf::st_geometry(basin_poly)) # union geometries into one
    basin_poly <- sf::st_sf(geometry = u_study, crs = sf::st_crs(basin_poly)) # rebuild sf with one row
  } # end multi-piece basin union
  
  # ---------------------------------------------------------------------------
  # Rivers: filter fast by attributes, then spatially crop to basin polygon
  # ---------------------------------------------------------------------------
  
  # Apply attribute filters first (continent + ord_stra).
  rivers_crop <- world_rivers_ffr %>% # start with all rivers
    dplyr::filter(.data$continent == continent_val, .data$ord_stra >= ord_stra_min) %>% # attribute filter
    sf::st_filter(y = basin_poly, .predicate = sf::st_intersects) # spatial crop to basin polygon
  
  # Stop if no rivers remain after cropping.
  if (nrow(rivers_crop) == 0L) { # no river reaches
    stop( # stop with guidance
      "No river reaches intersect the basin after continent = \"", continent_val, # show continent
      "\" and ord_stra >= ", ord_stra_min, # show threshold
      ". Check continent spelling, hybas_id, and CRS.", # hint
      call. = FALSE # no call stack
    ) # end stop
  } # end no-rivers check
  
  # Ensure river geometries are LINESTRING and drop empty features.
  rivers_crop <- suppressWarnings(sf::st_cast(rivers_crop, "LINESTRING")) # cast to lines
  rivers_crop <- rivers_crop[!sf::st_is_empty(rivers_crop), , drop = FALSE] # drop empty geometries
  
  # Stop if all rivers became empty after casting.
  if (nrow(rivers_crop) == 0L) { # no lines left
    stop("No LINESTRING reaches left after casting; check geometries.", call. = FALSE) # stop
  } # end empty-after-cast check
  
  # ---------------------------------------------------------------------------
  # Current dams: validate and crop to basin polygon
  # ---------------------------------------------------------------------------
  
  # Ensure GRanD id column exists so we can create `dam_id`.
  if (!"grand_id" %in% names(current_dams)) { # check column exists
    stop("current_dams must contain column `grand_id` (GRanD).", call. = FALSE) # stop if missing
  } # end current_dams id check
  
  # Remember the original CRS of current dams (used as CRS for future dams lon/lat conversion).
  crs_dam_points <- sf::st_crs(current_dams) # store CRS
  
  # Make current dam geometries valid and transform them to river CRS.
  current_dams <- sf::st_make_valid(current_dams) %>% # fix invalid geometry
    sf::st_transform(crs_r) # transform to river CRS
  
  # ---------------------------------------------------------------------------
  # Future dams: ensure table, build sf points from lon/lat, then crop
  # ---------------------------------------------------------------------------
  
  # Reject future_dams if it is already sf (we expect a data.frame).
  if (inherits(future_dams, "sf")) { # future_dams should not be sf here
    stop("future_dams must be a data.frame (FHReD table), not sf.", call. = FALSE) # stop
  } # end sf rejection
  
  # Ensure future_dams is a data.frame.
  if (!is.data.frame(future_dams)) { # validate type
    stop("future_dams must be a data.frame.", call. = FALSE) # stop
  } # end df check
  
  # Ensure lon/lat columns exist.
  if (!"lon_cleaned" %in% names(future_dams) || !"lat_cleaned" %in% names(future_dams)) { # require coords
    stop("future_dams must contain columns `lon_cleaned` and `lat_cleaned`.", call. = FALSE) # stop
  } # end coords check
  
  # Ensure `dam_id` exists on the future dams table (create if missing).
  if (!"dam_id" %in% names(future_dams)) { # no dam_id provided
    if ("grand_id" %in% names(future_dams)) { # try grand_id
      future_dams <- future_dams %>% dplyr::rename(dam_id = grand_id) # rename to dam_id
    } else if ("id" %in% names(future_dams)) { # try generic id
      future_dams <- future_dams %>% dplyr::rename(dam_id = id) # rename to dam_id
    } else { # no usable id column
      future_dams <- future_dams %>% dplyr::mutate(dam_id = as.character(dplyr::row_number())) # create synthetic id
    } # end id fallback
  } # end dam_id creation
  
  # Convert future_dams lon/lat into an sf point layer.
  future_dams <- sf::st_as_sf( # make sf
    future_dams, # input table
    coords = c("lon_cleaned", "lat_cleaned"), # coordinate columns
    crs = crs_dam_points, # CRS to interpret lon/lat in (matches your input dam file CRS)
    remove = FALSE # keep lon/lat columns
  ) %>% # continue pipeline
    sf::st_make_valid() %>% # fix invalid geometries
    sf::st_transform(crs_r) # transform to river CRS
  
  # Crop current dams to the basin polygon and rename grand_id -> dam_id.
  dams_current_cropped <- sf::st_intersection(basin_poly, current_dams) %>% # intersect with basin
    dplyr::rename(dam_id = grand_id) %>% # rename id
    dplyr::mutate(dam_id = as.character(.data$dam_id)) # make id character
  
  # Crop future dams to the basin polygon (dam_id already exists).
  dams_future_cropped <- sf::st_intersection(basin_poly, future_dams) %>% # intersect with basin
    dplyr::mutate(dam_id = as.character(.data$dam_id)) # make id character
  
  # ---------------------------------------------------------------------------
  # Tag current vs future for blending
  # ---------------------------------------------------------------------------
  
  # If current dams were returned with non-POINT geometries, convert to representative points.
  if (nrow(dams_current_cropped) > 0L) { # only if there are current dams
    gt <- unique(as.character(sf::st_geometry_type(dams_current_cropped))) # geometry types present
    if (!identical(gt, "POINT")) { # not already points
      dams_current_cropped <- sf::st_point_on_surface(dams_current_cropped) # convert to points
    } # end geometry coercion
  } # end current dam geometry fix
  
  # Add `is_current_dam = TRUE` for current dams, or keep empty logical column if none.
  if (nrow(dams_current_cropped) > 0L) { # if we have any current dams
    dams_current_cropped <- dams_current_cropped %>% dplyr::mutate(is_current_dam = TRUE) # tag current
  } else { # no current dams
    dams_current_cropped <- dams_current_cropped %>% dplyr::mutate(is_current_dam = logical()) # empty logical
  } # end current tag
  
  # Add `is_current_dam = FALSE` for future dams, or keep empty logical column if none.
  if (nrow(dams_future_cropped) > 0L) { # if we have any future dams
    dams_future_cropped <- dams_future_cropped %>% dplyr::mutate(is_current_dam = FALSE) # tag future
  } else { # no future dams
    dams_future_cropped <- dams_future_cropped %>% dplyr::mutate(is_current_dam = logical()) # empty logical
  } # end future tag
  
  # Combine current and future dams into one sf object with minimal columns for blending.
  all_dams <- dplyr::bind_rows( # stack the two sf objects
    dams_current_cropped %>% dplyr::select(dam_id, is_current_dam), # current dam_id + flag
    dams_future_cropped %>% dplyr::select(dam_id, is_current_dam) # future dam_id + flag
  ) # end bind_rows
  
  # ---------------------------------------------------------------------------
  # Build directed sfnetwork and blend dams (if any)
  # ---------------------------------------------------------------------------
  
  # Build a directed river network.
  net <- sfnetworks::as_sfnetwork(rivers_crop, directed = TRUE) %>% # create directed sfnetwork
    tidygraph::activate("edges") %>% # switch to edge table
    dplyr::mutate(weight = sfnetworks::edge_length()) # compute edge length in metres
  
  # Blend dams onto the network if any dam points exist.
  if (nrow(all_dams) > 0L) { # only blend when there are dam points
    net_with_dams <- net %>% sfnetworks::st_network_blend(all_dams, tolerance = blend_tolerance_m) # snap dams
  } else { # no dams to blend
    net_with_dams <- net # keep the rivers-only network
  } # end blend
  
  # Return a named list of outputs used by downstream steps.
  list( # return outputs
    net_with_dams = net_with_dams, # sfnetwork with optional dam nodes
    basin_poly = basin_poly, # basin polygon used for cropping
    rivers_cropped = rivers_crop, # cropped rivers sf
    dams_current_cropped = dams_current_cropped, # cropped current dams sf
    dams_future_cropped = dams_future_cropped, # cropped future dams sf
    all_dams = all_dams # combined dam sf used for blending
  ) # end return
}

