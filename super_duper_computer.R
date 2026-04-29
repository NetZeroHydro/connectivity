install.packages("sfnetworks")
install.packages("sf")
install.packages("tidygraph")
install.packages("dpylr")
install.packages("readr")
install.packages("janitor")

library(sf)
library(sfnetworks)
library(tidygraph)
library(dplyr)
library(readr)
library(janitor)

# =============================================================================
# Load Data
# =============================================================================

# Global Dam Watch
gdw <- st_read("home/lucian/hpc/netzerohydro/GDW/GDW_v1_0.gdb", layer = "GDW_barriers_v1_0") |>
  clean_names()

# HydroRIVERS and FFR
world_rivers_ffr <- readRDS("home/lucian/hpc/netzerohydro/world_rivers_ffr.rds")

# FHReD 
future_dams_raw <- read_csv("home/lucian/hpc/netzerohydro/FHReD_2015_future_dams_Zarfl_et_al_beta_version.csv") |>
  clean_names()

# =============================================================================
# Filter rivers (uneeded)
# =============================================================================

# Keep only river reaches with Strahler order >.
ord_stra_min <- 1L

# Optional subset for testing:
# - NULL = all rivers
main_riv_filter <- 41392598 

# Filter HydroRIVERS by Strahler order
rivers_filtered <- world_rivers_ffr %>% 
  dplyr::filter(.data$ord_stra > ord_stra_min)

# Optional: Filter by main_riv
if (!is.null(main_riv_filter)) {
  rivers_filtered <- rivers_filtered %>%
    dplyr::filter(as.numeric(.data$main_riv) == as.numeric(main_riv_filter))
}

# Ensure rivers are valid LINESTRING features for sfnetwork conversion.
rivers_filtered <- suppressWarnings(sf::st_cast(rivers_filtered, "LINESTRING"))
rivers_filtered <- rivers_filtered[!sf::st_is_empty(rivers_filtered), , drop = FALSE]

if (nrow(rivers_filtered) == 0L) {
  stop("No rivers left after filters. Check ord_stra/main_bas settings.", call. = FALSE)
}

# Find distance of each main_riv in KM
main_riv_dist <- rivers_filtered %>%
  dplyr::mutate(main_riv = as.character(.data$main_riv)) %>%
  dplyr::mutate(seg_len_m = as.numeric(sf::st_length(sf::st_geometry(.)))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(.data$main_riv) %>%
  dplyr::summarise(
    main_riv_dist = sum(.data$seg_len_m, na.rm = TRUE) / 1000,
    .groups = "drop"
  )

# =============================================================================
# Clean Future and Current Dams data
# =============================================================================
crs_r <- sf::st_crs(rivers_filtered)

current_dams_gdw <- gdw %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_r) %>%
  dplyr::mutate(
    dam_id = as.character(.data$gdw_id),
    is_current_dam = TRUE
  )

future_dams_sf <- future_dams_raw %>%
  sf::st_as_sf(
    coords = c("lon_cleaned", "lat_cleaned"),
    crs = 4326,
    remove = FALSE
  ) %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_r) %>%
  dplyr::mutate(
    dam_id = as.character(.data$dam_id),
    is_current_dam = FALSE
  )

# =============================================================================
# Handle duplicate Future and current dams
# =============================================================================
# Distance (meters) used to remove FHReD dams that are already represented
fhred_overlap_tolerance_m <- 500

# Keep only GDW rows explicitly labeled as FHReD-origin references.
gdw_fhred_ref <- current_dams_gdw 

#%>% dplyr::filter(.data$orig_src == "FHReD")

n_fhred_before <- nrow(future_dams_sf)

# st_is_within_distance() returns, for each FHReD row, indices of GDW ref points
# within the given tolerance. Any non-empty match means "already represented".
if (nrow(gdw_fhred_ref) > 0L && n_fhred_before > 0L) {
  overlap_idx <- sf::st_is_within_distance(
    future_dams_sf,
    gdw_fhred_ref,
    dist = fhred_overlap_tolerance_m
  )
  is_overlap <- lengths(overlap_idx) > 0L
  dams_removed <- sum(is_overlap)
  dams_future_fhred_clean <- future_dams_sf[!is_overlap, , drop = FALSE]
} else {
  dams_removed <- 0L
  dams_future_fhred_clean <- future_dams_sf
}

# =============================================================================
# Build river network (HEAVY SECTION)
# =============================================================================
# Snap tolerance in meters for st_network_blend().
blend_tolerance_m <- 500

all_dams <- dplyr::bind_rows(
  dams_future_fhred_clean %>%
    dplyr::transmute(dam_id = paste0("FHRED_", dam_id), is_current_dam, geometry = sf::st_geometry(.), project_name, capacity_mw),
  current_dams_gdw %>%
    dplyr::transmute(dam_id = paste0("GDW_", dam_id), is_current_dam, geometry = sf::st_geometry(.))
) %>% sf::st_as_sf()

# as_sfnetwork() converts river lines to a directed graph.
net <- sfnetworks::as_sfnetwork(rivers_filtered, directed = TRUE) %>%
  tidygraph::activate("edges") %>%
  # edge_length() computes geometric edge lengths; used as weighted distance input.
  dplyr::mutate(weight = sfnetworks::edge_length())

# st_network_blend() inserts/snap-aligns dam points onto nearest network edges.
net_with_dams <- sfnetworks::st_network_blend(net, all_dams, tolerance = blend_tolerance_m)


# =============================================================================
# Saving nodes and re-adding HydroRIVERS attributes
# =============================================================================
# Extract snapped network nodes and keep only dam nodes.
nodes_tbl <- net_with_dams %>%
  tidygraph::activate("nodes") %>%
  sf::st_as_sf() %>%
  dplyr::mutate(node_id = dplyr::row_number())

dams_snapped <- nodes_tbl %>%
  dplyr::filter(!is.na(.data$dam_id))

# Assign each snapped dam node to nearest filtered river reach.
# st_nearest_feature() returns row index of nearest river geometry per point.
nearest_riv_idx <- sf::st_nearest_feature(dams_snapped, rivers_filtered)

river_attr_keep <- c("hyriv_id", "bb_id", "ord_stra", "csi", "bas_name")
river_attr_keep <- river_attr_keep[river_attr_keep %in% names(rivers_filtered)]

dams_snapped <- dams_snapped %>%
  dplyr::bind_cols(
    rivers_filtered %>%
      sf::st_drop_geometry() %>%
      dplyr::slice(nearest_riv_idx) %>%
      dplyr::select(dplyr::all_of(river_attr_keep))
  )

# =============================================================================
# Save output
# =============================================================================
out <- list(
  net_with_dams = net_with_dams,
  rivers_filtered = rivers_filtered,
  dams_current_gdw = current_dams_gdw,
  dams_future_fhred_clean = dams_future_fhred_clean,
  all_dams = all_dams,
  dams_snapped = dams_snapped,
  main_riv_dist = main_riv_dist
)

saveRDS(out, file = "home/lucian/hpc/netzerohydro/out.RDS")
