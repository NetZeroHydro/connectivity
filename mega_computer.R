# install.packages("sfnetworks")
# install.packages("sf")
# install.packages("tidygraph")
# install.packages("dpylr")
# install.packages("readr")
# install.packages("janitor")
# install.packages("igraph")


  library(sf)
  library(sfnetworks)
  library(tidygraph)
  library(dplyr)
  library(readr)
  library(janitor)
  library(igraph)

# =============================================================================
# CONFIG
# =============================================================================

# Input paths (absolute paths with leading /)
path_gdw_gdb <- "/capstone/netzerohydro/data/raw/GDW/GDW_v1_0.gdb"
layer_gdw <- "GDW_barriers_v1_0"
path_rivers_rds <- "/capstone/netzerohydro/data/cleaned_v1/world_rivers_ffr.rds"
path_fhred_csv <- "/capstone/netzerohydro/data/raw/FHReD_2015_future_dams/FHReD_2015_future_dams_Zarfl_et_al_beta_version.csv"

# Optional filter.
main_riv_filter <- NULL
ord_stra_min <- 1L

# Dedupe / blending tolerances (meters)
fhred_overlap_tolerance_m <- 500
blend_tolerance_m <- 500

# Thresholds (numeric)
thresholds_km <- c(10, 20, 50, 100, 250, Inf)

# Enrichment edge attrs
edge_attr_cols <- c("csi", "bas_name", "main_riv", "hyriv_id")

# Output
output_dir <- "/capstone/netzerohydro/data/processed_data"
out_net_path <- file.path(output_dir, "out_world_net_with_dams.rds")
write_combined <- TRUE

# Safety test (how many main_rivs to run) NULL for world
max_main_riv <- 5

# =============================================================================
# FUNCTIONS
# =============================================================================
source("connectivity_function_v3.R")
source("add_ffr_attr.R")

# =============================================================================
# Load Data
# =============================================================================

message("Loading GDW...")
gdw <- sf::st_read(path_gdw_gdb, layer = layer_gdw, quiet = TRUE) |>
  janitor::clean_names()

message("Loading world_rivers_ffr...")
world_rivers_ffr <- readRDS(path_rivers_rds)

message("Loading FHReD...")
future_dams_raw <- readr::read_csv(path_fhred_csv, show_col_types = FALSE) |>
  janitor::clean_names()

# =============================================================================
# Filter rivers
# =============================================================================

rivers_filtered <- world_rivers_ffr %>%
  dplyr::filter(.data$ord_stra > ord_stra_min)

if (!is.null(main_riv_filter)) {
  rivers_filtered <- rivers_filtered %>%
    dplyr::filter(as.numeric(.data$main_riv) == as.numeric(main_riv_filter))
}

rivers_filtered <- suppressWarnings(sf::st_cast(rivers_filtered, "LINESTRING"))
rivers_filtered <- rivers_filtered[!sf::st_is_empty(rivers_filtered), , drop = FALSE]

if (nrow(rivers_filtered) == 0L) {
  stop("No rivers left after filters. Check ord_stra/main_riv settings.", call. = FALSE)
}

if (!"main_riv" %in% names(rivers_filtered)) {
  stop("`main_riv` is required on rivers for world-scale partitioning.", call. = FALSE)
}

main_riv_dist <- rivers_filtered %>%
  dplyr::mutate(main_riv = as.character(.data$main_riv)) %>%
  dplyr::mutate(seg_len_m = as.numeric(sf::st_length(sf::st_geometry(.)))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(.data$main_riv) %>%
  dplyr::summarise(main_riv_dist = sum(.data$seg_len_m, na.rm = TRUE) / 1000, .groups = "drop")

# =============================================================================
# Clean Future and Current Dams data
# =============================================================================

crs_r <- sf::st_crs(rivers_filtered)

current_dams_gdw <- gdw %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_r) %>%
  dplyr::mutate(
    dam_id_raw = as.character(.data$gdw_id),
    dam_id = paste0("GDW_", .data$dam_id_raw),
    is_current_dam = TRUE
  )

future_dams_sf <- future_dams_raw %>%
  sf::st_as_sf(coords = c("lon_cleaned", "lat_cleaned"), crs = 4326, remove = FALSE) %>%
  sf::st_make_valid() %>%
  sf::st_transform(crs_r) %>%
  dplyr::mutate(
    dam_id_raw = as.character(.data$dam_id),
    dam_id = paste0("FHRED_", .data$dam_id_raw),
    is_current_dam = FALSE
  )

# =============================================================================
# Handle duplicate Future and current dams
# =============================================================================

gdw_fhred_ref <- current_dams_gdw
n_fhred_before <- nrow(future_dams_sf)

if (nrow(gdw_fhred_ref) > 0L && n_fhred_before > 0L) {
  overlap_idx <- sf::st_is_within_distance(future_dams_sf, gdw_fhred_ref, dist = fhred_overlap_tolerance_m)
  is_overlap <- lengths(overlap_idx) > 0L
  dams_removed <- sum(is_overlap)
  dams_future_fhred_clean <- future_dams_sf[!is_overlap, , drop = FALSE]
} else {
  dams_removed <- 0L
  dams_future_fhred_clean <- future_dams_sf
}

message("FHReD removed by overlap: ", dams_removed)

# =============================================================================
# Build dam extra table used by add_ffr_attr()
# =============================================================================

fhred_extra <- future_dams_raw %>%
  dplyr::transmute(
    dam_id = paste0("FHRED_", as.character(dam_id)),
    dam_name = .data$project_name,
    capacity_mw = as.numeric(.data$capacity_mw),
    main_river = .data$main_river,
    major_basin = dplyr::coalesce(.data$major_basin, NA_character_),
    hybas_l12 = NA_character_
  )

gdw_extra <- gdw %>%
  sf::st_drop_geometry() %>%
  dplyr::transmute(
    dam_id = paste0("GDW_", as.character(.data$gdw_id)),
    dam_name = .data$dam_name,
    capacity_mw = as.numeric(.data$power_mw),
    main_river = NA_character_,
    major_basin = dplyr::coalesce(.data$main_basin, NA_character_),
    hybas_l12 = as.character(.data$hybas_l12)
  )

dam_extra_tbl <- dplyr::bind_rows(fhred_extra, gdw_extra) %>%
  dplyr::mutate(dam_id = trimws(as.character(.data$dam_id))) %>%
  dplyr::distinct(.data$dam_id, .keep_all = TRUE)

# =============================================================================
# Iterate by main_bas + thresholds
# =============================================================================

main_riv_values <- rivers_filtered %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(.data$main_riv) %>%
  dplyr::pull(.data$main_riv)

main_riv_values <- main_riv_values[!is.na(main_riv_values)]
if (!is.null(max_main_riv)) {
  main_riv_values <- head(main_riv_values, max_main_riv)
}

message("Total main_riv partitions: ", length(main_riv_values))

results_by_threshold <- stats::setNames(vector("list", length(thresholds_km)), as.character(thresholds_km))

for (r in main_riv_values) {
  message("Processing main_riv: ", r)
  
  rivers_b <- rivers_filtered %>% dplyr::filter(.data$main_riv == r)
  if (nrow(rivers_b) == 0L) next
  
  net_b <- sfnetworks::as_sfnetwork(rivers_b, directed = TRUE) %>%
    tidygraph::activate("edges") %>%
    dplyr::mutate(weight = sfnetworks::edge_length())
  
  # subset dams near this basin's rivers
  idx_cur <- sf::st_nearest_feature(current_dams_gdw, rivers_b)
  idx_fut <- sf::st_nearest_feature(dams_future_fhred_clean, rivers_b)
  
  nearest_cur <- suppressWarnings(sf::st_distance(current_dams_gdw, rivers_b[idx_cur, ], by_element = TRUE))
  nearest_fut <- suppressWarnings(sf::st_distance(dams_future_fhred_clean, rivers_b[idx_fut, ], by_element = TRUE))
  
  dams_cur_b <- current_dams_gdw[as.numeric(nearest_cur) <= blend_tolerance_m, , drop = FALSE]
  dams_fut_b <- dams_future_fhred_clean[as.numeric(nearest_fut) <= blend_tolerance_m, , drop = FALSE]
  
  all_dams_b <- dplyr::bind_rows(
    dams_fut_b %>% dplyr::transmute(dam_id = .data$dam_id, is_current_dam = .data$is_current_dam, geometry = sf::st_geometry(.)),
    dams_cur_b %>% dplyr::transmute(dam_id = .data$dam_id, is_current_dam = .data$is_current_dam, geometry = sf::st_geometry(.))
  )
  
  if (nrow(all_dams_b) == 0L) next
  
  all_dams_b <- sf::st_as_sf(all_dams_b, sf_column_name = "geometry", crs = crs_r)
  net_with_dams_b <- sfnetworks::st_network_blend(net_b, all_dams_b, tolerance = blend_tolerance_m)
  
  for (t in thresholds_km) {
    t_label <- as.character(t)
    message("  threshold km: ", t_label)
    
    conn_b <- connectivity_from_network(
      net_with_dams = net_with_dams_b,
      threshold_downstream_km = t,
      threshold_upstream_km = t,
      threshold_cascade_km = t
    )
    
    reach_b <- conn_b$reach_df
    
    reach_enriched_b <- add_ffr_attr(
      net_with_dams = net_with_dams_b,
      reach_df = reach_b,
      edge_attr_cols = edge_attr_cols,
      future_extra = dam_extra_tbl
    ) %>%
      dplyr::mutate(main_riv = r, threshold_label = t_label)
    
    results_by_threshold[[t_label]] <- dplyr::bind_rows(results_by_threshold[[t_label]], reach_enriched_b)
  }
}

# =============================================================================
# Save output
# =============================================================================

out <- list(
  rivers_filtered = rivers_filtered,
  dams_current_gdw = current_dams_gdw,
  dams_future_fhred_clean = dams_future_fhred_clean,
  dam_extra_tbl = dam_extra_tbl,
  main_riv_dist = main_riv_dist,
  results_by_threshold = results_by_threshold
)

saveRDS(out, file = out_net_path)
message("Saved network/workflow bundle: ", out_net_path)

for (nm in names(results_by_threshold)) {
  df <- results_by_threshold[[nm]]
  out_rds <- file.path(output_dir, paste0("reach_enriched_threshold_", nm, "km.rds"))
  saveRDS(df, out_rds)
  message("Saved: ", out_rds, " (rows=", nrow(df), ")")
}

if (isTRUE(write_combined)) {
  combined <- dplyr::bind_rows(results_by_threshold)
  combined_path <- file.path(output_dir, "reach_enriched_all_thresholds.rds")
  saveRDS(combined, combined_path)
  message("Saved combined: ", combined_path, " (rows=", nrow(combined), ")")
}

message("Done.")
