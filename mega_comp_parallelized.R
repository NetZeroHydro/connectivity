# install.packages("sfnetworks")
# install.packages("sf")
# install.packages("tidygraph")
# install.packages("dplyr")
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
# LOGGING HELPERS
# =============================================================================
t0 <- Sys.time()
log_msg <- function(...) {
  elapsed_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  message(sprintf("[+%s min] ", elapsed_min), paste0(..., collapse = ""))
}

log_msg("=== START world workflow ===")
log_msg("R version: ", R.version.string)

# =============================================================================
# CONFIG
# =============================================================================
path_gdw_gdb <- "/capstone/netzerohydro/data/raw/GDW/GDW_v1_0.gdb"
layer_gdw <- "GDW_barriers_v1_0"
path_rivers_rds <- "/capstone/netzerohydro/data/cleaned_v1/world_rivers_ffr.rds"
path_fhred_csv <- "/capstone/netzerohydro/data/raw/FHReD_2015_future_dams/FHReD_2015_future_dams_Zarfl_et_al_beta_version.csv"

main_riv_filter <- NULL
ord_stra_min <- 0L

fhred_overlap_tolerance_m <- 500
blend_tolerance_m <- 500

thresholds_km <- c(10, 20, 50, 100, 250, Inf)
edge_attr_cols <- c("csi", "bas_name", "main_riv", "hyriv_id")

output_dir <- "/capstone/netzerohydro/data/processed_data"
out_net_path <- file.path(output_dir, "out_world_net_with_dams.rds")
write_combined <- TRUE
max_main_riv <- NULL

source("connectivity_function_v3.R")
source("add_ffr_attr.R")

# =============================================================================
# Load Data
# =============================================================================
log_msg("[1/8] Loading GDW...")
gdw <- sf::st_read(path_gdw_gdb, layer = layer_gdw, quiet = TRUE) |>
  janitor::clean_names()
log_msg("GDW rows: ", nrow(gdw), " cols: ", ncol(gdw))

log_msg("[2/8] Loading world_rivers_ffr...")
world_rivers_ffr <- readRDS(path_rivers_rds)
log_msg("world_rivers_ffr rows: ", nrow(world_rivers_ffr), " cols: ", ncol(world_rivers_ffr))

log_msg("[3/8] Loading FHReD...")
future_dams_raw <- readr::read_csv(path_fhred_csv, show_col_types = FALSE) |>
  janitor::clean_names()
log_msg("FHReD rows: ", nrow(future_dams_raw), " cols: ", ncol(future_dams_raw))

# =============================================================================
# Filter rivers
# =============================================================================
log_msg("[4/8] Filtering rivers...")
rivers_filtered <- world_rivers_ffr %>%
  dplyr::filter(.data$ord_stra > ord_stra_min)

if (!is.null(main_riv_filter)) {
  rivers_filtered <- rivers_filtered %>%
    dplyr::filter(as.numeric(.data$main_riv) == as.numeric(main_riv_filter))
  log_msg("Applied main_riv_filter = ", main_riv_filter)
}

rivers_filtered <- suppressWarnings(sf::st_cast(rivers_filtered, "LINESTRING"))
rivers_filtered <- rivers_filtered[!sf::st_is_empty(rivers_filtered), , drop = FALSE]

if (nrow(rivers_filtered) == 0L) {
  stop("No rivers left after filters. Check ord_stra/main_riv settings.", call. = FALSE)
}
if (!"main_riv" %in% names(rivers_filtered)) {
  stop("`main_riv` is required on rivers for world-scale partitioning.", call. = FALSE)
}
log_msg("rivers_filtered rows: ", nrow(rivers_filtered))

main_riv_dist <- rivers_filtered %>%
  dplyr::mutate(main_riv = as.character(.data$main_riv)) %>%
  dplyr::mutate(seg_len_m = as.numeric(sf::st_length(sf::st_geometry(.)))) %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(.data$main_riv) %>%
  dplyr::summarise(main_riv_dist = sum(.data$seg_len_m, na.rm = TRUE) / 1000, .groups = "drop")
log_msg("main_riv_dist groups: ", nrow(main_riv_dist))

# =============================================================================
# Clean Future and Current Dams data
# =============================================================================
log_msg("[5/8] Preparing dam layers...")
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

log_msg("current_dams_gdw rows: ", nrow(current_dams_gdw))
log_msg("future_dams_sf rows: ", nrow(future_dams_sf))

# =============================================================================
# Handle duplicate Future and current dams
# =============================================================================
log_msg("[6/8] De-duplicating FHReD near GDW...")
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
log_msg("FHReD removed by overlap: ", dams_removed, " | remaining: ", nrow(dams_future_fhred_clean))

# =============================================================================
# Build dam extra table used by add_ffr_attr()
# =============================================================================
log_msg("[7/8] Building dam_extra_tbl...")
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
log_msg("dam_extra_tbl rows: ", nrow(dam_extra_tbl))

# =============================================================================
# Optional global plotting network (very expensive)
# =============================================================================
log_msg("[optional] Building global plotting network...")
all_dams_plot <- dplyr::bind_rows(
  dams_future_fhred_clean %>%
    dplyr::transmute(dam_id = .data$dam_id, is_current_dam = .data$is_current_dam, geometry = sf::st_geometry(.)),
  current_dams_gdw %>%
    dplyr::transmute(dam_id = .data$dam_id, is_current_dam = .data$is_current_dam, geometry = sf::st_geometry(.))
) %>%
  sf::st_as_sf(sf_column_name = "geometry", crs = sf::st_crs(rivers_filtered))

net_plot_world <- sfnetworks::as_sfnetwork(rivers_filtered, directed = TRUE) %>%
  tidygraph::activate("edges") %>%
  dplyr::mutate(weight = sfnetworks::edge_length()) %>%
  sfnetworks::st_network_blend(all_dams_plot, tolerance = blend_tolerance_m)

saveRDS(net_plot_world, file.path(output_dir, "net_with_dams_plotting_world.rds"))
log_msg("Saved net_with_dams_plotting_world.rds")

# =============================================================================
# Iterate by main_riv + thresholds
# =============================================================================
main_riv_values <- rivers_filtered %>%
  sf::st_drop_geometry() %>%
  dplyr::distinct(.data$main_riv) %>%
  dplyr::pull(.data$main_riv)

main_riv_values <- main_riv_values[!is.na(main_riv_values)]
if (!is.null(max_main_riv)) main_riv_values <- head(main_riv_values, max_main_riv)

n_part <- length(main_riv_values)
log_msg("[8/8] Total main_riv partitions: ", n_part)

results_by_threshold <- stats::setNames(vector("list", length(thresholds_km)), as.character(thresholds_km))

# counters
n_ok <- 0L
n_skip_small <- 0L
n_skip_net <- 0L
n_skip_no_dams <- 0L
n_skip_blend <- 0L

for (i in seq_along(main_riv_values)) {
  r <- main_riv_values[i]
  log_msg(sprintf("[main_riv %d/%d] %s", i, n_part, as.character(r)))
  
  rivers_b <- rivers_filtered %>% dplyr::filter(.data$main_riv == r)
  log_msg("  rivers_b rows: ", nrow(rivers_b))
  if (nrow(rivers_b) == 0L) next
  
  if (nrow(rivers_b) < 2L) {
    n_skip_small <- n_skip_small + 1L
    log_msg("  skip: too few river rows")
    next
  }
  
  net_b <- tryCatch(
    sfnetworks::as_sfnetwork(rivers_b, directed = TRUE) %>%
      tidygraph::activate("edges") %>%
      dplyr::mutate(weight = sfnetworks::edge_length()),
    error = function(e) {
      n_skip_net <<- n_skip_net + 1L
      log_msg("  skip: as_sfnetwork failed -> ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(net_b)) next
  
  idx_cur <- sf::st_nearest_feature(current_dams_gdw, rivers_b)
  idx_fut <- sf::st_nearest_feature(dams_future_fhred_clean, rivers_b)
  
  nearest_cur <- suppressWarnings(sf::st_distance(current_dams_gdw, rivers_b[idx_cur, ], by_element = TRUE))
  nearest_fut <- suppressWarnings(sf::st_distance(dams_future_fhred_clean, rivers_b[idx_fut, ], by_element = TRUE))
  
  dams_cur_b <- current_dams_gdw[as.numeric(nearest_cur) <= blend_tolerance_m, , drop = FALSE]
  dams_fut_b <- dams_future_fhred_clean[as.numeric(nearest_fut) <= blend_tolerance_m, , drop = FALSE]
  log_msg("  dams_cur_b rows: ", nrow(dams_cur_b), " | dams_fut_b rows: ", nrow(dams_fut_b))
  
  dams_fut_b_small <- dams_fut_b %>%
    dplyr::transmute(dam_id = as.character(.data$dam_id), is_current_dam = .data$is_current_dam) %>%
    sf::st_as_sf()
  
  dams_cur_b_small <- dams_cur_b %>%
    dplyr::transmute(dam_id = as.character(.data$dam_id), is_current_dam = .data$is_current_dam) %>%
    sf::st_as_sf()
  
  parts <- Filter(function(x) nrow(x) > 0, list(dams_fut_b_small, dams_cur_b_small))
  if (length(parts) == 0L) {
    n_skip_no_dams <- n_skip_no_dams + 1L
    log_msg("  skip: no dams within tolerance")
    next
  }
  
  all_dams_b <- do.call(sf::rbind, parts)
  
  net_with_dams_b <- tryCatch(
    sfnetworks::st_network_blend(net_b, all_dams_b, tolerance = blend_tolerance_m),
    error = function(e) {
      n_skip_blend <<- n_skip_blend + 1L
      log_msg("  skip: st_network_blend failed -> ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(net_with_dams_b)) next
  
  for (t in thresholds_km) {
    t_label <- as.character(t)
    log_msg("    threshold: ", t_label, " km")
    
    conn_b <- connectivity_from_network(
      net_with_dams = net_with_dams_b,
      threshold_downstream_km = t,
      threshold_upstream_km = t,
      threshold_cascade_km = t
    )
    
    reach_b <- conn_b$reach_df
    log_msg("    reach_b rows: ", nrow(reach_b))
    
    reach_enriched_b <- add_ffr_attr(
      net_with_dams = net_with_dams_b,
      reach_df = reach_b,
      edge_attr_cols = edge_attr_cols,
      future_extra = dam_extra_tbl
    ) %>%
      dplyr::mutate(main_riv = r, threshold_label = t_label)
    
    log_msg("    reach_enriched_b rows: ", nrow(reach_enriched_b))
    
    results_by_threshold[[t_label]] <- dplyr::bind_rows(results_by_threshold[[t_label]], reach_enriched_b)
  }
  
  n_ok <- n_ok + 1L
  
  if (i %% 25 == 0L) {
    log_msg(sprintf(
      "[checkpoint %d/%d] ok=%d small=%d no_dams=%d net_fail=%d blend_fail=%d",
      i, n_part, n_ok, n_skip_small, n_skip_no_dams, n_skip_net, n_skip_blend
    ))
  }
}

log_msg("Loop summary: ok=", n_ok,
        " | skip_small=", n_skip_small,
        " | skip_no_dams=", n_skip_no_dams,
        " | skip_net=", n_skip_net,
        " | skip_blend=", n_skip_blend)

# =============================================================================
# Save output
# =============================================================================
log_msg("Saving outputs...")

out <- list(
  rivers_filtered = rivers_filtered,
  dams_current_gdw = current_dams_gdw,
  dams_future_fhred_clean = dams_future_fhred_clean,
  dam_extra_tbl = dam_extra_tbl,
  main_riv_dist = main_riv_dist,
  results_by_threshold = results_by_threshold
)

saveRDS(out, file = out_net_path)
log_msg("Saved bundle: ", out_net_path)

for (nm in names(results_by_threshold)) {
  df <- results_by_threshold[[nm]]
  out_rds <- file.path(output_dir, paste0("reach_enriched_threshold_", nm, "km.rds"))
  saveRDS(df, out_rds)
  log_msg("Saved threshold ", nm, " -> rows=", nrow(df), " file=", out_rds)
}

if (isTRUE(write_combined)) {
  combined <- dplyr::bind_rows(results_by_threshold)
  combined_path <- file.path(output_dir, "reach_enriched_all_thresholds.rds")
  saveRDS(combined, combined_path)
  log_msg("Saved combined -> rows=", nrow(combined), " file=", combined_path)
}

log_msg("=== DONE ===")