#!/usr/bin/env Rscript

# =============================================================================
# world_connectivity_workflow.R
# =============================================================================
# PURPOSE
#   Run world-scale GDW + FHReD connectivity in main_riv partitions, write
#   threshold-specific enriched outputs, and save global plotting layers as sf.
#
# NOTES
#   - No global sfnetwork blend is built (expensive); connectivity is computed
#     per main_riv partition.
#   - Final plotting artifact is an sf bundle (global_edges_sf, global_dams_sf),
#     not sfnetwork.
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(sfnetworks)
  library(tidygraph)
  library(dplyr)
  library(readr)
  library(janitor)
  library(igraph)
  library(future.apply)
})

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

# Input paths
path_gdw_gdb <- "/capstone/netzerohydro/data/raw/GDW/GDW_v1_0.gdb"
layer_gdw <- "GDW_barriers_v1_0"
path_rivers_rds <- "/capstone/netzerohydro/data/cleaned_v1/world_rivers_ffr.rds"
path_fhred_csv <- "/capstone/netzerohydro/data/raw/FHReD_2015_future_dams/FHReD_2015_future_dams_Zarfl_et_al_beta_version.csv"

# Filters
main_riv_filter <- NULL
ord_stra_min <- 0L

# Distances (meters)
fhred_overlap_tolerance_m <- 500
blend_tolerance_m <- 500

# Thresholds
thresholds_km <- c(10, 20, 50, 100, 250)
include_inf_threshold <- FALSE
if (isTRUE(include_inf_threshold)) thresholds_km <- c(thresholds_km, Inf)

# Enrichment edge attrs
edge_attr_cols <- c("csi", "bas_name", "main_riv", "hyriv_id")

# Output
output_dir <- "/capstone/netzerohydro/data/processed_data/global"
out_net_path <- file.path(output_dir, "out_world_net_with_dams.rds")
plot_layers_path <- file.path(output_dir, "global_plot_layers.rds")
write_combined <- TRUE

# Runtime knobs
max_main_riv <- NULL
run_parallel <- FALSE
parallel_workers <- max(1L, parallel::detectCores() - 1L)
checkpoint_every <- 25L

# Functions
source("~/MEDS/capstone/connectivity/connectivity_function_v3.R")
source("~/MEDS/capstone/connectivity/add_ffr_attr.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# LOAD DATA
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
# FILTER RIVERS
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
# PREP DAM LAYERS (single canonical block)
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
# DEDUPE FHRED NEAR GDW
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
# BUILD DAM EXTRA TABLE USED BY add_ffr_attr()
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
# PARTITIONS
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
plot_edges_parts <- vector("list", n_part)
plot_dams_parts <- vector("list", n_part)

# counters
n_ok <- 0L
n_skip_small <- 0L
n_skip_no_dams <- 0L
n_skip_net <- 0L
n_skip_blend <- 0L

# =============================================================================
# PARTITION WORKER
# =============================================================================

process_one_main_riv <- function(r) {
  rivers_b <- rivers_filtered %>% dplyr::filter(.data$main_riv == r)
  if (nrow(rivers_b) < 2L) return(list(status = "skip_small"))
  
  net_b <- tryCatch(
    sfnetworks::as_sfnetwork(rivers_b, directed = TRUE) %>%
      tidygraph::activate("edges") %>%
      dplyr::mutate(weight = sfnetworks::edge_length()),
    error = function(e) NULL
  )
  if (is.null(net_b)) return(list(status = "skip_net"))
  
  idx_cur <- sf::st_nearest_feature(current_dams_gdw, rivers_b)
  idx_fut <- sf::st_nearest_feature(dams_future_fhred_clean, rivers_b)
  
  nearest_cur <- suppressWarnings(sf::st_distance(current_dams_gdw, rivers_b[idx_cur, ], by_element = TRUE))
  nearest_fut <- suppressWarnings(sf::st_distance(dams_future_fhred_clean, rivers_b[idx_fut, ], by_element = TRUE))
  
  dams_cur_b <- current_dams_gdw[as.numeric(nearest_cur) <= blend_tolerance_m, , drop = FALSE]
  dams_fut_b <- dams_future_fhred_clean[as.numeric(nearest_fut) <= blend_tolerance_m, , drop = FALSE]
  
  dams_fut_b_small <- dams_fut_b %>%
    dplyr::transmute(
      dam_id = as.character(.data$dam_id),
      is_current_dam = as.logical(.data$is_current_dam),
      geometry = sf::st_geometry(.)
    ) %>%
    sf::st_as_sf(sf_column_name = "geometry", crs = sf::st_crs(rivers_filtered))
  
  dams_cur_b_small <- dams_cur_b %>%
    dplyr::transmute(
      dam_id = as.character(.data$dam_id),
      is_current_dam = as.logical(.data$is_current_dam),
      geometry = sf::st_geometry(.)
    ) %>%
    sf::st_as_sf(sf_column_name = "geometry", crs = sf::st_crs(rivers_filtered))
  
  parts <- Filter(function(x) nrow(x) > 0, list(dams_fut_b_small, dams_cur_b_small))
  if (length(parts) == 0L) return(list(status = "skip_no_dams"))
  
  # guaranteed identical names/order
  parts <- lapply(parts, function(x) x[, c("dam_id", "is_current_dam", "geometry"), drop = FALSE])
  
  all_dams_b <- tryCatch(do.call(rbind, parts), error = function(e) NULL)
  if (is.null(all_dams_b)) return(list(status = "skip_blend"))
  
  net_with_dams_b <- tryCatch(
    sfnetworks::st_network_blend(net_b, all_dams_b, tolerance = blend_tolerance_m),
    error = function(e) NULL
  )
  if (is.null(net_with_dams_b)) return(list(status = "skip_blend"))
  
  threshold_out <- stats::setNames(vector("list", length(thresholds_km)), as.character(thresholds_km))
  
  for (t in thresholds_km) {
    conn_b <- connectivity_from_network(
      net_with_dams = net_with_dams_b,
      threshold_downstream_km = t,
      threshold_upstream_km = t,
      threshold_cascade_km = t
    )
    
    threshold_out[[as.character(t)]] <- add_ffr_attr(
      net_with_dams = net_with_dams_b,
      reach_df = conn_b$reach_df,
      edge_attr_cols = edge_attr_cols,
      future_extra = dam_extra_tbl
    ) %>% dplyr::mutate(main_riv = r, threshold_label = as.character(t))
  }
  
  edges_sf <- net_with_dams_b %>%
    tidygraph::activate("edges") %>%
    sf::st_as_sf() %>%
    dplyr::mutate(main_riv = r)
  
  dams_sf <- net_with_dams_b %>%
    tidygraph::activate("nodes") %>%
    sf::st_as_sf() %>%
    dplyr::filter(!is.na(.data$dam_id)) %>%
    dplyr::mutate(main_riv = r)
  
  list(
    status = "ok",
    threshold_out = threshold_out,
    edges_sf = edges_sf,
    dams_sf = dams_sf,
    rivers_rows = nrow(rivers_b),
    dam_rows = nrow(all_dams_b)
  )
}

# =============================================================================
# RUN PARTITIONS (serial or parallel)
# =============================================================================

if (isTRUE(run_parallel)) {
  log_msg("Running in parallel with workers=", parallel_workers)
  old_plan <- future::plan(future::multisession, workers = parallel_workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  part_results <- future.apply::future_lapply(main_riv_values, process_one_main_riv, future.seed = TRUE)
  
  for (i in seq_along(part_results)) {
    r <- main_riv_values[[i]]
    res <- part_results[[i]]
    
    log_msg(sprintf("[main_riv %d/%d] %s status=%s", i, n_part, as.character(r), res$status))
    
    if (identical(res$status, "skip_small")) { n_skip_small <- n_skip_small + 1L; next }
    if (identical(res$status, "skip_no_dams")) { n_skip_no_dams <- n_skip_no_dams + 1L; next }
    if (identical(res$status, "skip_net")) { n_skip_net <- n_skip_net + 1L; next }
    if (identical(res$status, "skip_blend")) { n_skip_blend <- n_skip_blend + 1L; next }
    
    n_ok <- n_ok + 1L
    plot_edges_parts[[i]] <- res$edges_sf
    plot_dams_parts[[i]] <- res$dams_sf
    
    for (t in as.character(thresholds_km)) {
      df_t <- res$threshold_out[[t]]
      log_msg("  threshold ", t, " rows=", nrow(df_t))
      results_by_threshold[[t]] <- dplyr::bind_rows(results_by_threshold[[t]], df_t)
    }
    
    if (i %% checkpoint_every == 0L) {
      log_msg(sprintf("[checkpoint %d/%d] ok=%d small=%d no_dams=%d net_fail=%d blend_fail=%d", i, n_part, n_ok, n_skip_small, n_skip_no_dams, n_skip_net, n_skip_blend))
    }
  }
  
} else {
  log_msg("Running in serial mode")
  
  for (i in seq_along(main_riv_values)) {
    r <- main_riv_values[[i]]
    log_msg(sprintf("[main_riv %d/%d] %s", i, n_part, as.character(r)))
    
    res <- process_one_main_riv(r)
    
    if (identical(res$status, "skip_small")) { n_skip_small <- n_skip_small + 1L; log_msg("  skip: too few river rows"); next }
    if (identical(res$status, "skip_no_dams")) { n_skip_no_dams <- n_skip_no_dams + 1L; log_msg("  skip: no dams within tolerance"); next }
    if (identical(res$status, "skip_net")) { n_skip_net <- n_skip_net + 1L; log_msg("  skip: as_sfnetwork failed"); next }
    if (identical(res$status, "skip_blend")) { n_skip_blend <- n_skip_blend + 1L; log_msg("  skip: st_network_blend failed"); next }
    
    n_ok <- n_ok + 1L
    log_msg("  rivers rows: ", res$rivers_rows, " | blended dams rows: ", res$dam_rows)
    
    plot_edges_parts[[i]] <- res$edges_sf
    plot_dams_parts[[i]] <- res$dams_sf
    
    for (t in as.character(thresholds_km)) {
      df_t <- res$threshold_out[[t]]
      log_msg("  threshold ", t, " rows=", nrow(df_t))
      results_by_threshold[[t]] <- dplyr::bind_rows(results_by_threshold[[t]], df_t)
    }
    
    if (i %% checkpoint_every == 0L) {
      log_msg(sprintf("[checkpoint %d/%d] ok=%d small=%d no_dams=%d net_fail=%d blend_fail=%d", i, n_part, n_ok, n_skip_small, n_skip_no_dams, n_skip_net, n_skip_blend))
    }
  }
}

log_msg("Loop summary: ok=", n_ok,
        " | skip_small=", n_skip_small,
        " | skip_no_dams=", n_skip_no_dams,
        " | skip_net=", n_skip_net,
        " | skip_blend=", n_skip_blend)

# =============================================================================
# GLOBAL PLOT LAYERS (non-sfnetwork)
# =============================================================================

global_edges_sf <- dplyr::bind_rows(Filter(Negate(is.null), plot_edges_parts))
global_dams_sf <- dplyr::bind_rows(Filter(Negate(is.null), plot_dams_parts))

global_plot_layers <- list(
  global_edges_sf = global_edges_sf,
  global_dams_sf = global_dams_sf
)

saveRDS(global_plot_layers, plot_layers_path)
log_msg("Saved global plot layers: ", plot_layers_path,
        " | edges=", nrow(global_edges_sf),
        " | dams=", nrow(global_dams_sf))

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

log_msg("Saving outputs...")

out <- list(
  rivers_filtered = rivers_filtered,
  dams_current_gdw = current_dams_gdw,
  dams_future_fhred_clean = dams_future_fhred_clean,
  dam_extra_tbl = dam_extra_tbl,
  main_riv_dist = main_riv_dist,
  results_by_threshold = results_by_threshold,
  plot_layers_path = plot_layers_path
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
