# Assigning River Network Connectivity Status for Future Dams Based on Spatial Relationship to Existing Dams

**Authors:** Leela Dixit, Megan Hessel, Aakriti Poudel, and Lucian Scher

## Purpose

This repository implements a **connectivity classification workflow** for future hydropower dams on a directed river network. Given a blended network of existing (“current”) and planned (“future”) dam points snapped to river reaches, it assigns each future dam a **connectivity category** (e.g. cascade, upstream-only, downstream-only, undammed) based on directed river distance and trunk-hop relationships to the nearest current dams.

The outputs are intended for downstream multi-criteria analysis: one row per dam with neighbor IDs, distances, hop counts, cascade level, and optional river/dam attributes. This repo does **not** run the full MCDA model—it supplies the connectivity layer and enrichment steps only.

## Repository layout

```
├── 0_setup.qmd                      # Packages, paths, and raw data loading (run first)
├── 1_net_with_dams_from_network.qmd # Build directed sfnetwork with blended dams
├── 2_run_connectivity_from_network.qmd  # Run connectivity classification
├── 3_enrich_out_conn.qmd            # Append FFR / river attributes to reach_df
├── alt_dam_data.qmd                 # Minimal guide: use connectivity on other data
├── LICENSE
├── R/
│   ├── add_ffr_attr.R               # Enrich reach_df with edge/dam columns
│   └── connectivity_from_network.R  # Core connectivity logic
└── README.md
```

| File | Role |
|------|------|
| `0_setup.qmd` | Defines editable `paths`, loads libraries, reads HydroRIVERS, FFR, GDW, and FHReD into the R session. |
| `1_net_with_dams_from_network.qmd` | Filters rivers, de-duplicates FHReD vs GDW, blends dams into an `sfnetwork`, produces `out`. |
| `2_run_connectivity_from_network.qmd` | Calls `connectivity_from_network()`; produces `out_conn` with `reach_df`. |
| `3_enrich_out_conn.qmd` | Calls `add_ffr_attr()`; produces `out_enriched`. |
| `alt_dam_data.qmd` | Short checklist for applying the functions to a custom network. |

## R functions

### `connectivity_from_network()`

**File:** `R/connectivity_from_network.R`

Takes a directed `sfnetwork` with dams already blended onto the graph (`net_with_dams`). For each **future** dam node, it finds the nearest **current** dam upstream and downstream along the river network (within km thresholds), using trunk-hop count then river km as tie-breakers. It then assigns `cascade_level` and `connectivity_category` (e.g. `cascade_classic`, `cascade1`, `cascade2+`, `upstream`, `downstream`, `undammed`).

**Produces:** A named list with:

- `reach_df` — one row per dam (future + current placeholders) with distances, neighbor dam IDs, hop counts, and category
- `debug` — intermediate tables for QA
- `threshold_used` — km thresholds applied

### `add_ffr_attr()`

**File:** `R/add_ffr_attr.R`

Appends columns to `reach_df` without recomputing connectivity. Joins river-edge attributes from `net_with_dams` by `bb_id`, and optionally extra dam-level columns by `dam_id`.

**Produces:** The same `reach_df` structure with additional attribute columns (e.g. `csi`, `bas_name`, `hyriv_id`).

## Workflow

Run the Quarto notebooks **in order in the same R session** (Render All, or run each interactively without restarting R):

1. `0_setup.qmd` — edit `paths` if needed; loads packages and raw data
2. `1_net_with_dams_from_network.qmd` — build network (`out`)
3. `2_run_connectivity_from_network.qmd` — connectivity (`out_conn`)
4. `3_enrich_out_conn.qmd` — enrichment (`out_enriched`)

Each notebook lists prerequisites at the top. To use different raw files, change the `paths` list in `0_setup.qmd` (default data root: `/capstone/netzerohydro/data`).

To run connectivity on **your own** river network and dam points, see [`alt_dam_data.qmd`](alt_dam_data.qmd).

## Data sources

Datasets used in the default pipeline (free online):

- [Global Dam Watch (GDW)](https://www.globaldamwatch.org/) — existing hydropower projects
- [Future Hydropower and Reservoir Data (FhRED)](https://www.globaldamwatch.org/directory) — planned projects ≥ 1 MW
- [HydroRIVERS](https://www.hydrosheds.org/products/hydrorivers) — global river network (~500 m)
- [Free-Flowing Rivers (FFR)](https://figshare.com/articles/dataset/Mapping_the_world_s_free-flowing_rivers_data_set_and_technical_documentation/7688801) — reach attributes on HydroSHEDS network

## Reproducibility

### Session info

Analysis was run with:

| Component | Version | Role |
|-----------|---------|------|
| R | 4.5.2 | — |
| sf | 1.0.14 | Spatial I/O, CRS, geometry operations |
| sfnetworks | 0.6.5 | Directed river graph and dam blending |
| tidygraph | 1.3.1 | Node/edge table access on the graph |
| igraph | 2.2.1 | Shortest-path distances and trunk-hop matrix |
| dplyr | 1.2.1 | Data wrangling |
| janitor | 2.2.0 | `clean_names()` on imported tables |
| readr | 2.1.4 | Tabular I/O where used |
| leaflet | 2.2.3 | Optional maps in notebook 1 |

Check versions on your machine:

```r
pkgs <- c("sf", "sfnetworks", "tidygraph", "igraph", "dplyr", "janitor", "readr", "leaflet")
sapply(pkgs, packageVersion)
```

### How to recreate the analysis

<img width="510" height="334" alt="image" src="https://github.com/user-attachments/assets/e46caeaf-61f5-444b-89e2-8a52db614212" />

1. **Install R** (4.5.2 or compatible) and the packages above (`install.packages()` for CRAN packages; `sf` may need system GDAL/GEOS libraries).

2. **Obtain data** from the links in [Data sources](#data-sources) and place files on disk (or use the Bren capstone paths if available).

3. **Clone this repo** and open the project in RStudio or VS Code with Quarto.

4. **Edit paths** in `0_setup.qmd` — set `data_root` and individual file paths to match your machine.

5. **Run notebooks 0 → 3 in one R session** without restarting R between steps. Objects (`world_rivers_ffr`, `out`, `out_conn`, etc.) persist in the global environment.

6. **Optional saves:** notebook 1 writes `out` to `paths$net_with_dams_rds`; notebook 2 writes `out_conn` under `paths$processed`. Adjust filenames in those notebooks if needed.

7. **Custom data:** if you are not using GDW/FHReD/HydroRIVERS, build your own `net_with_dams` following [`alt_dam_data.qmd`](alt_dam_data.qmd), then call `connectivity_from_network()` directly.

## Technical documentation

For the broader Bren project and MCDA context, see the [Bren project page](https://bren.ucsb.edu/projects/hydropowers-low-hanging-fruits-leveraging-least-impact-dams-power-net-zero-future) and technical documentation.
