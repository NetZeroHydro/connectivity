install.packages("sfnetworks")
install.packages("dplyr")
install.packages("tidygraph")
install.packages("sf")
install.packages("igraph")

library(sfnetworks)
library(dplyr)
library(tidygraph)
library(sf)
library(igraph)

out <- readRDS(file = "home/lucian/hpc/netzerohydro/out.RDS")
source("home/lucian/hpc/connectivity_from_network.R")

thresholds <- c("10", "20", "50", "100", "250", "Inf")


for(threshold in thresholds){
  out_conn <- connectivity_from_network(
    out$net_with_dams,
    threshold_upstream_km = threshold,
    threshold_downstream_km = threshold,
    threshold_cascade_km = threshold
  )
  
  saveRDS(out_conn, file = file.path("home/lucian/hpc/", paste0("out_conn_", threshold, ".rds")))
}
