# Assigning River Network Connectivity Status for Future Dams Based on Spatial Relationship to Exisitng Dams

## Network Connenctivity analysis

**Authors** Leela Dixit, Megan Hessel, Aakriti Poudel, and Lucian Scher

### Purpose

This repository contains three functions and tests to build a sfnetworks river network with dams, calculate connecvitiy and assign additonal attributes to be used in the MCDA model. 

### Data

The datasets used for current and future dams, river network geometries and additional attributes can be accessed for free online using the links below. 

-   [Global Dam Watch (GDW)](https://www.globaldamwatch.org/):Database that has all existing hydropwer projects.

-   [Future Hydropower and Reservoir Data (FhRED)](https://www.globaldamwatch.org/directory):Data set with planned and future hydropower projects with capacity of at least 1 MW.

-   [Global River Networks (HydroRIVERS Dataset)](https://www.hydrosheds.org/products/hydrorivers): Provides vectorized spatial global river network data at 15 arc-second resolution, approximately 500 m at the equator, and includes river attributes.

-   [Free-Flowing Rivers database (FFR)](https://www.hydrosheds.org/applications/free-flowing-rivers):Using HydroSHED river network, creates a underpinning hydrographic data to support the identification of free-flowing and at-risk rivers.

## File Structure

IN PROGRESS

### R 

DESCRIBE FUNCTIONS

### Packages used

Analysis was done in R v4.5.2 using the following packages:

- sf: Handles all spatial vector data workflows (points/lines, CRS transforms, geometry operations, nearest-feature matching, and spatial filtering).
- sfnetworks: Builds directed river sfnetwork objects and blends snapped dam points into the network topology.
- tidygraph: Switches between network node/edge tables and supports tidy manipulation of graph attributes.
- igraph: Computes directed weighted path distances and trunk-hop distances used in upstream/downstream and cascade logic.
- dplyr: Powers nearly all data wrangling (filtering, joining, grouping, summarizing, mutating, and selecting fields).

### Technical Documentation

To read more about the project and modeling processes, please refer to our [Bren project page](https://bren.ucsb.edu/projects/hydropowers-low-hanging-fruits-leveraging-least-impact-dams-power-net-zero-futurehttps://bren.ucsb.edu/projects/hydropowers-low-hanging-fruits-leveraging-least-impact-dams-power-net-zero-future) and technical documentation.
