# MPAS_utilities

This repository contains utilities for working with the native MPAS grid mesh.

plot_mpas_zoom_edges.ncl : This script uses CellFill and gsSegments to plot cell and edge fields from the native MPAS mesh. The edge segments are colored based on the value of the edge fields. In this specific example, the edge field plotted is the magnitude of the edge-normal theta gradient, which was precomputed and saved as a netcdf file.
