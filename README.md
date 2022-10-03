# MPAS_utilities

This repository contains utilities for working with the native MPAS grid mesh.

[**plot_mpas_zoom_edges.ncl**](https://github.com/theweathermanda/MPAS_utilities/blob/main/plot_mpas_zoom_edges.ncl) : This script uses CellFill and gsSegments to plot cell and edge fields from the native MPAS mesh. The edge segments are colored based on the value of the edge fields. 
   - _In this specific example, the edge field plotted is the magnitude of the edge-normal theta gradient, which was precomputed and saved as a netcdf file. The edge values are normalized by the maximum value, and colors are assigned to the normalized edge values from a specified colormap. Here, black values correspond to large gradient magnitudes, and white values correspond to small gradient magnitudes._

![alt text](https://github.com/theweathermanda/MPAS_utilities/blob/main/theta_gradient_edge_magnitude_lev7.png?raw=true)
