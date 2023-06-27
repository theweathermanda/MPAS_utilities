# MPAS_utilities

This repository contains various utilities for working with the native MPAS grid mesh.

[**plot_mpas_zoom_edges.ncl**](https://github.com/theweathermanda/MPAS_utilities/blob/main/plot_mpas_zoom_edges.ncl) : This script uses CellFill and gsSegments to plot cell and edge fields from the native MPAS mesh. The edge segments are colored based on the value of the edge fields. 
   - _In this specific example, the edge field plotted is the magnitude of the edge-normal theta gradient, which was precomputed and saved as a netcdf file. The edge values are normalized by the maximum value, and colors are assigned to the normalized edge values from a specified colormap. Here, black values correspond to large gradient magnitudes, and white values correspond to small gradient magnitudes._

![alt text](https://github.com/theweathermanda/MPAS_utilities/blob/main/theta_gradient_edge_magnitude_lev7.png?raw=true)

[**mpas_calc_operators.py**](https://github.com/theweathermanda/MPAS_utilities/blob/c3b35e42ea682a121fdd3a1ffd29d3ef59068666/mpas_calc_operators.py) : This script contains the following functions that operate on the native MPAS mesh:
- mpas_calc_edgesOnCell_sign: _populate edgesOnCell_sign variable if not contained within outfiles_
- mpas_calc_gradOnEdges: _calculate horizontal gradient of a field on each cell edge (e.g., potential temperature gradient)_
-- mpas_reconstruct_grad: _reconstruct horizontal edge gradients to the cell center (e.g., zonal potential temperature gradient)_
-- mpas_calc_curl: _calculate finite volume curl at the cell vertices using a vector on the cell edges (e.g., vertical vorticity)_
