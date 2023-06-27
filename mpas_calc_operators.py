'''
Script to compute the following on the native MPAS mesh:

-- mpas_calc_edgesOnCell_sign: populate edgesOnCell_sign variable if not contained within outfiles
-- mpas_calc_gradOnEdges: calculate horizontal gradient of a field on each cell edge (e.g., temperature gradient)
-- mpas_reconstruct_grad: reconstruct horizontal edge gradients to the cell center (e.g., zonal temperature gradient)
-- mpas_calc_curl: calculate finite volume curl at the cell vertices using a vector on the cell edges (e.g., vorticity)

Created by Manda Chasteen, NCAR/MMM
September 26, 2022

'''

import numpy as np
import xarray as xr

####

def mpas_calc_edgesOnCell_sign(edgesOnCell, cellsOnEdge):
    '''
    Given MPAS mesh information, calculate edgesOnCell_sign if it's not contained in the model output files. 
    
    Sign follows convention that "Positive u (normal) velocity is always defined as flow from cellsOnEdge[jEdge,0]
    to cellsOnEdge[jEdge,1] for edge jEdge" (MPAS tutorial 2019, with adaptations for Python's zero-based indexing).
    
    This variable is needed to compute gradients on the native MPAS mesh. 
    
    Parameters
    ----------
    edgesOnCell : 2D DataArray or NumPy array of dims [nCells, maxEdges], required 
        IDs of edges forming the boundary of a cell
    cellsOnEdge : 2D DataArray or NumPy of dims [nEdges, TWO], required
        IDs of cells divided by an edge
    
    Returns
    -------
    edgesOnCell_sign : 2D NumPy array of dims [nCells, maxEdges]
        Sign for edges surrounding a cell: convention is positive for positive outward normal velocity
    
    '''
    maxEdges = 10
    
    # Need number of cells. 
    if isinstance(edgesOnCell, xr.core.dataarray.DataArray):
        nCells = edgesOnCell.nCells
    elif isinstance(edgesOnCell, np.ndarray):
        nCells = np.arange(0,np.shape(edgesOnCell)[0],1)
    else:
         raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')

    # Need number of edges
    if isinstance(cellsOnEdge, xr.core.dataarray.DataArray):
        nEdges = cellsOnEdge.nEdges
    elif isinstance(cellsOnEdge, np.ndarray):
        nEdges = np.arange(0,np.shape(cellsOnEdge)[0],1)
    else:
         raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Initialize array for edgesOnCell_sign
    edgesOnCell_sign = np.zeros((len(nCells),maxEdges))
    
    # Loop over edges
    for j in np.arange(0,len(nEdges),1):
        # Add 1 to correct for 1-based indexing in Fortran and MPAS
        jEdge = j + 1 
        
        # Finds the correct index corresponding to each edge value for the cells on the edge
        edges_0 = edgesOnCell[cellsOnEdge[j,0].values-1].values
        edges_1 = edgesOnCell[cellsOnEdge[j,1].values-1].values

        index_j0 = np.where(edges_0 == jEdge)[0][0]
        index_j1 = np.where(edges_1 == jEdge)[0][0]
    
        edgesOnCell_sign[cellsOnEdge[j,0].values-1,index_j0] = 1
        edgesOnCell_sign[cellsOnEdge[j,1].values-1,index_j1] = -1
        
    return edgesOnCell_sign
    

def mpas_calc_gradOnEdges(cellVar, nEdgesOnCell, edgesOnCell, cellsOnEdge, edgesOnCell_sign, dcEdge):
    '''
    Given MPAS mesh information, calculate the horizontal gradient of a field on each cell edge using field values at the 
    adjacent primal cell centers as:
     
    varGrad[edgeInd,lev] = (cellVar[cellsOnEdge[edgeInd,1],k]-cellVar[cellsOnEdge[edgeInd,0],k]) / dcEdge[edgeInd]
    
    and then assign the correct sign based on its direction (i.e., into or out of the cell), the value of edgesOnCell_sign, and 
    the convention for the u winds: "Positive u (normal) velocity is always defined as flow from cellsOnEdge[jEdge,0]
    to cellsOnEdge[jEdge,1] for edge jEdge" (MPAS tutorial 2019, with adaptations for Python's zero-based indexing).

    The expression for calculating the gradient on each edge comes from Eq. 22 in Ringler et al. (2010), J. Comput. Phys.

    Parameters
    ----------
    cellVar : 1D or 2D DataArray or NumPy array of [nCells] or [nCells, nVertLevels]
        Variable at primal cell centers
    nEdgesOnCell : 1D or 2D DataArray or NumPy array of [nCells]
        Number of edges forming the boundary of a cell
    edgesOnCell : 2D DataArray or NumPy array of dims [nCells, maxEdges], required 
        IDs of edges forming the boundary of a cell
    cellsOnEdge : 2D DataArray or NumPy of dims [nEdges, TWO], required
        IDs of cells divided by an edge
    edgesOnCell_sign : 2D DataArray or NumPy array of dims [nCells, maxEdges]
        Sign for edges surrounding a cell: convention is positive for positive outward normal velocity
    dcEdge : 1D DataArray or NumPy array of dims [nEdges]
        Spherical distance between cells separated by an edge

    Returns
    -------
    varGrad : 1D or 2D DataArray or NumPy array of [nEdges] or [nEdges, nVertLevels]
        Gradient of variable on cell edges 
    
    '''
    maxEdges = 10
    
    # Check number of dimensions of cellVar
    if cellVar.ndim > 2:
        raise AttributeError('Variable has more than 2 dimensions. Exiting.')
    else:
        pass
    
    # Need number of cells. 
    if isinstance(edgesOnCell, xr.core.dataarray.DataArray):
        nCells = edgesOnCell.nCells
    elif isinstance(edgesOnCell, np.ndarray):
        nCells = np.arange(0,np.shape(edgesOnCell)[0],1)
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')

    # Need number of edges
    if isinstance(cellsOnEdge, xr.core.dataarray.DataArray):
        nEdges = cellsOnEdge.nEdges
    elif isinstance(cellsOnEdge, np.ndarray):
        nEdges = np.arange(0,np.shape(cellsOnEdge)[0],1)
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Need number of vertical levels
    if isinstance(cellVar, xr.core.dataarray.DataArray):
        if cellVar.ndim == 2:
            nLevs = cellVar.nVertLevels
        else:
            nLevs = [1]
    elif isinstance(cellVar, np.ndarray):
        if cellVar.ndim == 2:
            nLevs = np.arange(0,np.shape(cellVar)[1],1)
        else:
            nLevs = [1]
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Initialize array for varGrad
    try:
        varGrad = np.zeros((len(nEdges),len(nLevs)))
    except:
        varGrad = np.zeros((len(nEdges)))
    
    # Loop over cells, edges on cell, and vertical levels
    for i in np.arange(0,len(nCells),1):
        for j in np.arange(0,nEdgesOnCell[i].values,1):
            
            # Edge IDs, inds, and signs for j edge along parent i cell
            edgeUse = edgesOnCell[i,j]
            edgeInd = edgeUse - 1
            edgeSign = edgesOnCell_sign[i,j]

            # The indices of edgeUse likely differ in edgesOnCell array for each cell.
            # Need to find the correct indices and the sign of the normal vector for each edge in edgesOnCell
            # -- if sign_j0 > 0, normal vector points out of cellsOnEdge[edgeInd,0] 
            # -- if sign_j1 > 0, normal vector points out of cellsOnEdge[edgeInd,1] 
            
            edges_0 = edgesOnCell[cellsOnEdge[edgeInd,0].values-1].values
            edges_1 = edgesOnCell[cellsOnEdge[edgeInd,1].values-1].values

            index_j0 = np.where(edges_0 == edgeUse.values)[0][0]
            index_j1 = np.where(edges_1 == edgeUse.values)[0][0]
    
            sign_j0 = edgesOnCell_sign[cellsOnEdge[edgeInd,0].values-1, index_j0].values
            sign_j1 = edgesOnCell_sign[cellsOnEdge[edgeInd,1].values-1, index_j1].values
            
            # Loop over vertical levels and calculate gradient of field by taking the difference of the values
            # at the adjacent cell centers divided by the distance between the cells
            if len(nLevs) == 1:
                varGrad[edgeInd,0] = cellVar[cellsOnEdge[edgeInd,1].values-1] - cellVar[cellsOnEdge[edgeInd,0].values-1]
                varGrad[edgeInd,0] = varGrad[edgeInd,0] / dcEdge[edgeInd]
                
                # Ensure that the sign of the gradient is consistent with the convention for the u (normal winds). 
                # I think the signs are correct without doing this procedure, but I will keep it here to double check just in case. 
                if varGrad[edgeInd,0] > 0:
                    # Gradient vector points toward cellsOnEdge[edgeInd,1] -> should be directed inward for cellsOnEdge[edgeInd,1]
                    # What is sign of normal vector along edgeInd for each cell?  
                    if sign_j1 < 0:                # Normal vector points inward for cellsOnEdge[edgeInd,1] and outward for cellsOnEdge[edgeInd,0]
                        varGrad[edgeInd,0] = np.abs(varGrad[edgeInd,0])
                    else:
                        varGrad[edgeInd,0] = -np.abs(varGrad[edgeInd,0])
                        print('Changed sign for positive varGrad.')

                elif varGrad[edgeInd,0] < 0:
                    # Gradient vector points toward cellsOnEdge[edgeUse,0] -> should be directed inward for cellsOnEdge[edgeUse,0]
                    # What is sign of normal vector along edgeInd for each cell?  
                    if sign_j0 < 0:                # Normal vector points inward for cellsOnEdge[edgeInd,0] and outward for cellsOnEdge[edgeInd,1]
                        varGrad[edgeInd,0] = np.abs(varGrad[edgeInd,0])
                        print('Changed sign for negative varGrad.')
                    else:
                        varGrad[edgeInd,0] = -np.abs(varGrad[edgeInd,0])
            
            else:
                for k in np.arange(0,len(nLevs),1):
                    varGrad[edgeInd,k] = cellVar[cellsOnEdge[edgeInd,1].values-1,k] - cellVar[cellsOnEdge[edgeInd,0].values-1,k]
                    varGrad[edgeInd,k] = varGrad[edgeInd,k] / dcEdge[edgeInd]
                    
                    # Ensure that the sign of the gradient is consistent with the convention for the u (normal winds). 
                    # I think the signs are correct without doing this procedure, but I will keep it here to double check just in case. 
                    if varGrad[edgeInd,k] > 0:
                        # Gradient vector points toward cellsOnEdge[edgeInd,1] -> should be directed inward for cellsOnEdge[edgeInd,1]
                        # What is sign of normal vector along edgeInd for each cell?  
                        if sign_j1 < 0:                # Normal vector points inward for cellsOnEdge[edgeInd,1] and outward for cellsOnEdge[edgeInd,0]
                            varGrad[edgeInd,k] = np.abs(varGrad[edgeInd,k])
                        else:
                            varGrad[edgeInd,k] = -np.abs(varGrad[edgeInd,k])
                            print('Changed sign for positive varGrad.')
                            
                    elif varGrad[edgeInd,k] < 0:
                        # Gradient vector points toward cellsOnEdge[edgeUse,0] -> should be directed inward for cellsOnEdge[edgeUse,0]
                        # What is sign of normal vector along edgeInd for each cell?  
                        if sign_j0 < 0:                # Normal vector points inward for cellsOnEdge[edgeInd,0] and outward for cellsOnEdge[edgeInd,1]
                            varGrad[edgeInd,k] = np.abs(varGrad[edgeInd,k])
                            print('Changed sign for negative varGrad.')
                        else:
                            varGrad[edgeInd,k] = -np.abs(varGrad[edgeInd,k])
                                                
    return varGrad
    
def mpas_reconstruct_grad(gradEdge, coeffs_reconstruct, edgesOnCell, nEdgesOnCell, orientation='cartesian', latCell=None, lonCell=None):
    '''    
    Given MPAS mesh information, takes a horizontal gradient field valid on cell edges and reconstructs the horizontal gradient vectors at the cell 
    center in a manner analogous to the u reconstruction of mpas_reconstruct_2d in mpas_vector_reconstruction.F

    Parameters
    ----------
    gradEdge : 1D or 2D DataArray or NumPy array of [nEdges] or [nEdges, nVertLevels]
        Gradient of variable on cell edges 
    coeffs_reconstruct : 3D DataArray or NumPy array of [nCells, maxEdges, 3]
        Coefficients to reconstruct velocity vectors at cell centers
        Important note: while this variable exists in the grid static files, the values are all zeroes. This variable should be imported from either the 
        initial conditions or model output files
    nEdgesOnCell : 1D or 2D DataArray or NumPy array of [nCells]
        Number of edges forming the boundary of a cell
    edgesOnCell : 2D DataArray or NumPy array of dims [nCells, maxEdges], required 
        IDs of edges forming the boundary of a cell
    orientation : str, optional
        "cartesian" for x,y,z components, "directional" for zonal and meridional components, or "all" for all components. defaults to cartesian
    latCell : 1D DataArray or NumPy array of [nCells], required for orientation='directional' or orientation='all' 
        Latitude of cell centers
    lonCell : 1D DataArray or NumPy array of [nCells], required for orientation='directional' or orientation='all' 
        Longitude of cell centers
    
    Returns
    -------
    gradReconstruct : 3D NumPy array of [nCells, nVertLevels, nComponents]
        Components of horizontal gradient reconstructed to cell centers
        
        - nComponents and their ordering vary based on orientation:
        orientation = 'cartesian':
            nComponents = 3; d/dx, d/dy, d/dz
        orientation='directional':
            nComponents = 2; d/dx_zonal, d/dy_meridional
        orientation = 'all':
            nComponents = 5; d/dx, d/dy, d/dz, d/dx_zonal, d/dy_meridional    
    '''
    
    # Check number of dimensions of gradEdge
    if gradEdge.ndim > 2:
        raise AttributeError('Variable has more than 2 dimensions. Exiting.')
    else:
        pass
    
    # Need number of cells. 
    if isinstance(edgesOnCell, xr.core.dataarray.DataArray):
        nCells = edgesOnCell.nCells
    elif isinstance(edgesOnCell, np.ndarray):
        nCells = np.arange(0,np.shape(edgesOnCell)[0],1)
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Need number of vertical levels
    if isinstance(gradEdge, xr.core.dataarray.DataArray):
        if gradEdge.ndim == 2:
            nLevs = gradEdge.nVertLevels
        else:
            nLevs = [1]
    elif isinstance(gradEdge, np.ndarray):
        if gradEdge.ndim == 2:
            nLevs = np.arange(0,np.shape(gradEdge)[1],1)
        else:
            nLevs = [1]
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Check if latitude and longitude variables are provided for orientation='all' and orientation='directional'
    if (orientation == 'all') or (orientation == 'directional'):
        if latCell is None or lonCell is None:
             raise RuntimeError('Cell latitude and longitude arrays are needed for this orientation. Exiting.')
        else:
            pass
       
    # Initialize gradReconstruct
    if orientation == 'all':
        gradReconstruct = np.zeros((len(nCells),len(nLevs), 5))
    elif orientation == 'cartesian':
        gradReconstruct = np.zeros((len(nCells),len(nLevs), 3))
    elif orientation == 'directional':
        gradReconstruct = np.zeros((len(nCells),len(nLevs), 2))
    else:
        raise RuntimeError('Unrecognized orientation given. Exiting.')
    
    # Initialize individual component arrays
    gradReconstructX = np.zeros((len(nCells),len(nLevs)))
    gradReconstructY = np.zeros((len(nCells),len(nLevs)))
    gradReconstructZ = np.zeros((len(nCells),len(nLevs)))
    
    if (orientation == 'all') or (orientation == 'directional'):
        gradReconstructZonal = np.zeros((len(nCells),len(nLevs)))
        gradReconstructMeridional = np.zeros((len(nCells),len(nLevs)))
    
    # Loop over cells
    for i in np.arange(0,len(nCells),1):
        
        if (orientation == 'all') or (orientation == 'directional'):
            clat = np.cos(latCell[i].values)
            slat = np.sin(latCell[i].values)
            clon = np.cos(lonCell[i].values)
            slon = np.sin(lonCell[i].values)
            
        # Loop over edges on cell
        for j in np.arange(0,nEdgesOnCell[i].values,1):

            # Edge IDs, inds for j edge along parent i cell
            edgeUse = edgesOnCell[i,j]
            edgeInd = edgeUse - 1
            
            if len(nLevs) == 1:
                gradReconstructX[i,0] =  gradReconstructX[i,0] + coeffs_reconstruct[i,j,0].values * gradEdge[edgeInd]
                gradReconstructY[i,0] =  gradReconstructY[i,0] + coeffs_reconstruct[i,j,1].values * gradEdge[edgeInd]
                gradReconstructZ[i,0] =  gradReconstructZ[i,0] + coeffs_reconstruct[i,j,2].values * gradEdge[edgeInd]
               
            else:
                gradReconstructX[i,:] =  gradReconstructX[i,:] + coeffs_reconstruct[i,j,0].values * gradEdge[edgeInd,:]
                gradReconstructY[i,:] =  gradReconstructY[i,:] + coeffs_reconstruct[i,j,1].values * gradEdge[edgeInd,:]
                gradReconstructZ[i,:] =  gradReconstructZ[i,:] + coeffs_reconstruct[i,j,2].values * gradEdge[edgeInd,:]

            # Store gradient components in gradReconstruct array
            if (orientation == 'all') or (orientation == 'cartesian'):
                gradReconstruct[i,:,0] = gradReconstructX[i,:]
                gradReconstruct[i,:,1] = gradReconstructY[i,:]
                gradReconstruct[i,:,2] = gradReconstructZ[i,:]

            # Calculate zonal/meridional components if requested
            if (orientation == 'all') or (orientation == 'directional'):
                gradReconstructZonal[i,:] = -gradReconstructX[i,:]*slon + gradReconstructY[i,:]*clon
                gradReconstructMeridional[i,:] = -(gradReconstructX[i,:]*clon + gradReconstructY[i,:]*slon)*slat + gradReconstructZ[i,:]*clat

                if (orientation == 'all'):
                    gradReconstruct[i,:,3] = gradReconstructZonal[i,:]
                    gradReconstruct[i,:,4] = gradReconstructMeridional[i,:]

                elif (orientation == 'directional'):
                    gradReconstruct[i,:,0] = gradReconstructZonal[i,:]
                    gradReconstruct[i,:,1] = gradReconstructMeridional[i,:]
                    
        else:
            pass
            
    return gradReconstruct


def mpas_calc_curl(edgeVar, dcEdge, areaTriangle, verticesOnEdge):
    '''    
    Given MPAS mesh information, calculates the finite volume curl at the vertices using a vector on the cell edges.

    For the u wind field, this yields the relative vertical vorticity at the cell vertices. Code is adapted from the computation
    of circulation and relative vorticity at each vertex in atm_compute_solve_diagnostics() in mpas_atm_time_integration.F 
    

    Parameters
    ----------
    edgeVar : 1D or 2D DataArray or NumPy array of [nEdges] or [nEdges, nVertLevels]
        Variable on cell edges 
    dcEdge : 1D DataArray or NumPy array of dims [nEdges]
        Spherical distance between cells separated by an edge
    areaTriangle : 1D DataArray or NumPy array of dims [nVertices]
        Spherical area of a Delaunay triangle
    verticesOnEdge : 2D DataArray or NumPy array of dims [nEdges, TWO]
        IDs of the two vertex endpoints of an edge
    
    Returns
    -------
    curlVert : 1D or 2D NumPy array of [nVertices] or [nVertices, nVertLevels]
        Curl of edgeVar at cell vertices
    '''
    
    # Need number of edges
    if isinstance(dcEdge, xr.core.dataarray.DataArray):
        nEdges = dcEdge.nEdges
    elif isinstance(dcEdge, np.ndarray):
        nEdges = np.arange(0,np.shape(dcEdge)[0],1)
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    # Need number of vertical levels
    if isinstance(edgeVar, xr.core.dataarray.DataArray):
        if edgeVar.ndim == 2:
            nLevs = edgeVar.nVertLevels
        else:
            nLevs = [1]
    elif isinstance(edgeVar, np.ndarray):
        if edgeVar.ndim == 2:
            nLevs = np.arange(0,np.shape(edgeVar)[1],1)
        else:
            nLevs = [1]
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
        
    # Need number of vertices
    if isinstance(areaTriangle, xr.core.dataarray.DataArray):
        nVerts = areaTriangle.nVertices
    elif isinstance(areaTriangle, np.ndarray):
        nVerts = np.arange(0,np.shape(areaTriangle)[0],1)
    else:
        raise TypeError('Variable type is neither an Xarray DataArray nor a NumPy array. Exiting.')
    
    #Initialize array for curlVert
    curlVert = np.zeros((len(nVerts),len(nLevs)))
    
    # Loop over cell edges
    if len(nLevs) == 1:
        for j in np.arange(0,len(nEdges),1):
            curlVert[verticesOnEdge[j,0].values-1,0] = curlVert[verticesOnEdge[j,0].values-1,0] - dcEdge[j] * edgeVar[j]
            curlVert[verticesOnEdge[j,1].values-1,0] = curlVert[verticesOnEdge[j,1].values-1,0] + dcEdge[j] * edgeVar[j]
    
        # Now loop over vertices after all edge contributions have been taken into account and divide by areaTriangle
        for v in np.arange(0,len(nVerts),1):
            curlVert[v,0] = curlVert[v,0] / areaTriangle[v].values
            
    else:
        for j in np.arange(0,len(nEdges),1):
            curlVert[verticesOnEdge[j,0].values-1,:] = curlVert[verticesOnEdge[j,0].values-1,:] - dcEdge[j] * edgeVar[j]
            curlVert[verticesOnEdge[j,1].values-1,:] = curlVert[verticesOnEdge[j,1].values-1,:] + dcEdge[j] * edgeVar[j]
    
        # Now loop over vertices after all edge contributions have been taken into account and divide by areaTriangle
        for v in np.arange(0,len(nVerts),1):
            curlVert[v,:] = curlVert[v,:] / areaTriangle[v].values
        
    return curlVert
