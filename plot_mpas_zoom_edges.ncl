load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
;load "./gsn_csm.ncl"
;=============================================================================================
;----------------------------------------------------------------------
; Script to plot global MPAS cell and edge fields on native grid over a
; zoomed-in region of interest and color the cell edges based on value
;
; Manda Chasteen - NCAR/MMM 
; Sept 14, 2022
;----------------------------------------------------------------------
begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DEFINE WORKSTATION:

outfile = "theta_gradient_edge_magnitude_lev7"

wks = gsn_open_wks("x11",outfile)

setvalues NhlGetWorkspaceObjectId()
     "wsMaximumSize": 400000000
end setvalues
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; INPUT FILES:

f_hist = addfile("history.2017-03-05_00.01.15.nc","r")
f_grid = addfile("x1.2621442.static.nc","r")

; SPECIFIC VARIABLES: 
;;; I have output the magnitude of the theta difference between
;;; adjacent cells on each cell edge: edge_theta_gradient(nVertLevels, nEdges) 

fe = addfile("tempDiff_edges.nc","r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SPECIFY DESIRED DOMAIN BOUNDS:

lat_min = 35.
lat_max = 37.
lon_min = 360.-102.
lon_max = 360.-97.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; IMPORT VARIABLES:

zlev = 6  

; Grid Vars:
ind_edgeID     = f_grid->indexToEdgeID
ind_cellID     = f_grid->indexToCellID
ind_vertID     = f_grid->indexToVertexID

r2d            = get_r2d("double")
latCell        = f_grid->latCell * r2d
lonCell        = f_grid->lonCell * r2d
lonVertex      = f_grid->lonVertex * r2d 
latVertex      = f_grid->latVertex * r2d
latEdge        = f_grid->latEdge * r2d
lonEdge        = f_grid->lonEdge * r2d

edgesOnCell    = f_grid->edgesOnCell
verticesOnCell = f_grid->verticesOnCell
verticesOnEdge = f_grid->verticesOnEdge
maxEdges       = dimsizes(edgesOnCell(0,:))

nEdgesOnCell    = f_grid->nEdgesOnCell
nEdges          = dimsizes(latEdge)
nVertices       = dimsizes(latVertex)
nCells          = dimsizes(latCell)

; Note: for some reason, lonVertex is from -+180, while lonCell and lonEdge are from 0-360
lonVertex       =  where(lonVertex.lt.0,lonVertex+360,lonVertex)

; Meteorological Vars:
theta           = f_hist->theta(0,:,zlev)          ; dims are nCells
thetaGrad       = fe->edge_theta_gradient(zlev,:)  ; dims are nEdges

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SUBSET DATA WE WISH TO PLOT

;;; cell indices for cells within lat/lon box
map_ind    =  ind((latCell.ge.lat_min-2).and.(latCell.le.lat_max+2).and.\
                 (lonCell.ge.lon_min-2).and.(lonCell.le.lon_max+2))

;;; edge indices for edges within lat/lon box
map_ind_e  =  ind((latEdge.ge.lat_min-2).and.(latEdge.le.lat_max+2).and.\
                 (lonEdge.ge.lon_min-2).and.(lonEdge.le.lon_max+2))

;;; vertex indices for vertices within lat/lon box
map_ind_v  = ind((latVertex.ge.lat_min-2).and.(latVertex.le.lat_max+2).and.\
                 (lonVertex.ge.lon_min-2).and.(lonVertex.le.lon_max+2))

theta_subset          = theta(map_ind)
thetaGrad_subset      = thetaGrad(map_ind_e)

latCell_subset        = latCell(map_ind)
lonCell_subset        = lonCell(map_ind)
latEdge_subset        = latEdge(map_ind_e)
lonEdge_subset        = lonEdge(map_ind_e)
;latVertex_subset     = latVertex(map_ind_v)
;lonVertex_subset     = lonVertex(map_ind_v)

edgesOnCell_subset    = edgesOnCell(map_ind,:)         ; edgesOnCell for cell centers within lat/lon bounds
verticesOnEdge_subset = verticesOnEdge(map_ind_e,:)
verticesOnCell_subset = verticesOnCell(map_ind,:)

ncells_subset         = dimsizes(latCell_subset)
nedges_subset         = dimsizes(latEdge_subset)
;nverts_subset        = dimsizes(latVertex_subset)

lon_verts = new((/ncells_subset, maxEdges/), double)
lat_verts = new((/ncells_subset, maxEdges/), double)

do in=0,dimsizes(map_ind)-1
   verts    = verticesOnCell(map_ind(in),:)
   vert_ind = verts - 1

   do v=0,maxEdges-1;      
      v_use = vert_ind(v)
      if (v_use.gt.0) then
         lon_use         = lonVertex(v_use)
         lat_use         = latVertex(v_use)
         lat_verts(in,v) = lat_use
         lon_verts(in,v) = lon_use
         vkeep = v
     else 
         lat_verts(in,v) = lat_verts(in,vkeep)
         lon_verts(in,v) = lon_verts(in,vkeep)
     end if           
   end do
end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PRINT SOME INFORMATION

print("Original # cells       : " + nCells)
print("Subsetted # cells      : " + ncells_subset)
print("Original # edges       : " + nEdges)
print("Subsetted # edges      : " + nedges_subset)
;print("Original # vertices    : " + nVertices)
;print("Subsetted # vertices   : " + nverts_subset)

print("Original theta min/max  : " + min(theta) + "/" + max(theta))
print("Subsetted theta min/max : " + min(theta_subset) + "/" + max(theta_subset))

print("Original lat min/max   : " + min(latCell) + "/" + max(latCell))
print("Subsetted lat min/max  : " + min(latCell_subset) + "/" + max(latCell_subset))

print("Original lon min/max   : " + min(lonCell) + "/" + max(lonCell))
print("Subsetted lon min/max  : " + min(lonCell_subset) + "/" + max(lonCell_subset))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SET PLOTTING RESOURCES:

res = True
res@gsnDraw                     = False
res@gsnFrame                    = False
res@gsnMaximize                 = True
res@gsnAddCyclic                = False 
res@tmYROn                      = False
res@tmXTOn                      = False

res@mpDataBaseVersion           = "MediumRes" 
res@mpOutlineBoundarySets       = "GeophysicalAndUSStates"
res@mpNationalLineThicknessF    = 3.0      
res@mpGeophysicalLineThicknessF = 3.0
res@mpUSStateLineThicknessF     = 3.0
res@mpUSStateLineColor          = "navyblue"

res@mpLimitMode                 =  "LatLon"
res@mpMinLatF                   =  lat_min
res@mpMaxLatF                   =  lat_max
res@mpMinLonF                   =  lon_min
res@mpMaxLonF                   =  lon_max
res@mpCenterLonF                =  (lon_max+lon_min)*.5

res@sfXArray                    = lonCell_subset
res@sfYArray                    = latCell_subset
res@sfXCellBounds               = lon_verts  ; necessary for CellFill
res@sfYCellBounds               = lat_verts 

res@cnFillMode    	        = "CellFill"      ; Options: CellFill, RasterFill, AreaFill, PolygonFill
res@cnFillOn      	        = True
res@cnLinesOn     	        = False
res@cnInfoLabelOn 	        = False
res@cnFillPalette               = "MPL_Set1" 
res@gsnSpreadColors             = True
res@gsnSpreadColorStart         = 15
res@gsnSpreadColorEnd           = 122
res@cnLevelSelectionMode        = "ManualLevels"
res@cnMinLevelValF  	        = 285.0
res@cnMaxLevelValF  	        = 305.0
res@cnLevelSpacingF 	        = 0.5
res@gsnLeftString    	        = "Potential temperature and edge gradient magnitudes on model level "+tostring(zlev+1) 
res@gsnRightString  	        = "K"

; Draw cell filled theta
plot = gsn_csm_contour_map(wks,theta_subset,res)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; For coloring edges by edge gradient magnitude values, normalize colors in colormap to values of thetaGrad_subset

gradMax = max(thetaGrad_subset)
gradMin = 0.0

thetaGrad_norm = thetaGrad_subset/gradMax

cmap = read_colormap_file("MPL_gist_gray")
cmap_r    = cmap(::-1,:)                        ; reverse colormap so black = higher values
cmap_u    = cmap_r(5:,:)
cdims_all = dimsizes(cmap_u)
cdims     = cdims_all(0)-1                      ; ncolors in trimmed cmap

colors_grad      = cdims * thetaGrad_norm       ; array indicating where normalized gradient values lie on colormap index range
colors_grad_inds = round(colors_grad,3)         ; array of indices to be applied to colormap

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Initialize array containing all colors to use for cell edges
;;; 
;;; Note: I'm not entirely sure why, but the array needs to be twice the size of the 
;;; number of edges for NCL to plot the edge values in the correct location using
;;; lnres@gsSegments = ispan(0,nedges_subset*2-1,2)  
;;; The do loop below will apply colors corresponding to the first edge to colors_use(0:1,:), 
;;; the second edge to colors_use(2:3,:), etc. 

colors_use = new((/dimsizes(colors_grad_inds)*2+1,4/),float);, "integer")

do col=0,dimsizes(colors_grad_inds)-1   
   col0                    = 2*col
   col1                    = col0 + 2
   colors_use(col0:col1,0) = cmap_u(colors_grad_inds(col),0)
   colors_use(col0:col1,1) = cmap_u(colors_grad_inds(col),1)
   colors_use(col0:col1,2) = cmap_u(colors_grad_inds(col),2)
   colors_use(col0:col1,3) = cmap_u(colors_grad_inds(col),3)
end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;---Code to plot the MPAS edge lines on the existing map
lnres                      = True
;lnres@gsLineColor         = "NavyBlue"       ; set to plot cell edges using single color
lnres@gsLineThicknessF     = 3.75

ecx = new((/nedges_subset,2/),double)
ecy = new((/nedges_subset,2/),double)

ecx(:,0) = lonVertex(verticesOnEdge(map_ind_e,0)-1)
ecx(:,1) = lonVertex(verticesOnEdge(map_ind_e,1)-1)
ecy(:,0) = latVertex(verticesOnEdge(map_ind_e,0)-1)
ecy(:,1) = latVertex(verticesOnEdge(map_ind_e,1)-1)

;;; I think the following code would only be important for plotting areas near prime meridian
;ii0 = ind((abs(ecx(:,0)-ecx(:,1)).gt.180.and.(ecx(:,0).gt.ecx(:,1))))
;ii1 = ind((abs(ecx(:,0)-ecx(:,1)).gt.180.and.(ecx(:,0).lt.ecx(:,1))))
;ecx(ii0,0) = ecx(ii0,0) - 360.0
;ecx(ii1,1) = ecx(ii1,1) - 360.0

print("Attaching the polylines...")
lnres@gsSegments   = ispan(0,nedges_subset*2-1,2) 
lnres@gsColors     = colors_use
  
poly = gsn_add_polyline(wks,plot,ndtooned(ecx),ndtooned(ecy),lnres)

draw(plot)
frame(wks)


end

;;; Other resources: https://www.ncl.ucar.edu/Applications/Scripts/mpas_cell_3.ncl

;=============================================================================================
