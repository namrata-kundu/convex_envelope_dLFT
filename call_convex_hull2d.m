function new_values = call_convex_hull2d(values, xgrid, ygrid,symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb)
    parameters = containers.Map();
    parameters('symmetric')=symmetric;
    parameters('doubleDLFT')=doubleDLFT;
    parameters('xy_order')=xy_order;
    parameters('dual_grid_type')=dual_grid_type;
    parameters('dg_ext_slopes')=dg_ext_slopes;
    parameters('dg_inject')=dg_inject;
    parameters('dg_inject_value')=dg_inject_value;
    parameters('dg_slices')=dg_slices;
    parameters('dg_only_slice')=dg_only_slice;
    parameters('dg_size')=dg_size;
    parameters('dg_no_auto')=dg_no_auto;
    parameters('dg_lwb')=dg_lwb;
    parameters('dg_upb')=dg_upb;
    
    new_values = ConvexHull2D.convex_hull2d(values, xgrid, ygrid, parameters)
 
end

