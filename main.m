% xgrid = [-1.5:0.1:1.5];
% ygrid=[-1.5:0.1:1.5];
xgrid = [-3:0.1:3];
ygrid=[-3:0.1:3];
% xgrid = [-2:1:2];
% ygrid=[-2:1:2];
values = function2(xgrid,ygrid);
mesh(xgrid,ygrid,values) 
symmetric = true;
% doubleDLFT = false;
doubleDLFT = true;
xy_order=true;
dual_grid_type=0;
dg_ext_slopes=false;
dg_inject=false;
dg_inject_value=0;
dg_slices=10;
dg_only_slice=0;
dg_size=1;
dg_no_auto=0;
dg_lwb=0;
dg_upb=1;
new_values = call_convex_hull2d(values, xgrid, ygrid,symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb)
mesh(xgrid,ygrid,new_values) 