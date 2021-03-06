% xgrid = [-1.5:0.1:1.5];
% ygrid=[-1.5:0.1:1.5];
% xgrid = [-1.5:0.003:1.5];
% ygrid=[-1.5:0.003:1.5];
% values = function1(xgrid,ygrid);

% xgrid = [-3:0.1:3];
% ygrid=[-3:0.1:3];
% values = function2(xgrid,ygrid);

% xgrid = [-2:1:2];
% ygrid=[-2:1:2];
% xgrid = [-1:0.05:1];
% ygrid = [-1:0.05:1];
xgrid = [-3:0.004:3];
ygrid = [-3:0.004:3];
values = function4(xgrid,ygrid);
mesh(xgrid,ygrid,values)
% scatter3(xgrid,ygrid,values)


symmetric = true; %max alternating algo
% symmetric = false; %alternating algo
doubleDLFT = false;
% doubleDLFT = true; %standard dlft
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