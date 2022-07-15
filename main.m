%%Running examples

xgrid = (-1.5:0.1:1.5);
ygrid=(-1.5:0.1:1.5);
run_example(xgrid, ygrid, @function1)


xgrid = (-3:0.1:3);
ygrid=(-3:0.1:3);
run_example(xgrid, ygrid, @function2)

xgrid = (-3:0.02:3);
ygrid = (-3:0.02:3);
run_example(xgrid, ygrid, @function4)



function run_example(xgrid,ygrid, function_to_run)

    values = function_to_run(xgrid,ygrid);
    mesh(xgrid,ygrid,values)
    
    s_dDLFT_new_values = run_s_dDLFT(xgrid,ygrid,values);
    mesh(xgrid,ygrid,s_dDLFT_new_values) 
    
    a_dDLFT_new_values = run_a_dDLFT(xgrid,ygrid,values);
    mesh(xgrid,ygrid,a_dDLFT_new_values) 

    ma_dDLFT_new_values = run_ma_dDLFT(xgrid,ygrid,values);
    mesh(xgrid,ygrid,ma_dDLFT_new_values) 
    
end


function s_dDLFT_new_values = run_s_dDLFT(xgrid,ygrid,values)
    symmetric = false;
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

    doubleDLFT = true; %standard dlft
    s_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
    s_dDLFT_new_values = s_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, s_dDLFT_obj.parameters);
end

function a_dDLFT_new_values = run_a_dDLFT(xgrid,ygrid,values)
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

    doubleDLFT = false;
    symmetric = false; %alternating algo
    a_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
    a_dDLFT_new_values = a_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, a_dDLFT_obj.parameters);
end

function ma_dDLFT_new_values = run_ma_dDLFT(xgrid,ygrid,values)
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

    doubleDLFT = false;
    symmetric = true; %max alternating algo
    ma_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
    ma_dDLFT_new_values = ma_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, ma_dDLFT_obj.parameters);
end


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
% xgrid = [-3:0.004:3];
% ygrid = [-3:0.004:3];
% values = function4(xgrid,ygrid);
% mesh(xgrid,ygrid,values)
% 
% % symmetric = true; %max alternating algo
% symmetric = false; %alternating algo
% % doubleDLFT = false;
% xy_order=true;
% dual_grid_type=0;
% dg_ext_slopes=false;
% dg_inject=false;
% dg_inject_value=0;
% dg_slices=10;
% dg_only_slice=0;
% dg_size=1;
% dg_no_auto=0;
% dg_lwb=0;
% dg_upb=1;
% 
% doubleDLFT = true; %standard dlft
% 
% s_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
% s_dDLFT_new_values = s_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, s_dDLFT_obj.parameters);
% mesh(xgrid,ygrid,s_dDLFT_new_values) 
% 
% doubleDLFT = false;
% symmetric = false; %alternating algo
% 
% a_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
% a_dDLFT_new_values = a_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, a_dDLFT_obj.parameters);
% mesh(xgrid,ygrid,a_dDLFT_new_values) 
% 
% doubleDLFT = false;
% symmetric = true; %max alternating algo
% 
% ma_dDLFT_obj = ConvexHull2D(symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb);
% ma_dDLFT_new_values = ma_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, ma_dDLFT_obj.parameters);
% mesh(xgrid,ygrid,ma_dDLFT_new_values) 

% new_values = call_convex_hull2d(values, xgrid, ygrid,symmetric,doubleDLFT,xy_order,dual_grid_type,dg_ext_slopes,dg_inject,dg_inject_value, dg_slices,dg_only_slice,dg_size,dg_no_auto,dg_lwb,dg_upb)
% mesh(xgrid,ygrid,new_values) 