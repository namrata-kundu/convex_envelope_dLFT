classdef ConvexHull2DTest < matlab.unittest.TestCase
    %CONVEXHULL2DTEST Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Test)
        function test_convex_hull2d_s_dDLFT(testCase) 
            xgrid = (-2:1:2);
            ygrid=(-2:1:2);
            values = function1(xgrid,ygrid);
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
            test_values = s_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, s_dDLFT_obj.parameters);
            actual_values = [49    16     9    16    49; ...
                            16     0     0     0    16; ...
                             9     0     0     0     9; ...
                            16     0     0     0    16; ...
                            49    16     9    16    49] ;

            testCase.verifyEqual(test_values,actual_values);
        end
        
        function test_convex_hull2d_a_dDLFT(testCase) 
            xgrid = (-2:1:2);
            ygrid=(-2:1:2);
            values = function1(xgrid,ygrid);
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
            test_values = a_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, a_dDLFT_obj.parameters);
            actual_values = [49    16     9    16    49; ...
                            16     0     0     0    16; ...
                             9     0     0     0     9; ...
                            16     0     0     0    16; ...
                            49    16     9    16    49] ;

            testCase.verifyEqual(test_values,actual_values);
        end
        
        function test_convex_hull2d_ma_dDLFT(testCase) 
            xgrid = (-2:1:2);
            ygrid=(-2:1:2);
            values = function1(xgrid,ygrid);
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
            test_values = ma_dDLFT_obj.convex_hull2d(values, xgrid, ygrid, ma_dDLFT_obj.parameters);
            actual_values = [49    16     9    16    49; ...
                            16     0     0     0    16; ...
                             9     0     0     0     9; ...
                            16     0     0     0    16; ...
                            49    16     9    16    49] ;

            testCase.verifyEqual(test_values,actual_values);
        end
        
    end
end

