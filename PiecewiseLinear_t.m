classdef PiecewiseLinear_t
    %PIECEWISELINEAR_T Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        grid,
        values,
        left_slope,
        right_slope
    end
    
    methods
%         function obj = PiecewiseLinear_t(inputArg1,inputArg2)
%             %PIECEWISELINEAR_T Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.grid = inputArg1 + inputArg2;
%         end
        
        function y = evaluate(obj,x)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            y = ConvexHull1D.piecewise_linear_interpolation(obj.grid, obj.values, obj.left_slope, obj.right_slope, x);
        end
        
    end
end

