function y = piecewise_linear_interpolation(grid, values, left_slope, right_slope, x)
%PIECEWISE_LINEAR_INTERPOLATION Summary of this function goes here
%   Detailed explanation goes here
    y= zeros(1, length(x));
    for i = 1: length(x)
        if (x(i) < grid(1)) % cannot use <= here since if slope=Inf and x-grid(1)=0 we have NaN
            y(i) = values(1) + left_slope * (x(i) - grid(1));
        elseif (x(i) > grid(length(grid))) 
            y(i) = values(length(grid)) + right_slope * (x(i) - grid(length(grid)));
        else
            left = bisect(grid, x(i));
            if (x(i) == grid(left)) 
                y(i) = values(left);
            else
                shift = (x(i) - grid(left)) / (grid(left + 1) - grid(left));
                y(i)  = values(left) + shift * (values(left + 1) - values(left)); % more precise than the version with two occurences of shift when the two values are equal
            end 
        end 
    end
end

