function slope_val = slope(grid, values, left_slope, right_slope, left, right)
    if (left == 0) 
        slope_val = left_slope;
    elseif (right == length(grid) + 1) 
        slope_val = right_slope;
    elseif (values(right) == inf && values(left) == inf)  %! when both values are Inf
        slope_val = inf;                                           %! we want to have an Inf slope
    else
        slope_val = (values(right) - values(left)) / (grid(right) - grid(left));
    end 
end

