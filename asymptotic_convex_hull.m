function [asymptotic_rates, fixation_time] = asymptotic_convex_hull(grid, values, rates, asymptotic_rates, fixation_time)

    hull_indices = zeros(1,length(grid));

    for j=1:length(grid)
        if (values(j) ~= Inf()) 
            break;
        end
    end 
    if (j == length(grid) + 1) % ! all values are Inf()
        asymptotic_rates = inf;
        fixation_time    = 0;
        return;
    end

    hull_indices(1) = j; %! we select initially the first non-Inf grid point
    j               = 1;
    fixation_time   = 0; %! fixation time of the first point is 0

    while true
        if (hull_indices(j) == length(grid)) 
            break;
        end

        hull_indices(j+1) = hull_indices(j) + 1;
        minvalue = (rates(hull_indices(j+1)) - rates(hull_indices(j))) / (grid(hull_indices(j+1)) - grid(hull_indices(j)));
        for i=hull_indices(j)+2:length(grid)
            newvalue = (rates(i) - rates(hull_indices(j))) / (grid(i) - grid(hull_indices(j)));
            if (newvalue < minvalue) 
                minvalue = newvalue;
                hull_indices(j+1) = i;
            elseif (newvalue == minvalue) 
                if ( (values(hull_indices(j+1)) - values(hull_indices(j))) / (grid(hull_indices(j+1)) - grid(hull_indices(j))) ...
                     >= (values(i) - values(hull_indices(j))) / (grid(i) - grid(hull_indices(j))) ) 
                    hull_indices(j+1) = i;
                end
            end 
        end 
        
        arr = diag( (   diag((values(hull_indices(j+1))  - values(hull_indices(j)))  ...
                 ./ (grid(hull_indices(j+1))    - grid(hull_indices(j))) )   ...
                 - diag((values(hull_indices(j)+1:end) - values(hull_indices(j)))  ...
                 ./ (grid(hull_indices(j)+1:end)   - grid(hull_indices(j))) )   ...
               )                                                           ...
             ./ (  diag( (rates(hull_indices(j)+1:end)  - rates(hull_indices(j)))   ...
                 ./ (grid(hull_indices(j)+1:end)   - grid(hull_indices(j))) )   ...
                 - diag ((rates(hull_indices(j+1))   - rates(hull_indices(j)))   ...
                 ./ (grid(hull_indices(j+1))    - grid(hull_indices(j))) )   ...
               ));

        fixation_time = max(fixation_time, maxval_except(  arr, except=hull_indices(j+1)-hull_indices(j)));
        j = j + 1;

    end

    asymptotic_rates = piecewise_linear_interpolation( grid(hull_indices(1:j)), rates((hull_indices(1:j))),  ...
                                                       -inf, inf, grid                                );

end

