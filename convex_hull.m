function [new_grid, new_values] = convex_hull(grid, values, left_slope, right_slope, new_grid, new_values)
%CONVEX_HULL Summary of this function goes here
%   Detailed explanation goes here

% ! NB we expect that if the first value is Inf, then left_slope=Inf()
%         !    and similarly on the right side (otherwise we should also output new external slopes)
%         
%         ! remove starting Infs
    i = 1;
    while true
        if (i > length(values)) 
            new_grid   = 0; %(/   0   /);
            new_values = inf; %(/ inf /);
            return;
        end
        if (values(i) ~= inf) 
            break;
        end
        i = i + 1;
    end 

%     ! remove ending Infs
    imax = length(values);
    while true
        if (values(imax) ~= inf) 
            break;
        end
        imax = imax - 1;
    end 
%     ! deal with the case in which we have only one point
%     ! not really needed
%     !if (i == imax) then
%     !    new_grid   = (/  grid(i)  /)
%     !    new_values = (/ values(i) /)
%     !    return
%     !end if

    sp = [ 0, i ];%(/ 0, i /); % ! selected_points
    i = i+1;         %! index inside grid
    j = 2;           %! index inside sp

    while true
        if (i ~= length(grid) + 1) 
            if (i > imax)  %! when imax == size(grid) we have to consider the right_slope
                break;
            end 

            if (values(i) == -inf) 
                new_grid   = 0; %(/   0    /);
                new_values = -inf; %(/ -inf /)
                return; %! the convex hull can be actually different, but the LFT treat every function which is 
            end      %! -Inf in at least one point in the same way. Then it is ok,
                       %! since we never need such a convex hull without a LFT (if the input data is > -Inf)
        end

        if ( slope(grid, values, left_slope, right_slope, sp(j - 1), sp(j)) >= slope(grid, values, left_slope, right_slope, sp(j), i)   ) 
            while true
                if (j == 2 && i == length(grid) + 1) 
                    break; %! if the two ext slopes are identical, we want at least one internal point to survive
                end
                j = j - 1;
                if (j == 1) 
                    break;
                end
                if ( slope(grid, values, left_slope, right_slope, sp(j - 1), sp(j)) < slope(grid, values, left_slope, right_slope, sp(j), i)  ) 
                    break;
                end 
            end 
        end 

        j = j + 1;
        sp(j) = i;
        i = i + 1;
    end 

    if (i ~= length(grid) + 2) 
        j = j + 1; %! pretend we added the right_slope
    end 
    new_grid   = grid(sp(2:j-1));
    new_values = values(sp(2:j-1));
end

