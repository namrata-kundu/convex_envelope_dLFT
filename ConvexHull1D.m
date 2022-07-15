classdef ConvexHull1D
    %CONVEXHULL1D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods(Static)
        
        function index = bisect(array, x)
            if (length(array) == 1)
                index = 1;
                return
            end

            if (x == array(length(array))) 
                index = length(array);
                return
            end

            a = 1;
            b = length(array);
            while true
                index = floor((a + b) / 2);
                if (array(index) <= x && array(index + 1) > x)
                    return;
                elseif (array(index) > x) 
                    b = index;
                else
                    a = index;
                end
            end
        end
        
        
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
                    left = ConvexHull1D.bisect(grid, x(i));
                    if (x(i) == grid(left)) 
                        y(i) = values(left);
                    else
                        shift = (x(i) - grid(left)) / (grid(left + 1) - grid(left));
                        y(i)  = values(left) + shift * (values(left + 1) - values(left)); % more precise than the version with two occurences of shift when the two values are equal
                    end 
                end 
            end
        end
        
        
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

                if ( ConvexHull1D.slope(grid, values, left_slope, right_slope, sp(j - 1), sp(j)) >= ConvexHull1D.slope(grid, values, left_slope, right_slope, sp(j), i)   ) 
                    while true
                        if (j == 2 && i == length(grid) + 1) 
                            break; %! if the two ext slopes are identical, we want at least one internal point to survive
                        end
                        j = j - 1;
                        if (j == 1) 
                            break;
                        end
                        if ( ConvexHull1D.slope(grid, values, left_slope, right_slope, sp(j - 1), sp(j)) < ConvexHull1D.slope(grid, values, left_slope, right_slope, sp(j), i)  ) 
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


        function pcw = legendre_conjugate(grid, values, left_slope, right_slope, pcw)
        %LEGENDRE_CONJUGATE Summary of this function goes here
        %   Detailed explanation goes here
        %     pcw = PiecewiseLinear_t();
            cgrid=zeros(1,10000);
            cvalues=zeros(10000,1);
            [cgrid, cvalues] = ConvexHull1D.convex_hull(grid, values, left_slope, right_slope, cgrid, cvalues); %%changed parameters
            if (cvalues(1) == inf) 
                pcw.grid   = 0; %(/   0    /);
                pcw.values = -inf; %(/ -Inf() /);
                pcw.left_slope  = inf;
                pcw.right_slope = -inf;
            elseif (cvalues(1) == -inf) 
                pcw.grid   = 0; % (/   0   /);
                pcw.values = inf; % (/ inf /);
                pcw.left_slope  = -inf;
                pcw.right_slope = inf;
        %     !else if (size(cgrid) == 1) then ! fails if there are external slopes
        %     !    pcw%grid   = (/      0      /)
        %     !    pcw%values = (/ -cvalues(1) /)
        %     !    pcw%left_slope  = cgrid(1)
        %     !    pcw%right_slope = cgrid(1)
            else                   
                new_size = length(cgrid) - 1;
                lb = 1;
                ub = new_size;
                if (left_slope > -inf)  
                    new_size = new_size + 1;
                    lb = 2;
                    ub = ub + 1;
                end 
                if (right_slope < Inf()) 
                    new_size = new_size + 1;
                end 

                if (new_size == 0)  % in this case we have simply a line
                    pcw.grid   = 0;% (/      0      /);
                    pcw.values = -cvalues(1); %(/ -cvalues(1) /);
                    pcw.left_slope  = cgrid(1);
                    pcw.right_slope = cgrid(1);
                else
        %             allocate(pcw%grid(new_size))
        %             allocate(pcw%values(new_size))
                    pcw.grid = zeros(1,new_size);
                    pcw.values = zeros(1,new_size);

                    intermediate_grid = (cvalues(2:end) - cvalues(1:length(cvalues) - 1)) ./ (cgrid(2:end)   - cgrid(1:length(cgrid) - 1)); %%%changed here
        %             pcw.grid(lb:ub) = intermediate_grid(:,1);
                    intermediate_grid_size = size(intermediate_grid);
                    if intermediate_grid_size(1) == 1
                        pcw.grid(lb:ub) = intermediate_grid;
                    else
                        pcw.grid(lb:ub) = diag(intermediate_grid);
                    end

        %             pcw.grid(lb:ub)     =   (cvalues(2:end) - cvalues(1:length(cvalues) - 1)) ./ (cgrid(2:end)   - cgrid(1:length(cgrid) - 1)); %%%changed here
        %             umm=pcw.grid';
                    intermediate_values = cgrid(1:length(cgrid) - 1) .* pcw.grid(lb:ub)  - cvalues(1:length(cvalues) - 1);
                    intermediate_values_size = size(intermediate_values);
                    if intermediate_values_size(1) == 1
                        pcw.values(lb:ub)   =   intermediate_values;
                    else
                        pcw.values(lb:ub)   =   diag(intermediate_values);
                    end

                    if (left_slope > -inf) 
                        pcw.grid(1)    = left_slope;
                        pcw.values(1)  = cgrid(1) * left_slope - cvalues(1);
                        pcw.left_slope = -inf;
                    else
                        pcw.left_slope = cgrid(1);
                    end 

                    if (right_slope < Inf()) 
                        pcw.grid(new_size)   = right_slope;
                        pcw.values(new_size) = cgrid(length(cgrid)) * right_slope - cvalues(length(cvalues));
                        pcw.right_slope      = inf;
                    else
                        pcw.right_slope      = cgrid(length(cgrid));
                    end
                end
            end


        end




        
    end
end

