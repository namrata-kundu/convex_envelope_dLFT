function pcw = legendre_conjugate(grid, values, left_slope, right_slope, pcw)
%LEGENDRE_CONJUGATE Summary of this function goes here
%   Detailed explanation goes here
%     pcw = PiecewiseLinear_t();
    cgrid=zeros(1,10000);
    cvalues=zeros(10000,1);
    [cgrid, cvalues] = convex_hull(grid, values, left_slope, right_slope, cgrid, cvalues); %%changed parameters
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

