function out_values = stage2(pcws, xgrid, ygrid, out_values,ugrid,left_slopes,right_slopes)

    pcw = PiecewiseLinear_t();

    for i = 1:length(ygrid)
        out_values(1:length(ugrid),i) = - pcws(i).evaluate(ugrid);
    end

    pcw.left_slope  = -inf;
    pcw.right_slope =  inf; %! needed so we can evaluate correctly all cases
    for i = 1:length(ugrid)
        [pcw.grid, pcw.values] = convex_hull(ygrid, out_values(i,:), -Inf(), Inf(), pcw.grid, pcw.values);
        out_values(i,:) = - pcw.evaluate(ygrid);
    end 

    if (length(left_slopes)) 
        for i = 1:length(ygrid)
            pcw = legendre_conjugate(ugrid, out_values(1:length(ugrid),i), left_slopes(i), -right_slopes(i), pcw);
%             ! we did not change sign for the left_slope since that was already done
%             ! in the asymptotic hull since "time" is inverted on the left
            out_values(:, i) = pcw.evaluate(xgrid);
        end
    else
        for i = 1:length(ygrid)
            pcw = legendre_conjugate(ugrid, out_values(1:length(ugrid),i), -inf, inf, pcw);
            out_values(:, i) = pcw.evaluate(xgrid);
        end
    end 
end

