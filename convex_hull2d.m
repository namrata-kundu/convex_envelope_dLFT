function new_values = convex_hull2d(values, xgrid, ygrid, parameters)
%CONVEX_HULL2D Summary of this function goes here
%   Detailed explanation goes here

    if parameters('doubleDLFT')
        new_values = convex_hull2d_doubleDLFT(values, xgrid, ygrid, parameters);
    elseif ~parameters('symmetric')
        if parameters('xy_order')
            new_values = convex_hull2d_asymmetric(values, xgrid, ygrid, parameters);
        else
            new_values_transpose = convex_hull2d_asymmetric(values', ygrid, xgrid, parameters);
            new_values = new_values_transpose';
        end
    else
        new_values = convex_hull2d_asymmetric(values, xgrid, ygrid, parameters);
        sym_values_transpose = convex_hull2d_asymmetric(values', ygrid, xgrid, parameters);
        sym_values = sym_values_transpose';
        new_values = max(new_values, sym_values);
    end
                
end

