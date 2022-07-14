function [lwb, upb] = modify_dual_grid_bounds(lwb, upb, parameters)
    if (parameters('dg_no_auto') == 1) 
        lwb = lwb +    parameters('dg_lwb')    * (upb - lwb);
        upb = upb + (parameters('dg_upb') - 1) * (upb - lwb) / (1 - parameters('dg_lwb'));
    elseif (parameters('dg_no_auto') == 2)
        lwb = parameters('dg_lwb');
        upb = parameters('dg_upb');
    end
%     last_lwb = lwb;
%     last_upb = upb;

end

