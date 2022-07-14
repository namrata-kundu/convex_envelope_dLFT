function dgrid = build_uniform_dual_grid(dgrid, lwb, upb, N, parameters)
    [lwb, upb] = modify_dual_grid_bounds(lwb, upb, parameters);
    if (lwb == upb)
        dgrid = zeros(1,1);
    else
        dgrid_size = (min(max(ceil(N*parameters('dg_size')),1),N));
        dgrid = zeros(1,dgrid_size);
    end

    for i = 1:length(dgrid)-1
        dgrid(i) = lwb + (i - 1) * (upb - lwb) / (length(dgrid) - 1);
    end
    dgrid(length(dgrid)) = upb;  %! in order to have the exact endpoint

    if (parameters('dg_inject'))
        dgrid = injector(dgrid, parameters('dg_inject_value'));
    end 

end

