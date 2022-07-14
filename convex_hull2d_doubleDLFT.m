function new_values = convex_hull2d_doubleDLFT(values, xgrid, ygrid, parameters)

    new_values = zeros(length(xgrid),length(ygrid));
    xdgrid = zeros(1,1000000);
    ydgrid = zeros(1,1000000);
    pcws_size = max(length(xgrid),length(ygrid));
    pcws(1:pcws_size) = PiecewiseLinear_t;
    pcw = PiecewiseLinear_t();
    
    lwb =  inf;
    upb = -inf ;


    for j = 1:length(ygrid)
        pcws(j) = legendre_conjugate(xgrid, values(:,j), -inf, inf, pcws(j));
        lwb = min(lwb, pcws(j).grid(1));
        upb = max(upb, pcws(j).grid(length(pcws(j).grid)));
    end 

%     !$omp single
    
    xdgrid = build_uniform_dual_grid(xdgrid, lwb, upb, length(xgrid), parameters);

    for j = 1:length(ygrid)
        new_values(1:length(xdgrid),j) = - pcws(j).evaluate(xdgrid);
    end
    
    
    lwb =  inf;
    upb = -inf;

%     !$omp do                       &
%     !$omp reduction(max : upb)     &
%     !$omp reduction(min : lwb)
    for j = 1:length(xdgrid)
        pcws(j) = legendre_conjugate(ygrid, new_values(j,:), -inf, inf, pcws(j));
        lwb = min(lwb, pcws(j).grid(1));
        upb = max(upb, pcws(j).grid(length(pcws(j).grid)));
    end
%     !$omp end do

    ydgrid= build_uniform_dual_grid(ydgrid, lwb, upb, length(ygrid), parameters);

    for j = 1:length(xdgrid)
        new_values(j,1:length(ydgrid)) = pcws(j).evaluate(ydgrid);
    end

    for j = 1:length(ydgrid)
        pcw= legendre_conjugate(xdgrid, new_values(1:length(xdgrid),j), -inf, inf, pcw);
        new_values(:,j) = - pcw.evaluate(xgrid);
    end
    
    for j = 1:length(xgrid)
        pcw = legendre_conjugate(ydgrid, new_values(j,1:length(ydgrid)), -inf, inf, pcw);
        new_values(j,:) = pcw.evaluate(ygrid);
    end 
    
end

