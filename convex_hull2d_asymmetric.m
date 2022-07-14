function new_values = convex_hull2d_asymmetric(values, xgrid, ygrid, parameters)

    %Initiate required variables
    pcws(1:length(ygrid)) = PiecewiseLinear_t;
    tfix=0;
    ugrid=[];
    new_values = zeros(length(xgrid),length(ygrid));
    left_slopes = [];
    right_slopes = [];
    
    lwb =  inf;
    upb = -inf;

%     !$omp default(shared)                                         &
%     !$omp private(i, j)
%         !$omp do                       &
%         !$omp reduction(max : upb)     &
%         !$omp reduction(min : lwb)
    for j = 1:length(ygrid)
        pcws(j) = legendre_conjugate(xgrid, values(:,j), -inf, inf, pcws(j));
        lwb = min(lwb, pcws(j).grid(1));
        upb = max(upb, pcws(j).grid(length(pcws(j).grid)));
    end 
%     !$omp end do
% 
%     !$omp single
    if (parameters('dg_ext_slopes')) 
        left_slopes = zeros(1, length(ygrid));
        right_slopes = zeros(1, length(ygrid));
        border_values = zeros(1, length(ygrid));

        for j = 1:length(ygrid) 
            border_values(j) = - pcws(j).values(1) + pcws(j).left_slope * (pcws(j).grid(1) - lwb);
        end
        [left_slopes, tfix]= asymptotic_convex_hull(ygrid, border_values, pcws(:).left_slope, left_slopes, tfix);
%         ! we did not change sign for the left_slope since that is already taken in account
%         ! by the fact that "time" is inverted
        lwb = lwb - tfix;

        for j = 1:length(ygrid) 
            border_values(j) = - pcws(j).values(length(pcws(j).values)) + pcws(j).right_slope * (pcws(j).grid(length(pcws(j).values)) - upb);
        end
        [right_slopes, tfix] = asymptotic_convex_hull(ygrid, border_values, -pcws(:).right_slope, right_slopes, tfix);
        upb = upb + tfix;

%         deallocate(border_values)
    end
%     !$omp end single

% % %     dual_grid_choice: 
    if (parameters('dual_grid_type') == 0)
        ugrid = build_uniform_dual_grid(ugrid, lwb, upb, length(xgrid), parameters);
        new_values = stage2(pcws, xgrid, ygrid, new_values,ugrid,left_slopes,right_slopes);
    elseif (parameters('dual_grid_type') == 1)
         if (parameters('dg_slices') > 0) 
             step = (length(ygrid) - 1) / (parameters('dg_slices') - 1);
         else
             step = - parameters('dg_slices');
         end
%          ! we could have used an associate block for step, but f2py does not like them
         if (parameters('dg_only_slice') == 0) 
             temp_values = zeros(length(xgrid), length(ygrid));             
         end

         if (parameters('dg_only_slice') == 0) 
             new_values = stage2(pcws, xgrid, ygrid, new_values,pcws(1).grid,left_slopes,right_slopes);
             if length(left_slopes)
                 left_slopes = [];
                 right_slopes = [];
%              if (allocated(left_slopes)) deallocate(left_slopes,right_slopes)
             end

             for j = 1 + step:step: length(ygrid)
                 temp_values = stage2(pcws, xgrid, ygrid, temp_values,pcws(j).grid,left_slopes,right_slopes);
                 new_values = max(new_values, temp_values);
             end

             if (j < size(ygrid) + step) % ! the interval was not divisible exactly
                 temp_values = stage2(pcws, xgrid, ygrid, temp_values,pcws(size(ygrid)).grid,left_slopes,right_slopes);
                 new_values = max(new_values, temp_values);
             end 
         else
             new_values = stage2(pcws, xgrid, ygrid, new_values, pcws(min(size(ygrid),1 + (parameters('dg_only_slice') - 1) * step)).grid, left_slopes,right_slopes);
         end 
     elseif (parameters('dual_grid_type') == 2) 
         temp_values = zeros(length(xgrid), length(ygrid));
         ugrid = zeros(1, length(xgrid));
         [lwb, upb] = modify_dual_grid_bounds(lwb, upb, parameters);

         %%%%%%%subgrids
         for i = 1: parameters('dg_slices')
             for j = 1:size(ugrid) 
                 ugrid(j) = lwb + ( (i-1) * size(ugrid) + j - 1 ) * (upb - lwb) / ( size(ugrid) * parameters('dg_slices') - 1.0 );
             end
             
             if (i==1) 
                 new_values = stage2(pcws, xgrid, ygrid, new_values,ugrid,left_slopes,right_slopes);
                 if (length(left_slopes)) 
                     left_slopes=[];
                     right_slopes =[];
                 end
             else
                 temp_values = stage2(pcws, xgrid, ygrid, temp_values,ugrid,left_slopes,right_slopes);
                 new_values = max(new_values, temp_values);
             end
         end         
    end 
% % %     dual_grid_choice       
end

