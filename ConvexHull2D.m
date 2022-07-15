classdef ConvexHull2D
    %CONVEXHULL2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods(Static)
%         function obj = ConvexHull2D(inputArg1,inputArg2)
%             %CONVEXHULL2D Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
        
        function array = injector(array, value)
            i = bisect(array, value);
            if (i < length(array))
                if (array(i+1) - value < value - array(i)) 
                    i = i + 1;
                end
            end
            array(i) = value;
        end

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

        function dgrid = build_uniform_dual_grid(dgrid, lwb, upb, N, parameters)
            [lwb, upb] = ConvexHull2D.modify_dual_grid_bounds(lwb, upb, parameters);
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
                dgrid = ConvexHull2D.injector(dgrid, parameters('dg_inject_value'));
            end 

        end

        function result = maxval_except(array, except)
            result = max(max(array(1:except-1)), max(array(except+1:end)));
        end

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

                fixation_time = max(fixation_time, ConvexHull2D.maxval_except(  arr, except=hull_indices(j+1)-hull_indices(j)));
                j = j + 1;

            end

            asymptotic_rates = ConvexHull1D.piecewise_linear_interpolation( grid(hull_indices(1:j)), rates((hull_indices(1:j))),  ...
                                                               -inf, inf, grid                                );

        end


        function new_values = convex_hull2d(values, xgrid, ygrid, parameters)
        %CONVEX_HULL2D Summary of this function goes here
        %   Detailed explanation goes here

            if parameters('doubleDLFT')
                new_values = ConvexHull2D.convex_hull2d_doubleDLFT(values, xgrid, ygrid, parameters);
            elseif ~parameters('symmetric')
                if parameters('xy_order')
                    new_values = ConvexHull2D.convex_hull2d_asymmetric(values, xgrid, ygrid, parameters);
                else
                    new_values_transpose = ConvexHull2D.convex_hull2d_asymmetric(values', ygrid, xgrid, parameters);
                    new_values = new_values_transpose';
                end
            else
                new_values = ConvexHull2D.convex_hull2d_asymmetric(values, xgrid, ygrid, parameters);
                sym_values_transpose = ConvexHull2D.convex_hull2d_asymmetric(values', ygrid, xgrid, parameters);
                sym_values = sym_values_transpose';
                new_values = max(new_values, sym_values);
            end

        end

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
                pcws(j) = ConvexHull1D.legendre_conjugate(xgrid, values(:,j), -inf, inf, pcws(j));
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
                [left_slopes, tfix]= ConvexHull2D.asymptotic_convex_hull(ygrid, border_values, [pcws.left_slope], left_slopes, tfix);
        %         ! we did not change sign for the left_slope since that is already taken in account
        %         ! by the fact that "time" is inverted
                lwb = lwb - tfix;

                for j = 1:length(ygrid) 
                    border_values(j) = - pcws(j).values(length(pcws(j).values)) + pcws(j).right_slope * (pcws(j).grid(length(pcws(j).values)) - upb);
                end
                [right_slopes, tfix] = ConvexHull2D.asymptotic_convex_hull(ygrid, border_values, -pcws(:).right_slope, right_slopes, tfix);
                upb = upb + tfix;

        %         deallocate(border_values)
            end
        %     !$omp end single

        % % %     dual_grid_choice: 
            if (parameters('dual_grid_type') == 0)
                ugrid = ConvexHull2D.build_uniform_dual_grid(ugrid, lwb, upb, length(xgrid), parameters);
                new_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, new_values,ugrid,left_slopes,right_slopes);
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
                     new_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, new_values,pcws(1).grid,left_slopes,right_slopes);
                     if length(left_slopes)
                         left_slopes = [];
                         right_slopes = [];
        %              if (allocated(left_slopes)) deallocate(left_slopes,right_slopes)
                     end

                     for j = 1 + step:step: length(ygrid)
                         temp_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, temp_values,pcws(j).grid,left_slopes,right_slopes);
                         new_values = max(new_values, temp_values);
                     end

                     if (j < size(ygrid) + step) % ! the interval was not divisible exactly
                         temp_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, temp_values,pcws(size(ygrid)).grid,left_slopes,right_slopes);
                         new_values = max(new_values, temp_values);
                     end 
                 else
                     new_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, new_values, pcws(min(size(ygrid),1 + (parameters('dg_only_slice') - 1) * step)).grid, left_slopes,right_slopes);
                 end 
             elseif (parameters('dual_grid_type') == 2) 
                 temp_values = zeros(length(xgrid), length(ygrid));
                 ugrid = zeros(1, length(xgrid));
                 [lwb, upb] = ConvexHull2D.modify_dual_grid_bounds(lwb, upb, parameters);

                 %%%%%%%subgrids
                 for i = 1: parameters('dg_slices')
                     for j = 1:length(ugrid) 
                         ugrid(j) = lwb + ( (i-1) * size(ugrid) + j - 1 ) * (upb - lwb) / ( size(ugrid) * parameters('dg_slices') - 1.0 );
                     end

                     if (i==1) 
                         new_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, new_values,ugrid,left_slopes,right_slopes);
                         if (length(left_slopes)) 
                             left_slopes=[];
                             right_slopes =[];
                         end
                     else
                         temp_values = ConvexHull2D.stage2(pcws, xgrid, ygrid, temp_values,ugrid,left_slopes,right_slopes);
                         new_values = max(new_values, temp_values);
                     end
                 end         
            end 
        % % %     dual_grid_choice       
        end

        function out_values = stage2(pcws, xgrid, ygrid, out_values,ugrid,left_slopes,right_slopes)

            pcw = PiecewiseLinear_t();

            for i = 1:length(ygrid)
                out_values(1:length(ugrid),i) = - pcws(i).evaluate(ugrid);
            end

            pcw.left_slope  = -inf;
            pcw.right_slope =  inf; %! needed so we can evaluate correctly all cases
            for i = 1:length(ugrid)
                [pcw.grid, pcw.values] = ConvexHull1D.convex_hull(ygrid, out_values(i,:), -Inf(), Inf(), pcw.grid, pcw.values);
                out_values(i,:) = - pcw.evaluate(ygrid);
            end 

            if (length(left_slopes)) 
                for i = 1:length(ygrid)
                    pcw = ConvexHull1D.legendre_conjugate(ugrid, out_values(1:length(ugrid),i), left_slopes(i), -right_slopes(i), pcw);
        %             ! we did not change sign for the left_slope since that was already done
        %             ! in the asymptotic hull since "time" is inverted on the left
                    out_values(:, i) = pcw.evaluate(xgrid);
                end
            else
                for i = 1:length(ygrid)
                    pcw = ConvexHull1D.legendre_conjugate(ugrid, out_values(1:length(ugrid),i), -inf, inf, pcw);
                    out_values(:, i) = pcw.evaluate(xgrid);
                end
            end 
        end

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
                pcws(j) = ConvexHull1D.legendre_conjugate(xgrid, values(:,j), -inf, inf, pcws(j));
                lwb = min(lwb, pcws(j).grid(1));
                upb = max(upb, pcws(j).grid(length(pcws(j).grid)));
            end 

        %     !$omp single

            xdgrid = ConvexHull2D.build_uniform_dual_grid(xdgrid, lwb, upb, length(xgrid), parameters);

            for j = 1:length(ygrid)
                new_values(1:length(xdgrid),j) = - pcws(j).evaluate(xdgrid);
            end


            lwb =  inf;
            upb = -inf;

        %     !$omp do                       &
        %     !$omp reduction(max : upb)     &
        %     !$omp reduction(min : lwb)
            for j = 1:length(xdgrid)
                pcws(j) = ConvexHull1D.legendre_conjugate(ygrid, new_values(j,:), -inf, inf, pcws(j));
                lwb = min(lwb, pcws(j).grid(1));
                upb = max(upb, pcws(j).grid(length(pcws(j).grid)));
            end
        %     !$omp end do

            ydgrid= ConvexHull2D.build_uniform_dual_grid(ydgrid, lwb, upb, length(ygrid), parameters);

            for j = 1:length(xdgrid)
                new_values(j,1:length(ydgrid)) = pcws(j).evaluate(ydgrid);
            end

            for j = 1:length(ydgrid)
                pcw= ConvexHull1D.legendre_conjugate(xdgrid, new_values(1:length(xdgrid),j), -inf, inf, pcw);
                new_values(:,j) = - pcw.evaluate(xgrid);
            end

            for j = 1:length(xgrid)
                pcw = ConvexHull1D.legendre_conjugate(ydgrid, new_values(j,1:length(ydgrid)), -inf, inf, pcw);
                new_values(j,:) = pcw.evaluate(ygrid);
            end 

        end


        
    end
end

