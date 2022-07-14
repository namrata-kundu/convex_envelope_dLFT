function values = alternating_factorization_doubleDLFT(xgrid, ygrid, values)
%ALTERNATING_FACTORIZATION_DOUBLEDLFT Summary of this function goes here
%   Detailed explanation goes here

    %Compute discrete LFT along each column
    lwb = inf;
    upb = -inf;
    pcws = zeros(1,length(ygrid));
    
    for j=1:length(ygrid)
        pcws(j) = fast_dlft(xgrid, values(:,j));
        lwb = min(lwb, pcws(j).grid(1));
        grid_size = length(pcws(j).grid);
        upb = max(upb, pcws(j).grid(grid_size(2)));
    end
    
    %build the dual grid and evaluate each of the discrete LFTs on it
    cgrid = linspace(lwb , upb , length(xgrid));
    for j=1:length(ygrid)
        values(:,j) = -pcws(j)(cgrid)
    end
    
    %compute the convex envelope of each row
    for i=1:length(cgrid)
        pcw = convexenvelope_1d(ygrid , values(i,:));
        values(i,:) = -pcw(ygrid);
    end
    
    %compute the discrete LFT along each column
    for j=1:length(ygrid)
        pcw = fast_dlft(cgrid , values(:,j));
        values(:,j) = pcw(xgrid);
    end

end

