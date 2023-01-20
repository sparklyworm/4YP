function plateauMeans = plotHistogram(plateauSimMat)
    
    [~, col] = size(plateauSimMat);
    binWidth = 1;
    plateauMeans = zeros(col, 1);
    
    for c = 1:col
        X = plateauSimMat(:, c);
        plateauMeans(c) = mean(X);
        
        if c == 1
            binWidth = 1;
        elseif c == col
            binWidth = 10;
        else
            binWidth = 5;
        end
     
        [N,edges,bins] = histcounts(X,'BinWidth', binWidth);
        %edges = edges + 0.5; % shift edges 
        
        figure
        histogram(X,edges)
        title(sprintf("Histogram of Time to Plateau %d", c));
    
    end


end