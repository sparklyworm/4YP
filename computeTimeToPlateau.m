function times = computeTimeToPlateau(inflectionPoints)
    
    times = zeros(1, length(inflectionPoints)/2);
    pl_start = 1; 
    for i = 1:length(inflectionPoints)/2
        times(i) = inflectionPoints(pl_start+1) - inflectionPoints(pl_start);
        pl_start = pl_start + 2;
    end   
end