function plateauLengths = computePlateauLength(inflectionPoints)
    
    plateauLengths = zeros(1, length(inflectionPoints)/2);
    plateauLengths(1) = inflectionPoints(1);
    pl_start = 2; 
    for i = 2:length(inflectionPoints)/2
        plateauLengths(i) = inflectionPoints(pl_start+1) - inflectionPoints(pl_start);
        pl_start = pl_start + 2;
    end   
end