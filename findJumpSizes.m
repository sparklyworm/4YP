function jumps = findJumpSizes(inflectionPoints, y)

    jumps = zeros(1, length(inflectionPoints)/2);
    pl_start = 1; 
    for i = 1:length(inflectionPoints)/2
        jump_start = y(inflectionPoints(pl_start));
        jump_end = y(inflectionPoints(pl_start+1));
        jumps(i) =jump_end - jump_start;
        pl_start = pl_start + 2;
    end   
end