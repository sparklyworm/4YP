function conditionalDist = findConditional(plateau, col, plateauSimMat)
    
    % plateau: value length of previous plateau
    % col: which "step" on the "staircase" we were on for the previous plateau
    
    rows = find(plateauSimMat(:,col) == plateau);
    conditionalDist = plateauSimMat(rows, :);

end