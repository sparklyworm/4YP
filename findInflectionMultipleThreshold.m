function finalInflection = findInflectionMultipleThreshold(y, thresholds)
    
    allInflections = cell(1, length(thresholds));
    for i = 1:length(thresholds)
        [inflectionPoints,~] = findInflectionPoints(y, thresholds(i));
        allInflections{i} = inflectionPoints;
    end
    
    finalInflection = allInflections{1};
    for i = 2:length(thresholds)
        finalInflection = combineInflections(finalInflection, allInflections{i});
    end
    
end