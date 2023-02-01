function t = computeTimeToPercentageMax(y, percent)
    
    % percent: 0 - 1
    
    max_y = max(y);
    percent_y = max_y * percent;

    t = find(y>=percent_y, 1);

end