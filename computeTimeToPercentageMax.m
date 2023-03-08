function t = computeTimeToPercentageMax(y, percent)
    
    % percent: 0 - 1
    
    max_y = max(y);
    min_y = min(y);
    range_y = max_y - min_y;

    percent_y = range_y * percent;

    t = find(y>=percent_y + min_y, 1);

end