function total_e = measureWigglyness(y, startpt, endpt)
    % startpt, endpt: start and end index

    y_start = y(startpt);
    y_end = y(endpt);
    
    % find equation of straightline between start and end
    % y = mx + c;

    m = (y_end - y_start)/(endpt-startpt);
    c = y_end - m * endpt;
    
    % points along that line between start and end
    X = [startpt:1:endpt]';
    Y = m * X + c;
    
    % difference between curve and straight line
    error = Y - y(startpt:endpt);
    
    total_e = sum(abs(error));
end