function drawArrow(x1, y1, x2, y2, lineWidth)
    % plot arrow from (x1, y1) to (x2, y2)
    dx = x2 - x1;
    dy = y2 -y1;
    lineScale = 0.95;
    dx = dx * lineScale;
    dy = dy * lineScale;
    grey = 192/256;
    q = quiver(x1,y1, dx, dy,0, 'Marker', '*','LineWidth', lineWidth, 'Color', [grey, grey, grey]);
 
end