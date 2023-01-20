function plotGenotype(x, y, gene, scale_factor, cmap)
        % x, y: coordinates to plot (int, int)
        % gene: string "1001"
        % scale factor = fitness/maxFitness; 
        % to scale marker size and pick colors according to fitness value
        % cmap = colormap
        markerSizeMax = 50;
        markerSize = round(markerSizeMax * scale_factor) + 5;
        markerColor = min(round(256 * scale_factor) + 1, 256);
        plot(x, y,'o','MarkerSize', markerSize, 'MarkerFaceColor', cmap(markerColor, :), 'MarkerEdgeColor', cmap(markerColor, :));
        text(x-0.3, y, gene, 'FontSize',14,'FontWeight','bold');
end