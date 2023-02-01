function plotFitnessLandscapeGraph(gene_to_fitness)

    gene2genes = hammingMapFitness(gene_to_fitness);
    allPaths = findAllPaths(gene2genes);
    allPathFitness = pathToFitness(allPaths, gene_to_fitness);
    numPaths = length(allPaths);
    
    [global_optima_fitness, ind] = max(gene_to_fitness.values());
    genotypes = gene_to_fitness.keys();
    global_optima_gene = genotypes(ind);
    max_length = length(find(convertStringsToChars(global_optima_gene) == '1'));
   
    figure
    for i = 1:numPaths
        path = allPaths{i};
        pathFitness = allPathFitness{i};
        n = length(path);
        if n == max_length + 1
            lineColor = '#41B7FF';
        else
            lineColor = '#C5C5C5';
        end
        hold on
        plot(1:n, allPathFitness{i}, '-o', 'MarkerSize',8, 'LineWidth',1.5, 'Color', lineColor)
        for j = 1:length(path)
            text(j, pathFitness(j), "  " + path(j))
        end
    end
    
    hold on
    plot(max_length+1, global_optima_fitness, 'Marker','pentagram', 'MarkerSize', 15, 'MarkerFaceColor','red', 'MarkerEdgeColor','red')
    
end
