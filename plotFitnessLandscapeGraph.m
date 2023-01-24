function plotFitnessLandscapeGraph(gene_to_fitness)

    gene2genes = hammingMapFitness(gene_to_fitness);
    allPaths = findAllPaths(gene2genes);
    pathFitness = pathToFitness(allPaths, gene_to_fitness);

    numPaths = length(allPaths);
    
   
    figure
    for i = 1:numPaths
        hold on
        n = length(pathFitness{i});
        plot(1:n, pathFitness{i}, '-o', 'MarkerSize',10, 'LineWidth',2)
    end
    
end
