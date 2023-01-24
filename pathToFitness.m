function pathFitness = pathToFitness(allPaths, gene_to_fitness)
    numPaths = length(allPaths);
    pathFitness = cell(numPaths, 1);
    
    for i = 1:numPaths
        path = allPaths{i};
        path_len = length(path);
        fitness = zeros(path_len, 1);
        for j = 1:path_len
            gene = path(j);
            fitness(j) = gene_to_fitness(gene);
        end
        pathFitness{i} = fitness;
    end
end