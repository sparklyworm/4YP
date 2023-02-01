function gene_to_fitness = roughMountFuji(n, a)
    % assume global optima is genotype with maximum num of mutations 
    % eg. "1111"

    [genotypes_cell, ~] = generateGenotypes(n);
    global_optima = genotypes_cell{end};

    gene_to_fitness = dictionary(string.empty, []);
    for i = 1:2^n
        gene = genotypes_cell{i};
        geneStr = convertCharsToStrings(gene);
        d = sum(gene ~= global_optima); 
        gene_to_fitness(geneStr) = -a*d + normrnd(0, 1);
        
    end
    
    % shift fitness values so that they are all positive, and minima is 1
    fitness_minima = min(gene_to_fitness.values());
    fitness_values = gene_to_fitness.values() + abs(fitness_minima) + 1;
    genotypes = gene_to_fitness.keys();
    gene_to_fitness = dictionary(genotypes, fitness_values);

end