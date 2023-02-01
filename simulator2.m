function [fitness, gene_to_fitness, genotype_count_gen] = simulator2(mu, n, a, realFlag)
    % n = number of loci
    % mu = mutation rate

    % constant parameters 
    N = 10^8; % population size
    numGen = 600;
    selective_pressure = 1;
    
    if realFlag
        gene_to_fitness =  assignFitness(n, a);    % maps genotype to fitness value
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
    else
        gene_to_fitness = constructSimpleFitnessLandscape(n, a); 
        fitness = modelAdaptiveWalk(gene_to_fitness.values(), N, mu, numGen, selective_pressure);
    end
    
    %
    
end