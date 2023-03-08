function [fitness, gene_to_fitness, genotype_count_gen] = simulator(mu, n, a, realFlag)
    % a = random tune param
    % mu = mutation rate
    % n = number of loci
    

    % constant parameters
    N = 10^8; % population size
    numGen = 600;
    selective_pressure = 1;
    
    if realFlag
        gene_to_fitness =  assignFitness(n, a);    % maps genotype to fitness value
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
    else
        gene_to_fitness = constructSimpleFitnessLandscape(n, a); 
        [fitness, genotype_count_gen] = modelAdaptiveWalk(gene_to_fitness.values(), N, mu, numGen, selective_pressure);
    end
    
    %
    
end