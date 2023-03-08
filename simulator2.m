function [fitness, gene_to_fitness, genotype_count_gen] = simulator2(mu, n, a, effectSize, realFlag)
    % n = number of loci
    % mu = mutation rate

    % constant parameters 
    N = 10^8; % population size
    numGen = 600;
    selective_pressure = 1;
    
    if realFlag
        gene_to_fitness =  assignFitness2(n, a, effectSize);    % maps genotype to fitness value
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
    else
        gene_to_fitness = constructSimpleFL2(n, a, effectSize); 
        [fitness, genotype_count_gen] = modelAdaptiveWalk(gene_to_fitness.values(), N, mu, numGen, selective_pressure);
    end
    
    %
    
end