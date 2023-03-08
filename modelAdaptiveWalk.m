function [fitness, genotype_count_gen] =  modelAdaptiveWalk(fitnessLandscapeBones, N, mu, numGen, selective_pressure)
    
    % a simple model for adaptiveWalk, where we dunno the underlying
    % fitness landscape

    numGenotype = length(fitnessLandscapeBones);
    genotypes = [1:numGenotype]';
    gene_to_fitness = dictionary(genotypes, fitnessLandscapeBones); 
    genotype_p = dictionary(genotypes, zeros(numGenotype, 1)); 
    genotype_count = dictionary(genotypes, zeros(numGenotype, 1));
    wildtype = genotypes(1);
    genotype_count(wildtype) = N;
    genotype_p(wildtype) = 1;
    
    fitness = zeros(numGen, 1);
    genotype_count_gen = cell(numGen, 1);
    for gen = 1:numGen
        
        r_mean = calculateMeanFitness(genotype_count, gene_to_fitness);
    
        fitness(gen) = r_mean;  % record for plotting
        genotype_count_gen{gen} = genotype_count; % record for plotting 

        p = updateFrequency(genotype_p, gene_to_fitness, r_mean, selective_pressure); % apply selective pressure
        count = sampleNewGeneration(N, p); % multivariate normal distribution 
        genotype_count = dictionary(genotypes, count);
        genotype_count_new = genotype_count;
        for i = 1: length(genotypes)
            
            if i == length(genotypes) % the last fitness value has nowhere to go alrd
                 break
            end
    
            g = genotypes(i);
    
            if genotype_count(g) == 0
                continue
            end
            
            mutate_g = genotypes(i + 1); % can only mutate to the next one
            lambda = mu*genotype_count(g);
            new_gene = mutate_g;
            new_ind = poissrnd(lambda); % how many new individuals we get that mutate from g to new_gene
            
            genotype_count_new(new_gene) = genotype_count_new(new_gene) + new_ind;
            genotype_count_new(g) = genotype_count_new(g) - new_ind;
            
    
        end
        genotype_count = genotype_count_new;
        p = genotype_count_new.values()/sum(genotype_count_new.values());
        genotype_p = dictionary(genotypes, p);
        
    end
end