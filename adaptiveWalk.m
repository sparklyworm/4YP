function  [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure)
    
    genotypes = gene_to_fitness.keys(); 
    numGenotype = length(genotypes);
    genotype_count = dictionary(genotypes, zeros(numGenotype, 1)); % count of individuals of each genotype
    genotype_p = dictionary(genotypes, zeros(numGenotype, 1));     % frequency """
    gene2genes = hammingMapFitness(gene_to_fitness); % maps genotype to genotypes 1 mutation step away & are of higher fitness 
    
    wildtype = genotypes(1);
    genotype_count(wildtype) = N;
    genotype_p(wildtype) = 1;
    
    fitness = zeros(numGen, 1);
    genotype_count_gen = cell(numGen, 1);
    for gen = 1:numGen
    
        r_mean = calculateMeanFitness(genotype_count, gene_to_fitness);
    
        fitness(gen) = r_mean;                    % record for plotting
        genotype_count_gen{gen} = genotype_count; % record for plotting 
        
        p = updateFrequency(genotype_p, gene_to_fitness, r_mean, selective_pressure); % apply selective pressure
        count = sampleNewGeneration(N, p); % multivariate normal distribution 
        genotype_count = dictionary(genotypes, count);
        genotype_count_new = genotype_count;
    
        for i = 1: length(genotypes)
            
            g = genotypes(i);
            if genotype_count(g) == 0
                continue
            end
    
            genes = gene2genes(g);
            % mutate_g is the array of genotypes that are 1 mutation step away from gene
            if length(genes{1}) == 1
                mutate_g = genes{1};
            else
                mutate_g = string(genes{1}); % string array ["0110", "1011"]
            end
    
            lambda = mu*genotype_count(g);
            for k = 1:length(mutate_g)
                new_gene = mutate_g(k);
                new_ind = poissrnd(lambda);
                
                genotype_count_new(new_gene) = genotype_count_new(new_gene) + new_ind;
                genotype_count_new(g) = genotype_count_new(g) - new_ind;
            end
    
        end
        genotype_count = genotype_count_new;
        p = genotype_count_new.values()/sum(genotype_count_new.values());
        genotype_p = dictionary(genotypes, p);
        
    end

end