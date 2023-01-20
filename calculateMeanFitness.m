function r_mean = calculateMeanFitness(genotype_count, gene_to_fitness)
   

    genes = genotype_count.keys();
    total_population = genotype_count.values();
    N = sum(total_population);

    num = length(genes);
    
    r_mean = 0;
    for i = 1:num
        gene = genes(i);
        r_mean = r_mean + genotype_count(gene)/N * gene_to_fitness(gene);
        
    end


end

