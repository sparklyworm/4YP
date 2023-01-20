           function new_P = updateFrequency(genotype_p, gene_to_fitness, r_mean, selective_pressure)
    
    % new_P is a n x 1 vector  
    % selective pressure: tunning param between 0 and 1 

    P = genotype_p.values();
    genes = genotype_p.keys();
    n = length(genes);
    
    new_P = zeros(n, 1);

    for i = 1:n
        gene = genes(i);
        p = P(i)*gene_to_fitness(gene)/r_mean;
        d = p - P(i);
        new_P(i) = P(i) + selective_pressure * d;
    end

end