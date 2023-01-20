function gene_to_fitness =  assignFitness(n, a)

    % returns dictionary mapping genotype string to fitness value 
    % "1011" -> 0.12983
    % fitness value is proportional to number of mutations + randomness
    % (RMF model)

    % n = number of gene in 1 genotype (eg. if n = 4, genotype = 1011 )
    % a = parameter to tune size of randomness

    % r_a = 1 + s_a
    % wildtype fitness = 1 (lowest fitness) 
    % selection coefficient s_a
    % fitness of mutated r_a

    [genotypes_cell, ~] = generateGenotypes(n);

    gene_to_fitness = dictionary(string.empty, []);
    for i = 0: n
        for j = 1:2^n
            gene = genotypes_cell{j};
            geneStr = convertCharsToStrings(gene);
            if length(find(gene == '1')) == i
                if i == 0
                    s_a = 0; % ensure wildtype 0000 has fitness 1
                else
                    r = -1 + (1-(-1))* rand(1); % r is random number between -1 and 1
                    s_a = i/n * exp(a*r);   % selection coefficient is proportional to number of mutations + randomness
                end  
                gene_to_fitness(geneStr) = s_a + 1;
            end
        end
    end

end