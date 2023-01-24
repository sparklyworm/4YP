function gene2genes = hammingMapFitness(gene_to_fitness)
    
    % gene_to_fitness: dictionary mapping genotype to 1 fitness value
    % gene2genes: dictionary mapping genotype to genotypes within 
    % 1 hamming distance and of higher fitness

    % default assume NO backwards mutations -
    % genotype can only ACQUIRE mutations, not lose them 

  

    genotypes = gene_to_fitness.keys();
    num = length(genotypes);
    gene2genes = dictionary(genotypes, cell(num, 1));
    genotypes_groups = groupGenotypes(genotypes);
    
    for i = 1:length(genotypes_groups)-1
        gene_group = genotypes_groups{i};
        next_group = genotypes_groups{i+1};
        for j = 1:length(gene_group)
            gene1 = gene_group(j);
            for k = 1:length(next_group)
                gene2 = next_group(k);
                gene1_char = convertStringsToChars(gene1);
                gene2_char = convertStringsToChars(gene2);
                % within 1 hamming distance and has higher fitness
                if sum(gene1_char ~= gene2_char) == 1 && gene_to_fitness(gene2) > gene_to_fitness(gene1) 
                    c = gene2genes(gene1);
                    if isempty(c{1})
                        gene2genes(gene1) = {gene2};
                    elseif length(c{1}) == 1
                        newc = [c, {gene2}];
                        gene2genes(gene1) = {newc};
                    else
                        newc = [c{1}, {gene2}];
                        gene2genes(gene1) = {newc};
                    end
                end
            end
        end
    end

end


%{
visited = {};
    for i = 1:num
        gene1 = genotypes(i);
        for j = 1:num
            gene2 = genotypes(j);
            if any(strcmp(visited, gene2))
                continue
            end
            
            gene1_char = convertStringsToChars(gene1);
            gene2_char = convertStringsToChars(gene2);

            % if genotype within 1 Hamming distance 
            % and has 1 additional mutation
            if length(find(gene2_char == '1')) == length(find(gene1_char == '1')) + 1 && sum(gene1_char ~= gene2_char) == 1
                % and has higher fitness
                if gene_to_fitness(gene2) > gene_to_fitness(gene1)
                    c = gene2genes(gene1);
                    if isempty(c{1})
                        gene2genes(gene1) = {gene2};
                    elseif length(c{1}) == 1
                        newc = [c, {gene2}];
                        gene2genes(gene1) = {newc};
                    else
                        newc = [c{1}, {gene2}];
                        gene2genes(gene1) = {newc};
                    end
                   
                end
            end
        end
        visited{end + 1} = gene1;
    end







%}