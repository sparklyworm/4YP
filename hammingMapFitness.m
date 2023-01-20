function gene2genes = hammingMapFitness(gene_to_fitness, backwards)
    
    % gene_to_fitness: dictionary mapping genotype to 1 fitness value
    % gene2genes: dictionary mapping genotype to genotypes within 
    % 1 hamming distance and of higher fitness

    % default assume NO backwards mutations -
    % genotype can only ACQUIRE mutations, not lose them 

    if nargin < 2
        backwards = false;
    end

    gene_arr = gene_to_fitness.keys();
    num = length(gene_arr);
    gene2genes = dictionary(gene_arr, cell(num, 1));

    visited = {};
    for i = 1:num
        gene1 = gene_arr(i);
        for j = 1:num
            gene2 = gene_arr(j);
            if any(strcmp(visited, gene2))
                continue
            end
            
            gene1_char = convertStringsToChars(gene1);
            gene2_char = convertStringsToChars(gene2);

            % if genotype within 1 Hamming distance 
            % and has 1 additional mutation
            if length(find(gene2_char == '1')) == length(find(gene1_char == '1')) + 1
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


end