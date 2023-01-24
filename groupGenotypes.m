function genotypes_groups = groupGenotypes(genotypes)
   
    n = length(convertStringsToChars(genotypes(1)));
    genotypes_groups= cell(n+1, 1); % cell of genotypes groups by the number of mutations they have

    for j = 1:length(genotypes)
        gene = genotypes(j);
        gene_char = convertStringsToChars(gene);
        num_mutation = length(find(gene_char == '1'));
        genotypes_groups{num_mutation+1} = [genotypes_groups{num_mutation+1} gene];
    end

end