function allPaths = findAllPaths(gene2genes)
    
    
    genotypes = gene2genes.keys();
    n = length(convertStringsToChars(genotypes(1))); % legnth of genotype eg. 0000 -> 4 
    genotypes_groups = groupGenotypes(genotypes);

    allPaths = {[genotypes(1)]};
    for i = 1:length(genotypes_groups)-1
        gene_group = genotypes_groups{i};
        for j = 1:length(gene_group)
            g1 = gene_group(j);
            c = gene2genes(g1);
            next_genes = c{1};
            newPaths = {};
            for k = 1:length(allPaths)
                path = allPaths{k};
                if path(end) == g1
                    for j = 1:length(next_genes)
                        newPath = [path next_genes{j}];
                        newPaths{end + 1} = newPath;
                    end
                    if ~isempty(newPaths)
                        allPaths{k} = []; % remove "old path"
                    end
                end
            end
    
            % merge allPaths and newPaths;
            for s = 1:length(allPaths)
                if isempty(allPaths{s})
                    continue
                else
                    newPaths{end + 1} = allPaths{s};
                end
            end
            allPaths = newPaths;
        end
    end
end
