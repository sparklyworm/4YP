function plotFitnessLandscape(gene_to_fitness, n, forwardMutation)
    
    if nargin < 3
        forwardMutation = false;
    end

    figure
    maxFitness = max(values(gene_to_fitness));
    genotypes_cell = convertStringsToChars(gene_to_fitness.keys());
    %%
    x = 0;
    y = 0;
    gene_xcoord = dictionary(string.empty, []);
    gene_ycoord = dictionary(string.empty, []);
    
    cmap = colormap;
    for i = 0: n
        x = 0;
        for j = 1:2^n
            gene = genotypes_cell{j};
            geneStr = convertCharsToStrings(gene);
            if length(find(gene == '1')) == i
                
                fitness = gene_to_fitness(geneStr);
                    
                % PLOT 
                scale_factor = fitness/maxFitness;
                hold on 
                plotGenotype(x, y, gene, scale_factor, cmap);
    
                % save coordinates
                gene_xcoord(geneStr) = x;
                gene_ycoord(geneStr) = y;
                x = x + 2;
            end
        end
        y = y -1;
    end
    colorbar

     
    xlim([-1, 11])
    ylim([-4.5 0.5])
    
    %%
    % arrows between genotypes of 1 Hamming distance within each other
    % arrows from low fitness genotype to high fitness genotype
   
    visited = {};
    lwMax = 5; % maximum linewidth; linewidth varies with difference between fitness values
    for i = 1:2^n
        gene1 = genotypes_cell{i};
        for j = 1:2^n
            gene2 = genotypes_cell{j};
            if any(strcmp(visited, gene2))
                continue
            end
            % if genotype within 1 Hamming distance
            if sum(gene1 ~= gene2) == 1
                gene1_str = convertCharsToStrings(gene1);
                gene2_str = convertCharsToStrings(gene2);
                lw = lwMax * abs(gene_to_fitness(gene1_str) - gene_to_fitness(gene2_str))/maxFitness + 1;
                
                if forwardMutation
                    if length(find(gene2=='1')) > length(find(gene1=='1')) && gene_to_fitness(gene2_str) > gene_to_fitness(gene1_str)
                        drawArrowFlag = true;
                    else
                        drawArrowFlag = false;
                    end
                else
                    drawArrowFlag = true;
                end

                if gene_to_fitness(gene1_str) < gene_to_fitness(gene2_str) && drawArrowFlag
                    hold on
                    drawArrow(gene_xcoord(gene1), gene_ycoord(gene1), gene_xcoord(gene2), gene_ycoord(gene2), lw);
                
                elseif gene_to_fitness(gene1_str) > gene_to_fitness(gene2_str) && drawArrowFlag
                    hold on
                    drawArrow(gene_xcoord(gene2), gene_ycoord(gene2), gene_xcoord(gene1), gene_ycoord(gene1), lw);
                end
            end
        end
        visited{end + 1} = gene1;
    end
    hold off

end
