function all_genotypes = plotAdaptiveWalk(genotype_count_gen, mu, threshold, logScale)
    
    % genotype_count_gen (cell of dictionaries) 
    % == the population count of each genotype at each generation

    % logScale: boolean True or False (plot in logscale or not)
    % threshold: population count threshold for which to plot 
    % (ie. max population count must be > threshold to be plotted)

    genotypes = genotype_count_gen{1}.keys();
    n = length(convertStringsToChars(genotypes(1)));
    numGen = length(genotype_count_gen);
    all_genotypes = cell(2^n, 1);

    for i= 1:2^n
        all_genotypes{i} = zeros(numGen, 1);
    end
    
    for gen = 1:numGen
        g = genotype_count_gen{gen};
        count = g.values();
        for i = 1:2^n
            genotype = all_genotypes{i};
            genotype(gen) = count(i);
            all_genotypes{i} = genotype;
        end
    end
    
    
    figure
    labelGraph = {};
    %cmap = [hsv(2^(n-1));hsv(2^(n-1))];
    cmap = distinguishable_colors(2^n);
    for i= 1:2^n
        % only plotting genotypes that reach popoulation above 10^7
        if max(all_genotypes{i}) > threshold
            disp(max(all_genotypes{i}))
            if logScale
                plot([1:numGen],log(all_genotypes{i}), 'LineWidth', 3, 'Color', cmap(i, :))
            else
                plot([1:numGen],all_genotypes{i}, 'LineWidth', 3, 'Color', cmap(i, :))
            end
            
            labelGraph{end+1} = genotypes(i);
            hold on
        end
    end
    if logScale
        ytitle = "Population (Log scale)";
    else
        ytitle = "Population";
    end   
    ylabel(ytitle)
    xlabel("Generation")
    title(sprintf("Genotypes with population count more than %0.3g, mutation rate %0.3g", threshold, mu));
    legend(labelGraph)

end