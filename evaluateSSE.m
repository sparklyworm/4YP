function SSE = evaluateSSE(true_fitness, gene_to_fitness, N, mu_c, numGen, selective_pressure)
    % used in script bayesian_optimisation.m 
    
    num_samples = 10;
    samples = zeros(numGen, num_samples);

    % generate 10 samples
    for s = 1:num_samples
        [fitness, ~] = adaptiveWalk(gene_to_fitness, N, mu_c, numGen, selective_pressure);
        samples(:, s) = fitness;
    end
    avg_fitness = sum(samples, 2)/num_samples;
    SSE = sum((avg_fitness - true_fitness).^2);
    
end