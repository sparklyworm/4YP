% week 4
clc,clear
%% true data
true_a = 0.3;
true_mu = 5 * 10 ^-7;
true_n = 6; % number of genes 
true_effect_size = 0.5;
%%
gene_to_fitness =  assignFitness2(true_n, true_a, true_effect_size);    % maps genotype to fitness value
%%
plotFitnessLandscape(gene_to_fitness, true)
%% generate a few replicate experiments
N = 10^8;
numGen = 600;
selective_pressure = 1;
numExp = 3;
experiments = cell(numExp, 1);
figure
for i = 1:numExp
    [true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
    experiments{i} = true_fitness;
    hold on
    plot(1:numGen, true_fitness)
end
title_text = sprintf("mutation rate = %.2g, n = %i, effect size = %.2g", true_mu, true_n, true_effect_size)
title(title_text)
%% get distribution for each experiment
numExp = length(experiments);
mu_distributions = cell(numExp, 1);
n_distributions = cell(numExp, 1);
es_distributions = cell(numExp, 1);
numSamples = 100;
mu_prior = makedist('Loguniform','Lower',10^-9,'Upper',10^-4);
n_domain =1:1:10;
%%
tic
for j = 1:numExp
    fprintf("====== Experiment %i ======\n", j)
    true_fitness = experiments{j};
    true_t_01 = computeTimeToPercentageMax(true_fitness, 0.01);
    true_t_99 = computeTimeToPercentageMax(true_fitness, 0.99);
    true_t_01_to_99 = true_t_99 - true_t_01;
    inflectionPoints =  findInflectionPoints2ndOrder(true_fitness);
    mu_samples = zeros(numSamples, 1);
    n_samples = zeros(numSamples, 1);
    es_samples = zeros(numSamples, 1);
    total_samples = 0;
    % 1. compare time_1% to time_99%
    % 2. compare n (number of levels)

    t_threshold = 30;
    n_threshold = 1;
    f_threshold = 1;
    for i = 1:numSamples
        
        t_diff  = t_threshold + 1;
        n_diff = n_threshold + 1;
        f_diff = f_threshold + 1;

        while ~(t_diff <= t_threshold && n_diff <= n_threshold && f_diff <= f_threshold)
            mu_s = random(mu_prior,1);
            n_s = randsample(n_domain,1);
            es_s = rand(1);
            sim_fitness = simulator2(mu_s, n_s, true_a, es_s, false);
            t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
            t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
            t_01_to_99 = t_99 - t_01;
            inflpts =  findInflectionPoints2ndOrder(sim_fitness);
            num_levels = round((length(inflpts) - 1)/2);
    
            t_diff = abs(t_01_to_99 - true_t_01_to_99);
            n_diff = abs(num_levels - true_n);
            f_diff = sqrt(sum((sim_fitness(t_01:t_99) - true_fitness(true_t_01:true_t_01 + t_01_to_99)).^2)/t_01_to_99);
        end
        mu_samples(i) = mu_s;
        n_samples(i) = n_s;
        es_samples(i) = es_s;
        % disp progress
        if mod(i, numSamples/10) ==0
            fprintf("sample %i \n", i)
        end
    end
    
    % record data
    mu_distributions{j} = mu_samples;
    n_distributions{j} = n_samples;
    es_distributions{j} = es_samples;
    
    % plot
    data = [mu_samples, n_samples];
    [~,mu_edges] = histcounts(log10(mu_samples), 15);
    [~,n_edges] = histcounts(n_samples, 10);
    figure
    hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
    set(gca, 'xscale','log')
    xlabel("Mutation rate, mu")
    ylabel("Number of Loci, n")
    title_text = sprintf("Joint Distribution of mu and n for experiment %i. True mu = %.2g, true n =  %i",j, true_mu, true_n);
    title(title_text)
    colorbar
end


total_time = toc
%save("abc_3_data_simulator_3", "experiments", "gene_to_fitness", "mu_distributions", "n_distributions", "true_mu", "true_a","total_time")

%% combine all experiments
numSamples = length(mu_distributions{1});
mu_distribution = zeros(numSamples*numExp, 1);
n_distribution = zeros(numSamples*numExp, 1);
es_distribution = zeros(numSamples*numExp, 1);
for j = 1:numExp
    start = (j -1)*numSamples + 1;
    finish = start + numSamples -1;

    mu_distribution(start:finish) = mu_distributions{j};
    n_distribution(start:finish) = n_distributions{j};
    es_distribution(start:finish) = es_distributions{j};
end
%% mu and n 
data = [mu_distribution, n_distribution];
figure
hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
set(gca, 'xscale','log')
colorbar
xlabel("Mutation rate, mu")
ylabel("Number of Loci, n")
title_text = sprintf("Full Joint Distribution of mu and n. True mu = %.2g, true n =  %i", true_mu, true_n);
title(title_text)
%% es and n
data = [es_distribution, n_distribution];
es_samples = es_distributions{1};
[~,es_edges] = histcounts(es_samples, 10);
figure
hist3(data,'CdataMode','auto', 'Edges', {es_edges, n_edges})
colorbar
xlabel("Effect Size")
ylabel("Number of Loci, n")
title_text = sprintf("Full Joint Distribution of n and effect size. True n = %i, true effect size =  %.2g", true_n, true_effect_size);
title(title_text)

%% mu and es
data = [mu_distribution, es_distribution];
es_samples = es_distributions{1};
[~,es_edges] = histcounts(es_samples, 10);
figure
hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, es_edges})
set(gca, 'xscale','log')
colorbar
xlabel("Mutation rate, mu")
ylabel("Effect Size")
title_text = sprintf("Full Joint Distribution of mu and effect size. True mu = %.2g, true effect size =  %.2g", true_mu, true_effect_size);
title(title_text)


%% plot experiments and inflection points
for j = 1:numExp
    fitness = experiments{j};
    inflectionPoints =  findInflectionPoints2ndOrder(fitness);
    figure
    plot(1:numGen, fitness)
    hold on
    plot(inflectionPoints, fitness(inflectionPoints), 'o')
end

%% plot mu & n joint
for j = 1:numExp
    mu_samples = mu_distributions{j};
    n_samples = n_distributions{j};
    
    data = [mu_samples, n_samples];
    figure
    hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
    set(gca, 'xscale','log')
    colorbar
    xlabel("Mutation rate, mu")
    ylabel("Number of Loci, n")
    title_text = sprintf("Joint Distribution of mu and n, experiment %i", j);
    title(title_text)
end

%% plot mu and es joint
for j = 1:numExp
    mu_samples = mu_distributions{j};
    es_samples = es_distributions{j};

    [~,es_edges] = histcounts(es_samples, 10);
    data = [mu_samples, es_samples];
    figure
    hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, es_edges})
    set(gca, 'xscale','log')
    colorbar
    xlabel("Mutation rate, mu")
    ylabel("Effect Size")
    title_text = sprintf("Joint Distribution of mu and effect size, experiment %i", j);
    title(title_text)
end

%% plot n and es joint
for j = 1:numExp
    n_samples = n_distributions{j};
    es_samples = es_distributions{j};
    [~,n_edges] = histcounts(n_samples);
    [~,es_edges] = histcounts(es_samples, 10);
    data = [es_samples, n_samples];
    figure
    hist3(data,'CdataMode','auto', 'Edges', {es_edges, n_edges})
    colorbar
    xlabel("Effet Size")
    ylabel("Number of Loci, n")
    title_text = sprintf("Joint Distribution of effect size and n, experiment %i", j);
    title(title_text)
end


%%
figure
for j = 1:numExp
    fitness = experiments{j};
    hold on
    plot(1:numGen, fitness)
end

%%
[~,mu_edges] = histcounts(log10(mu_distribution), 15);
figure
histogram(mu_samples,10.^mu_edges)
set(gca, 'xscale','log')
%%
figure
histogram(n_distribution)
