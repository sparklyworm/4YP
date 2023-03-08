% week 5
% only 2 parameters: mu and n
clc,clear
%% true data
true_a = 0;
true_mu = 5 * 10 ^-8;
true_n = 4; % number of genes 

%%
gene_to_fitness =  assignFitness(true_n, true_a);    % maps genotype to fitness value
%%
plotFitnessLandscape(gene_to_fitness, true)
%% generate a few replicate experiments
N = 10^8;
numGen = 600;
selective_pressure = 1;
numExp = 1;
experiments = cell(numExp, 1);
figure
for i = 1:numExp
    [true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
    experiments{i} = true_fitness;
    hold on
    plot(1:numGen, true_fitness)
end
title_text = sprintf("mutation rate = %.2g, n = %i", true_mu, true_n);
title(title_text)
%% get distribution for each experiment
numExp = length(experiments);
mu_distributions = cell(numExp, 1);
n_distributions = cell(numExp, 1);

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
    true_inflpts =  findInflectionPoints2ndOrder(true_fitness);
    true_num_levels = round((length(true_inflpts) - 1)/2);
    mu_samples = zeros(numSamples, 1);
    n_samples = zeros(numSamples, 1);
   
    total_samples = 0;
    % 1. compare time_1% to time_99%
    % 2. compare n (number of levels)

    t_threshold = 30;
    n_threshold = 1;
    f_threshold = 0.5;
    obj_threshold = 1.8 ; %2.2 for 3 param, 1.2 for 2 params
    std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
    std_n = 0.8;

    for i = 1:numSamples
        
        t_diff  = t_threshold + 1;
        n_diff = n_threshold + 1;
        f_diff = f_threshold + 1;
        obj = obj_threshold + 1;

        while obj > obj_threshold
            mu_s = random(mu_prior,1);
            n_s = randsample(n_domain,1);
           
            sim_fitness = simulator(mu_s, n_s, true_a, false);
            t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
            t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
            t_01_to_99 = t_99 - t_01;
            inflpts =  findInflectionPoints2ndOrder(sim_fitness);
            num_levels = round((length(inflpts) - 1)/2);
    
            t_diff = abs(t_01_to_99 - true_t_01_to_99);
            n_diff = abs(num_levels - true_num_levels);
            f_diff = sqrt(sum((sim_fitness(t_01:t_99) - true_fitness(true_t_01:true_t_01 + t_01_to_99)).^2)/t_01_to_99);
            %obj = 1/std_n * n_diff + 1/std_time * t_diff + 1/f_threshold * f_diff;
            obj = 1/std_n * n_diff + 1/std_time * t_diff;
        end
        mu_samples(i) = mu_s;
        n_samples(i) = n_s;
        
        % disp progress
        if mod(i, numSamples/10) ==0
            fprintf("sample %i \n", i)
        end
    end
    
    % record data
    mu_distributions{j} = mu_samples;
    n_distributions{j} = n_samples;
   
    
    % plot
    %{
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
    %}
end


total_time = toc
%save("abc_3_data_simulator_3", "experiments", "gene_to_fitness", "mu_distributions", "n_distributions", "true_mu", "true_a","total_time")

%% combine all experiments
numSamples = length(mu_distributions{1});
mu_distribution = zeros(numSamples*numExp, 1);
n_distribution = zeros(numSamples*numExp, 1);

for j = 1:numExp
    start = (j -1)*numSamples + 1;
    finish = start + numSamples -1;

    mu_distribution(start:finish) = mu_distributions{j};
    n_distribution(start:finish) = n_distributions{j};
   
end
%% mu and n 
[~,mu_edges] = histcounts(log10(mu_distribution), 15);
[~,n_edges] = histcounts(n_samples);
data = [mu_distribution, n_distribution];
figure
hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
set(gca, 'xscale','log')
colorbar
xlabel("Mutation rate, mu")
ylabel("Number of Loci, n")
title_text = sprintf("Full Joint Distribution of mu and n. True mu = %.2g, true n =  %i", true_mu, true_n);
title(title_text)

% signal end
load train
sound(y,Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
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
histogram(mu_distribution,10.^mu_edges)
set(gca, 'xscale','log')
%}