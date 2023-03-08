% week 3
clc,clear
%% true data
true_a = 0.1;
true_mu = 1 * 10 ^-8;
true_n = 4; % number of genes 
%%
gene_to_fitness =  assignFitness(true_n, true_a);    % maps genotype to fitness value
plotFitnessLandscape(gene_to_fitness, true)
%% generate 10 replicate experiments
N = 10^8;
numGen = 600;
selective_pressure = 1;
numExp = 3;
experiments = cell(numExp, 1);
for i = 1:numExp
    [true_fitness, ~] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
    experiments{i} = true_fitness;
end
%%
% true_fitness = experiments{5};
% percent = 0.99;
% true_t = computeTimeToPercentageMax(true_fitness, percent);
% inflectionPoints =  findInflectionPoints2ndOrder(true_fitness);
% figure
% plot(1:numGen, true_fitness)
% hold on
% plot(true_t, true_fitness(true_t), 'o')
% hold on
% plot(inflectionPoints, true_fitness(inflectionPoints), 'go')
%% get distribution for each experiment
numExp = 3;
mu_distributions = cell(numExp, 1);
n_distributions = cell(numExp, 1);
percent = 0.99;
numSamples = 100;
mu_prior = makedist('Loguniform','Lower',10^-9,'Upper',10^-4);
n_domain =1:1:10;
%%
tic
for j = 1:numExp
    fprintf("====== Experiment %i ======\n", j)
    true_fitness = experiments{j};
    true_t = computeTimeToPercentageMax(true_fitness, percent);
    inflectionPoints =  findInflectionPoints2ndOrder(true_fitness);
    mu_samples = zeros(numSamples, 1);
    n_samples = zeros(numSamples, 1);
    total_samples = 0;
    % 1. time to 99% of max fitness
    % 2. compare fitness path within the same time frame
    % 3. compare n (number of levels)
    t_threshold = 30;
    f_threshold = 0.2;
    n_threshold = 1;
    
    for i = 1:numSamples
        
        t_diff  = t_threshold + 1;
        f_diff = f_threshold + 1;
        n_diff = n_threshold + 1;
        while ~(t_diff <= t_threshold && f_diff <= f_threshold && n_diff <= n_threshold)
            mu_s = random(mu_prior,1);
            n_s = randsample(n_domain,1);
    
            sim_fitness = simulator2(mu_s, n_s, true_a, false);
            t = computeTimeToPercentageMax(sim_fitness, 0.99);
            inflpts =  findInflectionPoints2ndOrder(sim_fitness);
            num_levels = round((length(inflpts) - 1)/2);
    
            t_diff = abs(t - true_t);
            f_diff = sqrt(sum((sim_fitness(1:t) - true_fitness(1:t)).^2)/t);
            n_diff = abs(num_levels - true_n);
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
    data = [mu_samples, n_samples];
    [~,mu_edges] = histcounts(log10(mu_samples), 15);
    [~,n_edges] = histcounts(n_samples);
    figure
    hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
    set(gca, 'xscale','log')
    colorbar
end


total_time = toc
 save("abc_3_data_simulator_3", "experiments", "gene_to_fitness", "mu_distributions", "n_distributions", "true_mu", "true_a","total_time")
%%
[~,mu_edges] = histcounts(log10(mu_samples), 15);
figure
histogram(mu_samples,10.^mu_edges)
set(gca, 'xscale','log')
%%
figure
[~,n_edges] = histcounts(n_samples);
histogram(n_samples)
%%
data = [mu_samples, n_samples];
figure
hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
set(gca, 'xscale','log')
colorbar

%%
for j = 1:numExp
    fitness = experiments{j};
    inflectionPoints =  findInflectionPoints2ndOrder(fitness);
    figure
    plot(1:numGen, fitness)
    hold on
    plot(inflectionPoints, fitness(inflectionPoints), 'o')
end

%%
for j = 1:numExp
    mu_samples = mu_distributions{j};
    n_samples = n_distributions{j};

%     figure
%     [~,mu_edges] = histcounts(log10(mu_samples), 15);
%     histogram(mu_samples,10.^mu_edges)
%     set(gca, 'xscale','log')
%     
%     figure
%     [~,n_edges] = histcounts(n_samples);
%     histogram(n_samples)

    data = [mu_samples, n_samples];
    figure
    hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, n_edges})
    set(gca, 'xscale','log')
    colorbar
end

%%
figure
for j = 1:numExp
    fitness = experiments{j};
    hold on
    plot(1:numGen, fitness)
end
