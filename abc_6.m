% week 7
% only 2 parameters: mu and n
clc,clear
%% true data
true_a = 0.3;
true_mu = 1 * 10 ^-8;
true_n = 4; % number of genes 

%%
gene_to_fitness =  assignFitness(true_n, true_a);    % maps genotype to fitness value
%%
plotFitnessLandscape(gene_to_fitness, true)
%% generate a few replicate experiments
N = 10^8;
numGen = 600;
selective_pressure = 1;
numExp = 5;
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
obj_threshold = 2.2 ; %2.2 for 3 param, 1.2 for 2 params
std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
std_n = 0.8;
std_wig = 2.5;
sim_a = 0;
tic
for j = 1:numExp
    fprintf("====== Experiment %i ======\n", j)
    true_fitness = experiments{j};
    true_t_01 = computeTimeToPercentageMax(true_fitness, 0.01);
    true_t_99 = computeTimeToPercentageMax(true_fitness, 0.99);
    true_t_01_to_99 = true_t_99 - true_t_01;
    true_inflpts =  findInflectionPoints2ndOrder(true_fitness);
    true_num_levels = round((length(true_inflpts) - 1)/2);

    true_t_20 = computeTimeToPercentageMax(true_fitness, 0.2);
    true_t_80 = computeTimeToPercentageMax(true_fitness, 0.8);
    true_wiggle = measureWigglyness(true_fitness, true_t_20, true_t_80);

    mu_samples = zeros(numSamples, 1);
    n_samples = zeros(numSamples, 1);

    for i = 1:numSamples
       
        obj = obj_threshold + 1;

        while obj > obj_threshold
            mu_s = random(mu_prior,1);
            n_s = randsample(n_domain,1);
           
            sim_fitness = simulator(mu_s, n_s, sim_a, false);
            t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
            t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
            t_01_to_99 = t_99 - t_01;
            inflpts =  findInflectionPoints2ndOrder(sim_fitness);
            num_levels = round((length(inflpts) - 1)/2);
            
            t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
            t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
            wiggle = measureWigglyness(sim_fitness, t_20, t_80);

            wig_diff = abs(wiggle - true_wiggle);
            t_diff = abs(t_01_to_99 - true_t_01_to_99);
            n_diff = abs(num_levels - true_num_levels);
            obj = 1/std_n * n_diff + 1/std_time * t_diff + 1/std_wig * wig_diff;
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