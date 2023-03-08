% week 7
% added wiggle
% 
% week 6
% plotting objective function
% only 2 parameters: mu and n
% for each experiment
% get value of obj for each pair of (mu, n)
clc,clear
%% true data
true_a = 0.3;
true_mu = 1 * 10^-8;
true_n = 4;

gene_to_fitness =  assignFitness(true_n, true_a); 
plotFitnessLandscape(gene_to_fitness, true)
%% generate a few replicate experiments
N = 10^8;
numGen = 600;
selective_pressure = 1;
numExp = 10 ;
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
%% mu and n values we wish to get values of obj 
mu_array = logspace(-9, -4, 15);
n_array = linspace(1, 10, 10);
obj_matrix = zeros(length(n_array), length(mu_array));
%%
exp_numlevels = zeros(numExp, 1);
exp_t_01_99 = zeros(numExp, 1);
exp_wiggles =  zeros(numExp, 1);
for j = 1:numExp
    true_fitness = experiments{j};
    true_t_01 = computeTimeToPercentageMax(true_fitness, 0.01);
    true_t_99 = computeTimeToPercentageMax(true_fitness, 0.99);
    true_t_01_to_99 = true_t_99 - true_t_01;
    true_inflpts =  findInflectionPoints2ndOrder(true_fitness);
    true_num_levels = round((length(true_inflpts) - 1)/2);

    true_t_20 = computeTimeToPercentageMax(true_fitness, 0.2);
    true_t_80 = computeTimeToPercentageMax(true_fitness, 0.8);
    true_wiggle = measureWigglyness(true_fitness, true_t_20, true_t_80);

    exp_wiggles(j) = true_wiggle;
    exp_numlevels(j) = true_num_levels;
    exp_t_01_99(j) = true_t_01_to_99;

end

%%
tic
numSamples = 10;
std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
std_n = 0.8;
std_wig = 2;
sim_a = 0;
for m = 1:length(mu_array)
    disp(m)
    mu = mu_array(m);
    for n_i = 1:length(n_array)
        n = n_array(n_i);
        obj_exp = zeros(numExp, 1);
        for j = 1:numExp
            obj_samples = zeros(numSamples, 1);
            true_t_01_to_99 = exp_t_01_99(j);
            true_num_levels = exp_numlevels(j);
            true_wiggle = exp_wiggles(j);
            for i = 1:numSamples
                sim_fitness = simulator(mu, n, sim_a, false);
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
                obj_samples(i) = obj;
            end
            avg_obj = mean(obj_samples);
            obj_exp(j) = avg_obj;
        end
        obj_matrix(n_i, m) = mean(obj_exp);
    end
end


total_time = toc
disp_obj_matrix = - obj_matrix;
% signal end
load train
sound(y,Fs)


%%

mu_edges = [10^-9 10^-4];
n_edges = [1 10];
title_text = sprintf("Objective function value for mu = %.2g, n = %i", true_mu, true_n);
figure
imagesc('XData', mu_edges, 'YData', n_edges, 'CData', disp_obj_matrix)
title(title_text)
set(gca,'YDir','normal')
set(gca, 'xscale','log')
ylim([0.5 10.5])
xlim([0.5*10^-9 1.5*10^-4])
