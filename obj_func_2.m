%%
% Week 9
% calculate the ratio of the 5 experiments which match the sample (obj<10)
clc,clear

%% true data
true_a = 0.2;

%%
true_mu = 5 * 10^-6;
true_n = 6;
%%
gene_to_fitness =  assignFitness(true_n, true_a); 
%%
plotFitnessLandscape(gene_to_fitness, true)
%% generate a few replicate experiments
N = 10^8;
numGen = 2800;
selective_pressure = 1;
numExp = 5 ;
%%
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
exp_t_01_95 = zeros(numExp, 1);
exp_wiggles =  zeros(numExp, 1);
for j = 1:numExp
    true_fitness = experiments{j};
   [true_t_01_to_95, true_wiggle, true_num_levels] = calcSummaryStats(true_fitness);
    exp_wiggles(j) = true_wiggle;
    exp_numlevels(j) = true_num_levels;
    exp_t_01_95(j) = true_t_01_to_95;

end

%%
tic
std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
std_n = 1;
std_wig = 2;
sim_a = 0;
obj_thres = 10;
numSamples = 30;
for m = 1:length(mu_array)
    disp_text = sprintf("Processing %i of %i", m, length(mu_array));
    disp(disp_text)
    mu = mu_array(m);
    for n_i = 1:length(n_array)
        n = n_array(n_i);
        obj_exp = zeros(numExp, 1);
        [t_samples, n_samples, wig_samples] = deal(zeros(numSamples, 1));
        for i = 1:numSamples
            sim_fitness = simulator3(mu, n, sim_a, numGen, false);
            [t_01_to_95, wiggle, num_levels] = calcSummaryStats(sim_fitness);
            t_samples(i) = t_01_to_95;
            n_samples(i) = num_levels;
            wig_samples(i) = wiggle;
        end

        for j = 1:numExp
            true_t_01_to_95 = exp_t_01_95(j);
            true_num_levels = exp_numlevels(j);
            true_wiggle = exp_wiggles(j);
            obj_samples = zeros(1, numSamples);
            for s = 1:numSamples
                t_01_to_95 = t_samples(s);
                num_levels = n_samples(s);
                wiggle = wig_samples(s);

                wig_diff = abs(wiggle - true_wiggle);
                t_diff = abs(t_01_to_95 - true_t_01_to_95);
                n_diff = abs(num_levels - true_num_levels);
                obj_samples(s) = 1/std_n * n_diff + 1/std_time * t_diff + 1/std_wig * wig_diff;
            end
            obj_exp(j) = mean(obj_samples);
        end
        obj_ind = find(obj_exp<obj_thres);
        obj_meet_thres = obj_exp(obj_ind);
        obj_matrix(n_i, m) = mean(obj_meet_thres)*length(obj_exp)/length(obj_ind);
    end
end


total_time = toc
% signal end
load train
sound(y,Fs)


%% plot obj function
disp_obj_matrix = -obj_matrix;
mu_edges = [10^-9 10^-4];
n_edges = [1 10];
title_text = sprintf("Objective function value for mu = %.2g, n = %i. Model a = %.1f", true_mu, true_n, sim_a);
figure
imagesc('XData', mu_edges, 'YData', n_edges, 'CData', disp_obj_matrix)
title(title_text)
set(gca,'YDir','normal')
set(gca, 'xscale','log')
ylim([0.5 10.5])
xlim([0.5*10^-9 1.5*10^-4])

%% plot log obj function
obj_log = log(obj_matrix);
disp_obj_log =  -obj_log;
mu_edges = [10^-9 10^-4];
n_edges = [1 10];
title_text = sprintf("Objective function value for mu = %.2g, n = %i. Model a = %.1f", true_mu, true_n, sim_a);
figure
imagesc('XData', mu_edges, 'YData', n_edges, 'CData', disp_obj_log)
title(title_text)
set(gca,'YDir','normal')
set(gca, 'xscale','log')
ylim([0.5 10.5])
xlim([0.5*10^-9 1.5*10^-4])