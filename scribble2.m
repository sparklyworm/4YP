a = 0;
n = 4;
mu = 1 * 10^-8;
effectSize = 1/n;

[fitness, avg_gene_to_fitness, sim_genotype_count_gen] = simulator2(mu, n, a, effectSize, true);
plotModelAdaptiveWalk(sim_genotype_count_gen, mu, 10^6, false);
plotFitnessLandscape(avg_gene_to_fitness, true);
%%
inflpts =  findInflectionPoints2ndOrder(fitness);
figure
hold on
plot(1:600, fitness)
hold on
plot(inflpts, fitness(inflpts), 'ro')
%%
grad_1st_order = gradient(fitness);
grad_2nd_order = gradient(grad_1st_order);

plot(1:600, fitness)
hold on
plot(1:600, grad_1st_order*100 + 1, 'g')
hold on
plot(1:600, grad_2nd_order*1000, 'r')


%% 1
a = 0;
n = 6;
mu = 5 * 10^-6;
numSamples = 10;
% figure
for i = 1:numSamples
    sim_fitness = simulator3(mu, n, a, numGen, false);
    hold on
    plot(1:numGen, sim_fitness, 'b');
end
%% 2
a = 0;
n = 4;
mu = 5 * 10^-6;
numSamples = 10;

for i = 1:numSamples
    sim_fitness = simulator3(mu, n, a, numGen, false);
    hold on
    plot(1:numGen, sim_fitness, 'r');
end
%%
a = 0;
n = 4;
mu = 1 * 10^-8;
numGen = 800;

sim_fitness = simulator3(mu, n, a, numGen, false);
inflpts =  findInflectionPoints2ndOrder(sim_fitness);
t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
t_95 = computeTimeToPercentageMax(sim_fitness, 0.95);
t_01_to_95 = t_95 - t_01
num_levels = round((length(inflpts) - 1)/2);

t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
total_e = measureWigglyness(sim_fitness, t_20, t_80)


title_text = sprintf("mu = %.2g, n = %i; num levels = %i, t-01-to-99 = %i", mu, n, num_levels, t_01_to_95);

figure
plot(1:numGen, sim_fitness, 'b')
title(title_text);
hold on
% plot(inflpts, sim_fitness(inflpts), 'ro')

plot([t_01 t_95], sim_fitness([t_01, t_95]), 'go')
%hold on
%plot([t_20 t_80], sim_fitness([t_20, t_80]), 'r')
%%
figure
numGen = length(experiments{1});
for i = 1:length(experiments)
    hold on
    plot(1:numGen, experiments{i})
end

%%
a = 0;
n = 10;
mu = 1 * 10^-5;
numGen = 2800;

sim_fitness = simulator3(mu, n, a, numGen, false);
hold on
plot(1:numGen, sim_fitness, 'b')
obj = calcObjFunc(experiments, sim_fitness)

%%
fitness = experiments{1};
inflpts =  findInflectionPoints2ndOrder(fitness);
figure
plot(1:length(fitness), fitness);
hold on
plot(inflpts, fitness(inflpts), 'ro')

%% plot obj
obj_1 = obj_matrix;
obj_1(obj_1 > 60) = 60;
obj_1_log = log(obj_1);
disp_obj_1_log = - obj_1_log;
mu_edges = [10^-9 10^-4];
n_edges = [1 10];
title_text = sprintf("Objective function value for mu = %.2g, n = %i", true_mu, true_n);
figure
imagesc('XData', mu_edges, 'YData', n_edges, 'CData', disp_obj_1_log)
title(title_text)
set(gca,'YDir','normal')
set(gca, 'xscale','log')
ylim([0.5 10.5])
xlim([0.5*10^-9 1.5*10^-4])