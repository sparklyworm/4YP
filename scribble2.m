
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
n = 4;
mu = 5 * 10^-8;
numSamples = 10;
figure
for i = 1:numSamples
    sim_fitness = simulator(mu, n, a, false);
    hold on
    plot(1:600, sim_fitness, 'b');
end
%% 2
a = 0;
n = 8;
mu = 5 * 10^-5;
numSamples = 10;

for i = 1:numSamples
    sim_fitness = simulator(mu, n, a, false);
    hold on
    plot(1:600, sim_fitness, 'r');
end
%%
a = 0;
n = 4;
mu = 1 * 10^-6;
sim_fitness = simulator(mu, n, a, false);
inflpts =  findInflectionPoints2ndOrder(sim_fitness);
t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
t_01_to_99 = t_99 - t_01;
num_levels = round((length(inflpts) - 1)/2);

t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
total_e = measureWigglyness(sim_fitness, t_20, t_80)

title_text = sprintf("mu = %.2g, n = %i; num levels = %i, t-01-to-99 = %i", mu, n, num_levels, t_01_to_99);
figure
hold on
plot(1:600, sim_fitness, 'b')
title(title_text);
hold on
% plot(inflpts, sim_fitness(inflpts), 'ro')
plot([t_20 t_80], sim_fitness([t_20, t_80]), 'go')
hold on
plot([t_20 t_80], sim_fitness([t_20, t_80]), 'r')

%%
clc, clear
x = [ 11 12 13 14 15];
y = [ 1 2 3 4 5];
z = zeros(5);
x_edges = [11 15];
y_edges = [1 5];
for i = 1:length(y)
    for j = 1:length(y)
        z(j, i) = j*i;
    end
end

imagesc('XData', x_edges, 'YData', y_edges, 'CData',z)
set(gca,'YDir','normal')

%%
mu = 1 * 10^-8;
a = 0;
n = 4 ;
std_n = 1;
std_time = 25;

mu_samples= [0.5 1 5] * 10^-8;
n_samples = [ 3 4 5];
% obj_funcs = zeros(length(mu_samples) * length(n_samples));
wiggles = zeros(length(mu_samples) * length(n_samples), 1);
i = 1;
for m = 1:length(mu_samples)
    mu = mu_samples(m);
    for n_i = 1:length(n_samples)
        n = n_samples(n_i);
        sim_fitness = simulator(mu, n, a, false);
        %{
        t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
        t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
        t_01_to_99 = t_99 - t_01;
        t_diff = abs(t_01_to_99 - true_t_01_to_99);
        obj_funcs(i) = 1/std_n * abs(num_levels - true_n) + 1/std_time * t_diff;
        %}
        t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
        t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
        wiggles(i) = measureWigglyness(sim_fitness, t_20, t_80);
        i = i + 1;
    end
end

%%
y = [1 1.2 1.5 3.5 5]';
total_e = measureWigglyness(y, 1, 5);
