
%%
clc, clear
[fitness, ~, ~] = simulator(1 * 10^-4, 0.3, true);
figure
plot(1:600, fitness)

%%
inflectionPoints2 =  findInflectionPoints2ndOrder(fitness);
%%
figure
hold on
plot(1:600, fitness)
%%
hold on
plot(inflectionPoints2, fitness(inflectionPoints2), 'ro')
%%
grad_1st_order = gradient(fitness);
grad_2nd_order = gradient(grad_1st_order);
figure
plot(1:600, fitness)
hold on
plot(1:600, grad_1st_order*100 + 1, 'g')
hold on
plot(1:600, grad_2nd_order*1000, 'r')


%%
a = 0.2;
n = 4;
mu = 1 * 10^-8;

sim_fitness = simulator2(mu, n, a, false);
fitness = simulator2(mu, n, a, true);

hold on
plot(1:600, sim_fitness, 'm')
hold on
plot(1:600, fitness, 'c')

%%
