% find the average numLevels and t_01_to_99
clc, clear
% true data
true_a = 0.5;
true_mu = 1 * 10 ^-6;
true_n = 8; 
% population parameters
N = 10^8;
numGen = 600;
selective_pressure = 1;

numSamples = 100;
fitness_curves = cell(numSamples, 1);
times = zeros(numSamples, 1);
numLevels = zeros(numSamples, 1);
for i = 1:numSamples
    [sim_fitness, ~, ~] = simulator(true_mu, true_n, true_a, false);
 
    inflpts =  findInflectionPoints2ndOrder(sim_fitness);
    num_levels = round((length(inflpts) - 1)/2);
    t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
    t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
    t_01_to_99 = t_99 - t_01;

    fitness_curves{i} = sim_fitness;
    numLevels(i) = num_levels;
    times(i) = t_01_to_99;
end

%%
avg_time = mean(times)
std_time = std(times)
avg_numLevels = mean(numLevels)
std_numLevels = std(numLevels)