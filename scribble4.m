N = 10^8;
true_a = 0;
true_n = 4;
true_mu = 1 * 10^-8;
%%
numGen2 = 2800;
[true2_fitness, ~] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen2, 1);
true2_inflpts =  findInflectionPoints2ndOrder(true2_fitness);
true2_t_01 = computeTimeToPercentageMax(true2_fitness, 0.01);
true2_t_99 = computeTimeToPercentageMax(true2_fitness, 0.99);
true2_t_01_to_99 = true2_t_99 - true2_t_01;
true2_num_levels = round((length(true2_inflpts) - 1)/2);

true2_t_20 = computeTimeToPercentageMax(true2_fitness, 0.2);
true2_t_80 = computeTimeToPercentageMax(true2_fitness, 0.8);
true2_wiggle = measureWigglyness(true2_fitness, true2_t_20, true2_t_80);

%%
sim_a = 0;
sim_n = 4;
sim_mu = 1 * 10^-8;
numGen2 = 2800;
numSamples = 20;
std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
std_n = 0.8;
std_wig = 2;
obj2_arr = zeros(numSamples, 1);
t2_arr = zeros(numSamples, 1);
n2_arr = zeros(numSamples, 1);
wig2_arr = zeros(numSamples, 1);
figure
for i = 1:numSamples
    sim_fitness = simulator3(sim_mu, sim_n, sim_a, numGen2, false);
    sim_inflpts =  findInflectionPoints2ndOrder(sim_fitness);
    sim_t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
    sim_t_99 = computeTimeToPercentageMax(sim_fitness, 0.99);
    sim_t_01_to_99 = sim_t_99 - sim_t_01;
    sim_num_levels = round((length(sim_inflpts) - 1)/2);
    
    sim_t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
    sim_t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
    sim_wiggle = measureWigglyness(sim_fitness, sim_t_20, sim_t_80);
    
    wig_diff = abs(sim_wiggle - true2_wiggle);
    t_diff = abs(sim_t_01_to_99 - true2_t_01_to_99);
    n_diff = abs(sim_num_levels - true2_num_levels);
    obj = 1/std_n * n_diff + 1/std_time * t_diff + 1/std_wig * wig_diff;

    t2_arr(i) = t_diff;
    n2_arr(i) = n_diff;
    wig2_arr(i) = wig_diff;
    obj2_arr(i) = obj;

    hold on
    plot(1:numGen2, sim_fitness);
end

avg2_obj = mean(obj2_arr);
avg2_t = mean(t2_arr);
avg2_n = mean(n2_arr);
avg2_wig = mean(wig2_arr);
