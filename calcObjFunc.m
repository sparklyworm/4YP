function [obj, obj_arr] = calcObjFunc(experiments, sim_fitness)
    numExp = length(experiments);
    sim_inflpts =  findInflectionPoints2ndOrder(sim_fitness);
    sim_t_01 = computeTimeToPercentageMax(sim_fitness, 0.01);
    sim_t_95 = computeTimeToPercentageMax(sim_fitness, 0.95);
    sim_t_01_to_95 = sim_t_95 - sim_t_01;
    sim_num_levels = round((length(sim_inflpts) - 1)/2);
    
    sim_t_20 = computeTimeToPercentageMax(sim_fitness, 0.2);
    sim_t_80 = computeTimeToPercentageMax(sim_fitness, 0.8);
    sim_wiggle = measureWigglyness(sim_fitness, sim_t_20, sim_t_80);

    std_time = 25; % 30: get slight overestimates for mu; 20 is nice but longer time 
    std_n = 0.8;
    std_wig = 2;
    obj_arr = zeros(numExp, 1);
    t_arr = zeros(numExp, 1);
    n_arr = zeros(numExp, 1);
    wig_arr = zeros(numExp, 1);
    
    
    exp_numlevels = zeros(numExp, 1);
    exp_t_01_99 = zeros(numExp, 1);
    exp_wiggles =  zeros(numExp, 1);
    for j = 1:numExp
        true_fitness = experiments{j};
        true_t_01 = computeTimeToPercentageMax(true_fitness, 0.01);
        true_t_95 = computeTimeToPercentageMax(true_fitness, 0.95);
        true_t_01_to_95 = true_t_95 - true_t_01;
        true_inflpts =  findInflectionPoints2ndOrder(true_fitness);
        true_num_levels = round((length(true_inflpts) - 1)/2);
    
        true_t_20 = computeTimeToPercentageMax(true_fitness, 0.2);
        true_t_80 = computeTimeToPercentageMax(true_fitness, 0.8);
        true_wiggle = measureWigglyness(true_fitness, true_t_20, true_t_80);
        
        wig_diff = abs(sim_wiggle - true_wiggle);
        t_diff = abs(sim_t_01_to_95 - true_t_01_to_95);
        n_diff = abs(sim_num_levels - true_num_levels);
        obj = 1/std_n * n_diff + 1/std_time * t_diff + 1/std_wig * wig_diff;
    
        t_arr(j) = t_diff;
        n_arr(j) = n_diff;
        wig_arr(j) = wig_diff;
        obj_arr(j) = obj;

        exp_wiggles(j) = true_wiggle;
        exp_numlevels(j) = true_num_levels;
        exp_t_01_99(j) = true_t_01_to_95;
    
    end
    
    obj = mean(obj_arr);

end