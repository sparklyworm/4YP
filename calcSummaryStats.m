function [t_01_to_95, wiggle, num_levels] = calcSummaryStats(fitness)
    inflpts =  findInflectionPoints2ndOrder(fitness);
    t_01 = computeTimeToPercentageMax(fitness, 0.01);
    t_95 = computeTimeToPercentageMax(fitness, 0.95);
    t_01_to_95 = t_95 - t_01;
    num_levels = round((length(inflpts) - 1)/2);
    
    t_20 = computeTimeToPercentageMax(fitness, 0.2);
    t_80 = computeTimeToPercentageMax(fitness, 0.8);
    wiggle = measureWigglyness(fitness, t_20, t_80);
end