function [TIME, W, NUMLEVELS] = calcSummaryStatsAll(experiments)
    numExp = length(experiments);
    [TIME, W, NUMLEVELS]  = deal(zeros(numExp, 1));
    for i = 1:numExp
        [t_01_to_95, wiggle, num_levels] = calcSummaryStats(experiments{i});
        TIME(i) = t_01_to_95;
        W(i) = wiggle;
        NUMLEVELS(i) = num_levels;
    end
end