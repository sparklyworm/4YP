function [jumps, plateauLengths, numLevels] = computeSummaryStats(fitness)
    
    inflectionPoints = findInflectionMultipleThreshold(fitness, [3 1.3]*10^-3);
    jumps = findJumpSizes(inflectionPoints, fitness);
    plateauLengths = computePlateauLength(inflectionPoints);
    numLevels = length(plateauLengths);
end