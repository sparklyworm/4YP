function [newInflectionSim, outlierIndex] = removeInflectionOutlier(numInflectionPoints, inflectionSim)

    [counts,lens] = groupcounts(numInflectionPoints);
    [count, index] = max(counts);
    avg_len = lens(index);
    newInflectionSim = cell(count, 1);
    outlierIndex = zeros(length(numInflectionPoints) - count, 1);
    j = 1; k = 1;
    for i = 1:length(inflectionSim)
        if length(inflectionSim{i}) == avg_len
            newInflectionSim{j} = inflectionSim{i};
            j = j + 1;
        else
            outlierIndex(k) = i;
            k = k + 1;
            continue
        end
    end

end