function inflectionPoints = combineInflections(inflection1, inflection2)
    % inflection1 from tighter threshold
    inflect = [0];
    inflect2visited = [0];
    for i = 1: length(inflection1)
        inflect(end+1) = inflection1(i);
        upperb = inflection1(i) + 10;
        lowerb = inflection1(i)-10;
        for j = 1:length(inflection2)
            if inflection2(j) >= lowerb && inflection2(j) <= upperb
                inflect2visited(end+1) = inflection2(j);
                break
            end
        end
    end

    inflec2remain = setdiff(inflection2, inflect2visited);
    inflect = sort([inflect, inflec2remain]);
    inflectionPoints = inflect(2:end);


end