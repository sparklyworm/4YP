function [inflectionPoints, grad] = findInflectionPoints(y, threshold)
    % inflectionPoints: an array of indexes for x of where the inflection points occur
    grad = gradient(y);
    zeroFlag = true;
    inflpts = [0];

    for i = 1: length(grad)
       
        if zeroFlag 
            if grad(i) < threshold
                continue
            else
                zeroFlag = false;
                inflpts(end+1) = i;
            end
        else
            if grad(i) >= threshold
                continue
            else
                zeroFlag = true;
                inflpts(end+1) = i;
            end
        end
    
    end

    remove = zeros(1, length(inflpts));
    lastpt = 0;
    for i = 1:length(inflpts)
        if lastpt < inflpts(i) && inflpts(i) <= lastpt + 3
            remove(i-1)=lastpt;
            remove(i) = inflpts(i);
        end
        lastpt = inflpts(i);
    end

    inflectionPoints = setdiff(inflpts, remove);
    
end