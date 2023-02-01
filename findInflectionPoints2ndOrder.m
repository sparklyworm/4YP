function inflectionPoints =  findInflectionPoints2ndOrder(y)
    % use 2nd order derivative to find inflection points
    % by finding cross over from postive to negative
    grad_1st_order = gradient(y);
    grad_2nd_order = gradient(grad_1st_order);
    
    negFlag = true;
    zeroFlag = false;
    inflpts = [0];  % x indices
    
    n1 = length(grad_1st_order);
    n2 = length(grad_2nd_order);
    thresholdPos = 0;
    thresholdNeg = -0.05 * 10^-3;
    for i = 1:n2
       
        grad2 = grad_2nd_order(i);
        if negFlag && grad2 >= thresholdPos && ~zeroFlag
            inflpts(end + 1) = i;
            negFlag = false;
            zeroFlag = true;
        elseif ~negFlag && grad2 <= thresholdNeg && zeroFlag
            inflpts(end + 1) = i;
            negFlag = true;
            zeroFlag = false;
        end
    end

    inflectionPoints = inflpts(2:end);

end