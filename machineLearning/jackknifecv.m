function cv = jackknifecv(x_data, y_data, lambda)

    n = length(x_data);
    cv = 0;
    for i = 1:n 
        if i == 1
            x_i = x_data(2:end);
            y_i = y_data(2:end);
        elseif i == n
            x_i = x_data(1:end-1);
            y_i = y_data(1:end-1);
        else
            x_i = [x_data(1:i-1); x_data(i+1:end)];
            y_i = [y_data(1:i-1); y_data(i+1:end)];
        end
        y = gaussianKernel([x_data(i)], x_i, y_i, lambda);
        cv = cv + (y_data(i) - y)^2;
    end
    
    cv = cv/n;
end
