function X = gaussianbasis(x, sigma, M)
    
    % x = column vector of n data points
    % M = number of basis functions
    
    n = length(x);
    X = ones(n, M);
    x_domain = max(x) - min(x);
    
    for i = 1:M
        mu_i = i/M * x_domain * ones(n, 1);
        X(:,i) = exp(-(x - mu_i).^2/sigma^2); 
    end

end