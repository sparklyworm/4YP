function X = polybasis(x, d)
    % d = dimension of poly basis
    % x = column vector of n data points
    n = length(x);
    X = ones(n, d+1);
    
    for i = 2:d+1
        X(:,i) = x.^(i-1); 
    end
end