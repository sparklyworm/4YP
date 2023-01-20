function f = gaussianKernel(x, x_data, y_data, lambda)
    
    n = length(x);
    N = length(x_data);
    f = ones(n, 1);
    for i = 1:n
        X = x(i) * ones(N, 1);
        normalise_w = sum(normalpdf(X, x_data,lambda))/lambda;
        f(i) = 1/(lambda) * y_data' * normalpdf(X, x_data, lambda) /normalise_w;
    end
   
    

end