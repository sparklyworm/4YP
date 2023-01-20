function a = acquisitionEI(best_so_far, mu, x, variance)
    % x is a column vector 
    % mu is a column vector same size as x 
    % a is a column vector with values of the acquisition function
    % evalutaed at each point of x 

    n = length(x);
    a = zeros(n, 1);
    for i = 1:n
        var = variance(i);
        a(i) = max((mu(i) - best_so_far) * normcdf(best_so_far, mu(i), var) + var * normpdf(best_so_far, mu(i), var), 0);
    end
end