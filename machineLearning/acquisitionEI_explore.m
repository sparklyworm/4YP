function a = acquisitionEI_explore(best_so_far, mu, x, variance, expl)
    % x is a column vector 
    % mu is a column vector same size as x 
    % a is a column vector with values of the acquisition function
    % evalutaed at each point of x 
    % expl is exploration factor!!!! how much exploration we want 
    n = length(x);
    a = zeros(n, 1);
    for i = 1:n
        var = variance(i);
        a(i) = max((mu(i) - best_so_far) * normcdf(best_so_far, mu(i), var) + expl*var * normpdf(best_so_far, mu(i), var), 0);
    end
end