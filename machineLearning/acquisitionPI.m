function a = acquisitionPI(best_so_far, mu, x, variance)
    % x is a column vector 
    % mu is a column vector same size as x 
    % a is a column vector with values of the acquisition function
    % evalutaed at each point of x 

    n = length(x);
    a = zeros(n, 1);
    for i = 1:n
        a(i) = 1- normcdf(best_so_far, mu(i), variance(i));
    end
    
 
end