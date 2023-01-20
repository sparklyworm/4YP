function K = covRBF(x1, x2, sigma, lengthScale)
% covariance function using RBF (Gaussian kernel)
% x1 and x2 are column vectors of data points
% if x2 = [], then compute covariance of x1 with itself
    
    if isempty(x2)
        n = length(x1);
        L = zeros(n);
    
        for i = 1:n
            for j = 1:i
                L(i,j) = sigma^2 * exp(-0.5 * (x1(i) - x1(j))^2/lengthScale^2);
            end
        end
        K = L + L';
        K = K - diag(diag(L));
    else
        n1 = length(x1);
        n2 = length(x2);
        K = zeros(n1,n2);
    
        for i = 1:n1
            for j = 1:n2
                K(i,j) = sigma^2 * exp(-0.5 * (x1(i) - x2(j))^2/lengthScale^2);
            end
        end
      
       
    end

   
end