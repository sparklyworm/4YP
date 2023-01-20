function F = normalcdf(x, mu, sigma)
    % x is column vector
    % F is a column vector
    n = length(x);
   
    syms z 
    f = normalpdf(z, mu, sigma);
    %F_1 = vpaintegral(f, z,-inf,x(1));
    F_1 = int(f, [-inf, x(1)]);
    F = repmat(F_1, n, 1);
    %F = double(F);
end