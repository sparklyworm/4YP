function finalE = expectedNormalOrder(r, n, mu, sigma)
    
    syms x
    f = normalpdf(x, mu, sigma);
    F = normalcdf(x, mu, sigma);
    G = x * ( 1 - F)^(r-1) * F^(n-r) * f;
    E = r * nchoosek(n, r)* vpaintegral(G, x, -inf, inf);
    try
        finalE = double(E); % numerial approximation 
    catch
        disp(E)
    end
end