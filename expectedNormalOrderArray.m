function E = expectedNormalOrderArray(R, n, mu, sigma)
    
    E = ones(length(R), 1);
    for i = 1:length(R)
    
        E(i) = expectedNormalOrder(R(i), n, mu, sigma);
    
    end


end