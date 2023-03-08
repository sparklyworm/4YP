function simpleFitnessLandscape = constructSimpleFL2(n, a, effectSize)
    genotypes = [1:1:n+1]';
    simpleFitnessLandscape = dictionary(genotypes, ones(n+1, 1));

    for i = 2:n+1
        num_mutated = i-1;
        r = -1 + (1-(-1))* rand(1); % r is random number between -1 and 1
        s = effectSize * num_mutated * exp(a*r); 
        simpleFitnessLandscape(i) = 1 + s;
    end

end