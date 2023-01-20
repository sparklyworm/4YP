function fitnessLandscapeBones = constructFitnessLandscapeBones(jumps)
    n = length(jumps);
    fitnessLandscapeBones = ones(n+1, 1);
    for i = 2:n+1 
        fitnessLandscapeBones(i) = fitnessLandscapeBones(i-1) + jumps(i-1);
    end
end