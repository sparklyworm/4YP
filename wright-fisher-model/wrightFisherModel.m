function dataSimulation = wrightFisherModel(N, p0, numGen, numSim)
    % N = population sizs
    % p0 = starting frequency of allele A
    % numGen = number of generations
    % numSim = number of simulations 
    
    dataSimulation = zeros(numSim, numGen);
    
    for j = 1: numSim
        p = p0;
        sim = zeros(1, numGen);
        for i = 1: numGen
            y = binornd(N, p);
            p = y/N;
            sim(i) = p;
        end
        dataSimulation(j,:) = sim;
    end
end




