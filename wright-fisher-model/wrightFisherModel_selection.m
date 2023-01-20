function dataSimulation = wrightFisherModel_selection(N, p0, s_a, mu_a, numGen, numSim)
    % A is the beneficial mutation 
    % N = population sizs
    % p0 = starting frequency of allele A
    % numGen = number of generations
    % numSim = number of simulations
    % mu_a = mutation rate of A 
    % s_a = selection coefficient of A relative to wildtype
    % r_a = relative fitness of A (wildtype fitness = 1) 
    % ref Avecilla 2022 

    r_a = 1 + s_a;

    
    dataSimulation = zeros(numSim, numGen);
    
    
    for j = 1: numSim
        p = p0;
        sim = zeros(1, numGen);
        for i = 1: numGen
            r_mean = r_a * p + 1*(1-p); % population mean fitness
            p = r_a * p/r_mean; % change in frequency due to natural selection 
            % y = binornd(N, p); 
            y = normrnd(N*p,sqrt(N*p*(1-p))); % change in frequency due to genetic drift
            % wt = N - y; % wild type 
            % new_a = wt * mu_a; 
            lambda = mu_a *N;
            new_a = poissrnd(lambda); % mutation from wild type poisson process
            y = y + new_a; 
            p = y/N;
            sim(i) = p;
        end
        dataSimulation(j,:) = sim;
        
    end
   
end




