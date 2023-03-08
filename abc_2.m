%% week 2
clc,clear
% true data
true_a = 0.3;
true_mu = 1 * 10 ^-8;
n = 4; % number of genes 
N = 10^8; % population size
numGen = 600;
selective_pressure = 1;

% gene_to_fitness =  assignFitness(n, true_a);    % maps genotype to fitness value
load("interesting-landscape\interesting-landscape-5.mat")
%%
[true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);
figure
plot(1:numGen, true_fitness)
% plotFitnessLandscape(gene_to_fitness, n, true)
[true_jumps, true_plateauLengths, true_numLevels] = computeSummaryStats(true_fitness);
%%
numSamples = 100;
mu_samples = zeros(numSamples, 1);
a_samples = zeros(numSamples,1 );
prior = makedist('Loguniform','Lower',10^-9,'Upper',10^-4);
total_samples = 0;
j_threshold = 5;
pl_threshold = 65;
tic
for i = 1:numSamples

    pl_d = pl_threshold + 1;
    j_d = j_threshold + 1;
    numLevels  = true_numLevels + 1;
    while ~(j_d < j_threshold && numLevels == true_numLevels && pl_d < pl_threshold)
        mu_s = random(prior,1);
        a_s = rand(1);
        sim_fitness = simulator(mu_s, a_s, true);
        try
            [jumps, plateauLengths, numLevels] = computeSummaryStats(sim_fitness);
            total_samples= total_samples+ 1;
        catch 
            continue
        end
        if numLevels ~= true_numLevels
            continue
        end
        j_dist = abs(jumps - true_jumps);
        j_d = sqrt(sum(j_dist.^2));
        pl_dist = abs(plateauLengths - true_plateauLengths);
        pl_d = sqrt(sum(pl_dist.^2));
        
    end
    mu_samples(i) = mu_s;
    a_samples(i) = a_s;
    % disp progress
    if mod(i, numSamples/10) ==0
        fprintf("sample %i \n", i)

    end
end
total_time = toc;
%%
[~,mu_edges] = histcounts(log10(mu_samples), 15);
figure
histogram(mu_samples,10.^mu_edges)
set(gca, 'xscale','log')

figure
[~,a_edges] = histcounts(a_samples, 15);
histogram(a_samples, 15)
%%
data = [mu_samples, a_samples];
figure
hist3(data,'CdataMode','auto', 'Edges', {10.^mu_edges, a_edges})
set(gca, 'xscale','log')
colorbar

%%
 plotFitnessLandscape(gene_to_fitness, n, true)