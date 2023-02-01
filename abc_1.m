clc,clear
% true data
true_a = 0.3;
true_mu = 1 * 10 ^-8;
n = 4; % number of genes 
N = 10^8; % population size
numGen = 600;
selective_pressure = 1;

gene_to_fitness =  assignFitness(n, true_a);    % maps genotype to fitness value
[true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, true_mu, numGen, selective_pressure);

%%
numSamples = 1000;
mu_samples = zeros(numSamples, 1);
a_samples = zeros(numSamples,1 );
threshold = 10;
prior = makedist('Loguniform','Lower',10^-9,'Upper',10^-4);
tic
total_samples = 0;
for i = 1:numSamples

    distance = threshold + 1;
    while distance > threshold
        mu_s = random(prior,1);
        a_s = rand(1);
        sim_fitness = simulator(mu_s, a_s, false);
        distance = sqrt(sum((sim_fitness - true_fitness).^2));
        total_samples= total_samples+ 1;
    end
    mu_samples(i) = mu_s;
    a_samples(i) = a_s;
    % disp progress
    if mod(i, numSamples/10) ==0
        fprintf("sample %i \n", i)

    end
end
toc
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
plotFitnessLandscape(gene_to_fitness, true)
plotFitnessLandscapeGraph(gene_to_fitness)
%%
figure
plot(1:numGen, true_fitness)