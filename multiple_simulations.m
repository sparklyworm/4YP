clc, clear

n = 1; % number of genes !!!!!!!!!!!!!!!!!!!!!! 
a = 0.5; % random tune param
N = 10^8; % population size
mu = 1.1 *10^-8; % mutation rate

gene_to_fitness =  assignFitness(n, a);    % maps genotype to fitness value
genotypes = gene_to_fitness.keys();        % string array of genotypes
plotFitnessLandscape(gene_to_fitness, n, true);
% save("interesting-landscape.mat")
%%
selective_pressure = 1;
MU = [1] * 10^-8; % [ 1 5 10] * 10^-8
numSimulation = 1;
numGen = 100; % 600
numPlot = 1;
figure
for mu = MU
    disp(numPlot)
    for i = 1:numSimulation
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
        subplot(1, 3, numPlot)
        plot([1:numGen], fitness);
        title(sprintf("mutation rate %g", mu))
        hold on
    end
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
    numPlot = numPlot + 1;
end

