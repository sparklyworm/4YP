clc, clear
load("interesting-landscape\interesting-landscape-5.mat")

% model parameters
numGen = 600;
mu = 1 * 10^-8;
selective_pressure = 1;


% simulation data
numSimulation = 10;
fitnessSim = cell(numSimulation,1);   % fitness curve of each simulation
inflectionSim = cell(numSimulation,1);% inflection points of each simulation
numInflPts = zeros(numSimulation, 1); % number of inflection points in each simulation
thresholds = [3 1.3] * 10^-3;         % gradient threshold for inflection points

for i = 1:numSimulation
        disp(i)
        [fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
        inflectionPoints = findInflectionMultipleThreshold(fitness, thresholds);
        % record
        numInflPts(i) = length(inflectionPoints);
        inflectionSim{i} = inflectionPoints;
        fitnessSim{i} = fitness;
end

%% tidy inflection points data
[inflectionSim_new, outlierIndex] = removeInflectionOutlier(numInflPts, inflectionSim);
inflectionSimMat = cell2mat(inflectionSim_new);

%% compute plateaus
plateauSim = cell(length(inflectionSim_new), 1);
for i = 1: length(inflectionSim_new)
    plateauLengths = computePlateauLength(inflectionSim_new{i});
    plateauSim{i} = plateauLengths;
end
plateauSimMat = cell2mat(plateauSim);

%% compute time to plateaus
time2plateauSim = cell(length(inflectionSim_new), 1);
for i = 1: length(inflectionSim_new)
    times = computeTimeToPlateau(inflectionSim_new{i});
    time2plateauSim{i} = times;
end
time2plateauSimMat = cell2mat(time2plateauSim);

%% plot histogram
plateauMeans = plotHistogram(plateauSimMat);

%% plot time to plateau
timeMeans = plotHistogram(time2plateauSimMat);

%% save variables
filename = sprintf("simulation_mu_%g_selective_%g.mat", mu, selective_pressure);
save(filename, "fitnessSim", "gene_to_fitness", "inflectionSim_new", "plateauSimMat", "plateauMeans", "mu", "selective_pressure");

