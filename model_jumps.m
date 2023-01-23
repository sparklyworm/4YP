clc, clear
load("interesting-landscape\interesting-landscape-5.mat")

%% setup true data; fit smooth spline; find inflectionPoints
selective_pressure = 0.8;
mu = 1 * 10^-8;
numGen = 600;
[true_fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);
noisyFitness = zeros(numGen, 1);
for i = 1:numGen
    noisyFitness(i) = true_fitness(i) + 0.05*(-1 + 2*rand(1));
end

% figure
% plot(1:numGen, noisyFitness, "g")
hold on
plot(1:numGen, true_fitness)

%% fit smoothing spline
hold on
[smoothingSpline, goodness, output] = fit([1:numGen]', noisyFitness,'smoothingspline','SmoothingParam',0.5, 'Weights', [repmat(0.0005,1,200), repmat(0.000025,1,400)]);
plot(smoothingSpline);

% get value from smoothingSpline
y_smoothingSpline = feval(smoothingSpline, 1:numGen);

%% inflection points
inflectionPoints = findInflectionMultipleThreshold(y_smoothingSpline, [3 1.3]*10^-3);
% inflectionPoints = findInflectionMultipleThreshold(true_fitness, [3 1.3]*10^-3);
hold on
plot(inflectionPoints, y_smoothingSpline(inflectionPoints), 'bo', 'LineWidth', 2);
% plot(inflectionPoints, true_fitness(inflectionPoints), 'bo', 'LineWidth', 2);
if rem(length(inflectionPoints), 2) ~= 0
    disp("more inflection points than expected, check and change!")
    disp("debug breakpoint")
end
%%
mu = 1.4 * 10^-8;
jumps = findJumpSizes(inflectionPoints, y_smoothingSpline);
fitnessLandscapeBones = constructFitnessLandscapeBones(jumps);
fitness =  modelAdaptiveWalk(fitnessLandscapeBones, N, mu, numGen, selective_pressure);
hold on 
plot(1:numGen, fitness, 'LineWidth', 2, 'Color', 'magenta', 'DisplayName','model')










