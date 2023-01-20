clc, clear
load("interesting-landscape\interesting-landscape-5.mat")
hold off
%%

selective_pressure = 1;
mu = 1 * 10^-8;
numGen = 600;
[fitness, genotype_count_gen] = adaptiveWalk(gene_to_fitness, N, mu, numGen, selective_pressure);

% figure
hold on
plot([1:numGen], fitness);


%%
noisyFitness = zeros(numGen, 1);
for i = 1:numGen
    noisyFitness(i) = fitness(i) + 0.05*(-1 + 2*rand(1));
end

hold on
figure
plot(1:numGen, noisyFitness, "g")
%%
hold on
[smoothingSpline, goodness, output] = fit([1:numGen]', noisyFitness,'smoothingspline','SmoothingParam',0.5, 'Weights', [repmat(0.0005,1,200), repmat(0.000025,1,400)]);
plot(smoothingSpline);
% get value from smoothingSpline
y_smoothingSpline = feval(smoothingSpline, 1:numGen);

%%
hold on
cubicspline = csaps([1:numGen]', noisyFitness, 0.6,[],[repmat(0.0005,1,200), repmat(0.000025,1,400)]);
points = fnplt(cubicspline);
x = unique(points(1,:));
y = unique(points(2, :));
plot(points(1, :), points(2,:), 'DisplayName','Cubic Spline')

%%
finalInflection = findInflectionMultipleThreshold(y_smoothingSpline, [3 1.3]*10^-3);
hold on
plot(finalInflection, y_smoothingSpline(finalInflection), 'bo', 'LineWidth', 2);