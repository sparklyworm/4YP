a = 0;
n = 6;
mu = 1 * 10^-6;
numGen = 2800;
sim_fitness = simulator3(mu, n, a, numGen, false);
% hold on
% plot(1:numGen, sim_fitness)
%%
figure
numGen = length(experiments{1});
for i = 1:length(experiments)
     hold on
     plot(1:numGen, experiments{i})
%      [t_01_to_95, wiggle, num_levels] = calcSummaryStats(experiments{i})
end

% in manual_landscape_n_6, experiments 1 and 4 and 5 are the ones stuck 

%%
[obj, obj_arr] = calcObjFunc(experiments, sim_fitness)
%%
[TIME, W, NUMLVL] = calcSummaryStatsAll(experiments)
%%
[t, w, nlvl] = calcSummaryStats(sim_fitness)
%%
a = 0;
n = 6;
mu = 5 * 10^-6;
numGen = 2800;
numSamples = 5;
[T, W, L] = deal(zeros(1, numSamples));
% figure
for i = 1:numSamples
    sim_fitness = simulator3(mu, n, a, numGen, false);
    [t_01_to_95, wiggle, num_levels] = calcSummaryStats(sim_fitness);
    T(i) = t_01_to_95;
    W(i) = wiggle;
    L(i) = num_levels;
%     hold on
%     plot(1:numGen, sim_fitness, 'b')
end
%%
figure

[R, C] = ndgrid(1:size(disp_obj_log,1), 1:size(disp_obj_log,2));
R = R(:); C = C(:) - 1/4;
imagesc(disp_obj_matrix)
set(gca,'YDir','normal')
vals = disp_obj_log(:);
mask = vals >= -0.6;
% text(C(mask), R(mask), string(round(vals(mask), 2)), 'color', 'w')
text(C(mask), R(mask), string(round(vals(mask),2)), 'color', 'k')

%%
title(title_text)
mu_edges = [10^-9 10^-4];
n_edges = [1 10];
imagesc('XData', mu_edges, 'YData', n_edges, 'CData', disp_obj_matrix)
set(gca,'YDir','normal')
set(gca, 'xscale','log')
ylim([0.5 10.5])
xlim([0.5*10^-9 1.5*10^-4])
