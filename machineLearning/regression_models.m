clc, clear 
N = 10;
x_data = linspace(1, 10, N)';
x_plot = [1:0.1:max(x_data)]';
true_y_plot = sin(2*pi/10*x_plot) + 5 ;
true_y = sin(2*pi/10*x_data) + 5 ;
sigma_noise = 0.5;
y_data = true_y + normrnd(0, sigma_noise, N, 1);
% x_data = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]';
% y_data = [-4.9, -3.5, -2.8, 0.8, 0.3, -1.6, -1.3, 0.5, 2.1, 2.9, 5.6]';


%% linear regression with basis function
% = MLE assuming underlying Gaussian noise
% https://cookieblues.github.io//guides/2021/03/08/bsmalea-notes-1a/
% polynomial basis (d = 4)
% weights W = [w0
%              w1 
%              w2 
%              w3
%              w4] column vector
% basis X = [1, x1, x1^2, x1^3, x1^4
%            1, x2, x2^2, x2^3, x2^4
%            ...
%            1, xN, xN^2, xN^3, xN^4]
% model f = X*W

d = 4;
X = polybasis(x_data, d);
W = (X'*X)\X'* y_data; % least squares 
model_y = X * W;
X_plot = polybasis(x_plot,4);
model_f = X_plot * W;
plotModel(x_data, y_data, true_y_plot, model_y, model_f)

% plot polynomial basis functions
PolyBasis = cell(d+1, 1);
P = zeros(d+1, 1);
L = cell(d+1, 1);
for i = 1:d+1
    PolyBasis{i} = x_plot.^(i-1);
    hold on 
    legend_title = sprintf('x^%d', i-1);
    P(i) = plot(x_plot, PolyBasis{i} * W(i), 'DisplayName', legend_title);
    L{i} = legend_title;
end
legend(P, L);
ylim([min(y_data) - 5, max(y_data) + 5])

%% gaussian basis
% https://cookieblues.github.io/guides/2021/03/22/bsmalea-notes-2/
M = 6;
sigma = 1;
X = gaussianbasis(x_data, sigma, M);
W = (X'*X)\X'* y_data; % least squares 
X_plot = gaussianbasis(x_plot, sigma, M);
model_y = X*W;
model_f = X_plot*W;
plotModel(x_data, y_data, true_y_plot, model_y, model_f)

% plot Gaussian basis functions
G = cell(M,1);
x_domain = max(x_plot) - min(x_plot);
for i = 1:M
    mu_i = i/M * x_domain;
    G{i} = exp(-(x_plot - mu_i).^2/sigma^2);
    hold on 
    plot(x_plot, G{i} * W(i));
end




%% non-parametric regression: gaussian kernel regression
% number of parameters grows with training data 

lambda = 0.6;
model_y = gaussianKernel(x_data, x_data, y_data, lambda);
model_f = gaussianKernel(x_plot, x_data, y_data, lambda);
plotModel(x_data, y_data, true_y_plot, model_y, model_f)

% plot Guassian kernels
G_K = cell(length(x_data));
G_K_weight = zeros(length(x_plot), 1);
for i = 1:length(x_data)
    G_K{i} = normalpdf(x_plot, x_data(i),lambda);
    G_K_weight = G_K_weight + G_K{i};
end
for i = 1:length(x_data)
    hold on 
    p = plot(x_plot, y_data(i) * G_K{i}./G_K_weight);
    hold on
    plot(x_plot,G_K{i}, 'Color', p.Color, 'LineStyle','--')
end

%% jackknife cross validation to find best lambda
% works well when we have more training data
m = 100;
L = linspace(0.1, 1, m);
C = zeros(m, 1);
for i = 1:m
    C(i) = jackknifecv(x_data, y_data, L(i));
end
plot(L, C)
