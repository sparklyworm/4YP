%{
Bishop 2006, Exercise 3.8 Bayesian Linear Regression

f = w0 + w1 * x
x ~ U(-1, 1) : x = unifrnd(-1,1)
W = [w0; w1]
mu_W = [0; 0]
cov_W = 1/alpha * eye(2)

%}

clc, clear
%%
% plot specs 
sample_line_color = "#ffdc00"; %#FF9500
true_line_color = 'r'; %#e45c00
sample_mean_color = 'b';
observed_point_color = 'g';
true_data_point_color = "#57625b";
sample_line_width = 1.5;
true_line_width = 2.2;
true_line_style = "-";
mean_line_style = '--';

%% prior distribution of W

alpha_w = 2;
mu_w0 = 0;
mu_w1 = 0;
sigma_w0 = 1/alpha_w;
sigma_w1 = 1/alpha_w;
sigma_w0w1 = 0;
w0_plot = [-1:0.1:1]';
w1_plot = [-1:0.1:1]';

Mu_W = [mu_w0; mu_w1];
Sigma_W = [sigma_w0 sigma_w0w1; sigma_w0w1 sigma_w1];
W = [w0_plot, w1_plot];

plotJointDistributionDual(W, Mu_W, Sigma_W)

%% true data
a0 = -0.3; a1 = 0.5; % true values for w0 and w1
n = 10;
x_data = linspace(-1, 1, n)';
true_y = a0 * ones(n, 1) + a1 * x_data;
sigma_noise = 0.2;
y_data = true_y + normrnd(0, sigma_noise, n, 1);
beta_noise = 1/sigma_noise^2;

%% data for plotting
N = 50;
x_plot = linspace(-1, 1, N)';
X_plot = [ones(N, 1) , x_plot];
%% prior: sample some w0 w1 values from prior
% and plot the lines (models) we would get using these w0 and w1
num_samples = 5;
prior_W = mvnrnd(Mu_W,Sigma_W, num_samples);
figure
for i = 1:num_samples
    w = prior_W(i, :)';
    y = X_plot * w;
    hold on
    plot(x_plot, y, 'Color', sample_line_color, 'LineWidth', sample_line_width)
end
hold on 
plot(x_data, true_y, 'Color', true_line_color, 'LineWidth', true_line_width);
xlabel('x')
ylabel('y')


%% basis functions
X = [ones(length(x_data), 1) x_data];

%% posterior P(w|y) after one/a few training data point
points = [4];
X1 = X(points, :);
Y1 = y_data(points);

inv_cov_w = alpha_w * eye(length(Sigma_W));
inv_cov_w_post = inv_cov_w + beta_noise * (X1'*X1);

cov_w_post = inv(inv_cov_w_post);
mu_w_post = inv_cov_w_post\(inv_cov_w * Mu_W + beta_noise * X1'*Y1 );

%% plot posterior distribution P(w|y)
plotJointDistributionDual(W, mu_w_post, cov_w_post)
hold on 
plot(a0, a1, 'r+', 'MarkerSize', 10, 'LineWidth', 3)
plot(mu_w_post(1), mu_w_post(2), 'b+', 'MarkerSize', 10, 'LineWidth', 3)


%% sample some w0 w1 values from posterior
% and plot the lines (models) we would get using these w0 and w1
num_samples = 5;
post_W = mvnrnd(mu_w_post,cov_w_post, num_samples);
figure
for i = 1:num_samples
    w = post_W(i, :)';
    y = X_plot * w;
    hold on
    plot(x_plot, y,'Color', sample_line_color, 'LineWidth', sample_line_width)
end
hold on 
plot(x_data, true_y, 'Color', true_line_color,'LineStyle', true_line_style,'LineWidth', true_line_width);
hold on
plot(x_data(points), Y1, 'Color', observed_point_color, 'Marker', ".", 'LineWidth', 3, 'MarkerSize',30)
xlabel('x')
ylabel('y')

%% plot mean of posterior predictive distribution using posterior mean of W
% = the mean of our estimate of y, using posterior mean of W
w = mu_w_post;
y_mean_post = X_plot * w; % mean of posterior predictive distribution
y_sigma_post = sqrt(1/beta_noise + diag(X_plot * cov_w_post * X_plot')) ;

top_sd = y_mean_post + y_sigma_post;
bottom_sd = y_mean_post - y_sigma_post;
hold on 
patch('Faces', linspace(1,2*N,2*N),'Vertices',[[x_plot; flip(x_plot)], [top_sd; flip(bottom_sd)]], 'FaceColor', 'yellow','EdgeColor', 'none','FaceAlpha', '0.4')
plot(x_plot, y_mean_post, 'Color', sample_mean_color, 'LineStyle', mean_line_style, 'LineWidth', true_line_width)

%% 
hold on
plot(x_data, y_data, 'Color', true_data_point_color, 'Marker', 'o', 'MarkerSize',6, 'LineWidth', 1, 'LineStyle','none')