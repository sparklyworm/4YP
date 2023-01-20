x_data = [0:1:4]';
m = 1;
y_true = m * x_data;
y_data = m * x_data + normrnd(0, 1, [length(x_data),1]);
% y_data = [0.4;1.7;2.8;2.3;4.4]
%%
figure
plot(x_data, y_data, 'rx', 'LineWidth', 2)
hold on 
plot(x_data, y_true, 'g', 'LineWidth',2)

%% 
M = [-0.5:0.05:1.5];
SSE = zeros(length(M), 1);
for i = 1: length(M)
    model_y = M(i) * x_data;
    SSE(i) = sum((model_y - y_data).^2);
end

%% ridge regression L2 
L = [0:10:50];
figure
for j = 1:length(L)
    SSE_L2 = zeros(length(M), 1);
    for i = 1: length(M)
        SSE_L2(i) = SSE(i) + L(j) * M(i)^2;
    end
    legend_title = sprintf("\\lambda = %i", L(j));
    plot(M, SSE_L2, "LineWidth", 2, "DisplayName",legend_title);
    hold on
end
legend
ylabel("SSE + L2 Ridge Regression")
xlabel("slope m")
%% lasso regression L1
figure
for j = 1:length(L)
    SSE_L1 = zeros(length(M), 1);
    for i = 1: length(M)
        SSE_L1(i) = SSE(i) + L(j) * abs(M(i));
    end
    legend_title = sprintf("\\lambda = %i", L(j));
    plot(M, SSE_L1, "LineWidth", 2, "DisplayName",legend_title);
    hold on
end
legend
ylabel("SSE + L1 Lasso Regression")
xlabel("slope m")