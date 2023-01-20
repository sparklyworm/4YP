clc, clear

x1_c1 = [6 7 5 4 8]';
x2_c1 = [1 2 3 1 3]';
X_c1 = [x1_c1, x2_c1];

x1_c2 = [1 2 3 2]';
x2_c2 = [5 6 5 3]';
X_c2 = [x1_c2, x2_c2];

t = [-1 * ones(length(x1_c1),1); ones(length(x1_c2),1)];
N = length(t);

%%
% lambdas
A = sym('a', [N 1]);
% weights
W = sym('w', [2 1]);
syms b
%lagrangian
