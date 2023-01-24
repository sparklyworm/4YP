% bivariate histogram
x = [ 1 1.4 2.3 2.2 3 3.1 3.4 3.2 3.1 4.44 4.2 5]';
y = [ 1 1 2.2 2.1 2.3 3.5 3 3.5 3.3 4 5 5]';
data = [x, y];
hist3(data,'CdataMode','auto', 'NBins', [5, 5], 'FaceColor','interp')
colorbar
% view(2)
%%
n = 4;
a = 0.5;
numSim = 1000;
simData = zeros(numSim, n+1);

for i = 1:numSim
    simpleFitnessLandscape = constructSimpleFitnessLandscape(n, a);
    simData(i, :) = simpleFitnessLandscape.values()';
end

%%
C = cell(4, 1);
c{1} = [C{1} "bbb"]
