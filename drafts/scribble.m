%%
arr = [ 1 2 3; 2 3 4];
A = [2 3 4; 1 2 3; 3 4 5; 2 3 3; 1 0 3];
[~, c] = ismember(arr, A, 'rows')

n = 2;
X = dec2bin(0:2^n-1)-'0';
X(X==0) = -1 ;


ncols = 4;
zloc = nchoosek(1:ncols,n)
m = size(zloc,1);
matn = zeros(m,ncols)
matn(sub2ind(size(matn),repmat((1:m)',[1,n]),zloc)) = 1

a = [ 1 0 1]
ind = find(a==1)

%%
% ploting binary at X axis
N = 4;
c = dec2bin(0:2^N-1);
s = cellstr(c)
C = categorical(s)
plot(C, [1: 16])
%%
% Calculate hamming distance between 2 binary strings
visited = {};
sum(s{1} ~= s{2})
if any(strcmp(s,'1000'))
    visited{end+1} = s{2};
    visited{end+1} = s{3};
end
%%
hold off
rectangle('Position',[0 0 4 2],'Curvature',0.9)
ylim([-1, 5])
xlim([-1, 5])
%%


