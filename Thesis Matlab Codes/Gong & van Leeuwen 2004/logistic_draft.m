%con = [0,0,1,0;0,1,0,1;1,0,0,1;0,1,1,0];
p = 0.8; % probability of connection
n = 100; % number of nodes
N = 100; % number of iterations

% random matrix of connectivity
A = 1*(p>=rand(n));

g = digraph(randi([0 n], n))

x = zeros(n,N);

for i = 1:N
    x(:,i+1) = logisticgong2004(x(:,i),A);
end

plot(x(:,:))
