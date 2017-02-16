function [ H ] = totalz( N )
% total magnetization operator

sz = [1 0;0 -1];
dim = 2^N;
H = zeros(dim);
for i = 1:N
    H = H + kron(kron(eye(2^(i-1)),sz),eye(2^(N-i)));
end

