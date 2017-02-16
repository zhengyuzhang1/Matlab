function [ v ] = RBMprob( X, N, M )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a = X(1:N);
b = X(N+1:N+M);
W = reshape(X(N+M+1:end),M,N);
n = 2^N;
v = zeros(n,1);
for i = 1:n
    state = flip(de2bi(i-1,N)');
    state = 1 - state;
    state(state==0) = -1;
    alpha = a.' * state; 
    theta = b + W * state;
    v(i) = Phi_log(alpha,theta);
end
v = v - min(real(v));
v = exp(v);
v = v / norm(v);
end

