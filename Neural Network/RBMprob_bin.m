function [ v ] = RBMprob_bin( w, c )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n = size(w,1);
N = 2^n;
c = c(:);
v = zeros(N,1);
for i = 1:N
    state = flip(de2bi(i-1,n)');
    theta = c + w' * state;
    v(i) = exp(Phi_log_bin(0,theta));
end
v = v / sum(v);
end

