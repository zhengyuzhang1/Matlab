function [ s ] = State( X, N, M )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

a = X(1);
b = X(2:M+1);
W = reshape(X(M+2:end),M,N);
s = zeros(2^N,1);
for ibase = 1:2^N
    walker = 1-2*decimalToBinaryVector(ibase-1,N)';
    s(ibase) = Phi(a*sum(walker),repmat(b,1,N) + W * TransInv(walker));
end
s = s/norm(s);

end

