function [ T ] = Trans1D( N )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

x = 0:2^N-1;
T = sparse(x+1, 2^(N-1)*mod(x,2)+floor(x/2)+1,1);

end

