function [ E ] = Efsf2d( N, gamma, lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = 0;
for i = 1:N
    for j = 1:N
        t = lambda - cos(2*pi*i/N) - cos(2*pi*j/N);
        delta = gamma * ( sin(2*pi*i/N) + sin(2*pi*i/N) );
        E = E - t - sqrt( t^2 + delta^2 );
    end
end
E = E / N^2;

end

