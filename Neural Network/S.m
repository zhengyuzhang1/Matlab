function [ Sx ] = S( x, conjO, meanconjO, lambda)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Nmont = size(conjO,2);
a = conjO' * x;
b = meanconjO' * x;
Sx = conjO * a / Nmont - b * meanconjO ...
    + lambda * (mean(abs(conjO).^2,2) - abs(meanconjO).^2) .* x;

end

