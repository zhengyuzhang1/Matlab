function [ phi ] = Phi( alpha, theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

phi = exp(alpha) * prod(cosh(theta(:)));

end

