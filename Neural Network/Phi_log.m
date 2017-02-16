function [ phi_log ] = Phi_log( alpha, theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x1 = zeros(size(theta));
mask = abs(real(theta))>12;
x1(mask) = abs(real(theta(mask))) - log(2);
x1(~mask) = log(cosh(real(theta(~mask))));
phi_log = alpha + sum(x1 + log( cos(imag(theta)) + 1i*tanh(real(theta)).*sin(imag(theta))) ) ;

end

