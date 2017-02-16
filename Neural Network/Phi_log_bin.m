function [ phi_log ] = Phi_log_bin( alpha, theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
phi_log = alpha;
for i = 1:length(theta)
    if real(theta(i)) > 0
        phi_log = phi_log + theta(i) + log( 1 + exp(-theta(i)) ) ;
    else
        phi_log = phi_log + log( 1 + exp(theta(i)) );
    end
end

end

