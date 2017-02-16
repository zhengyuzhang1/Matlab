function [ transS ] = TransInv( walker )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(walker);
transS = repmat(walker,1,N);
for i = 2:N
    transS(:,i) = circshift(transS(:,i),i-1);
end

end

