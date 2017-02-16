function [x,u0,u1,axes,angle] = unequaltimeCPMG(w, azx, t, m)
%w is the nuclear Larmor frequency; azx is a n*2 matrix with
%[azz1,azx1;azz2,azx2;...]; t is the list of free evolution time; m is
%the number of nuclear that we want to have conditional rotation;
%x is the .

tn = length(t);
n = size(azx,1);
u0 = zeros(2,2,n);
u1 = zeros(2,2,n);
axes = zeros(3,n,2);
angle = zeros(2,n);
x = 0;
for i = 1:n
    u0(:,:,i) = eye(2);
    u1(:,:,i) = eye(2);
    for j = 1:(tn-1)/2
        u0(:,:,i) = expm(-0.5i*t(2*j)*((azx(i,1)+w)*[1,0;0,-1]+azx(i,2)*[0,1;1,0]))*expm(-0.5i*w*t(2*j-1)*[1,0;0,-1])*u0(:,:,i);
        u1(:,:,i) = expm(-0.5i*w*t(2*j)*[1,0;0,-1])*expm(-0.5i*t(2*j-1)*((azx(i,1)+w)*[1,0;0,-1]+azx(i,2)*[0,1;1,0]))*u1(:,:,i);
    end
    u0(:,:,i) = expm(-0.5i*w*t(tn)*[1,0;0,-1])*u0(:,:,i);
    u1(:,:,i) = expm(-0.5i*t(tn)*((azx(i,1)+w)*[1,0;0,-1]+azx(i,2)*[0,1;1,0]))*u1(:,:,i);
    [axes(:,i,1),angle(1,i)] = axisangle(u0(:,:,i));
    [axes(:,i,2),angle(2,i)] = axisangle(u1(:,:,i));
    %axes are x and z 
%     if i == m
%         x = x + (4-(axes(1,i,1)-axes(1,i,2))^2);
%     else
%         x = x + (4-(axes(3,i,1)+axes(3,i,2))^2);
%     end
end
%axes can be any direction 
tmp = sum(axes(:,:,1).*axes(:,:,2)) - ones(1,n);
tmp(m) = tmp(m) + 2;
x = sum(abs(tmp))/n;    

end

