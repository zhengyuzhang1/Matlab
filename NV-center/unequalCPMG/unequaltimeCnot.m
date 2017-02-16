function [infid,u0,u1,axes,angle] = unequaltimeCnot(w, w1, t, m, n_repeat)
%w is the nuclear Larmor frequency; w1 is n*3 matrix with [sqrt((azz1+w)^2+azx1^2), azz1+w/sqrt(), azx1/sqrt();...]; 
%t is the list of free evolution time; m is
%the number of nuclear that we want to have conditional rotation;
%infid is the averaged infidelity, n_repeat is the number of repetation of
%the two-pulse unit;

tn = length(t);
n = size(w1,1);
u0 = zeros(2,2,n);
u1 = zeros(2,2,n);
axes = zeros(3,n,2);
angle = zeros(2,n);
infid = 0;
lastt = sum(t(2:2:end))-sum(t(1:2:end));
for i = 1:n
    u0(:,:,i) = eye(2);
    u1(:,:,i) = eye(2);
    for j = 1:tn/2
        t1 = t(2*j-1);
        t2 = t(2*j);
        u0(:,:,i) = (cos(w1(i,1)*t2/2)*eye(2)-1i*sin(w1(i,1)*t2/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) *...
            [exp(-0.5i*w*t1),0;0,exp(0.5i*w*t1)] * u0(:,:,i);
        u1(:,:,i) = [exp(-0.5i*w*t2),0;0,exp(0.5i*w*t2)] *...
            (cos(w1(i,1)*t1/2)*eye(2)-1i*sin(w1(i,1)*t1/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1(:,:,i);
    end
    u0(:,:,i) = [exp(-0.5i*w*lastt),0;0,exp(0.5i*w*lastt)] * u0(:,:,i);
    u1(:,:,i) = (cos(w1(i,1)*lastt/2)*eye(2)-1i*sin(w1(i,1)*lastt/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1(:,:,i);
    u0(:,:,i) = u0(:,:,i)^n_repeat;
    u1(:,:,i) = u1(:,:,i)^n_repeat;
    [axes(:,i,1),angle(1,i)] = axisangle(u0(:,:,i));
    [axes(:,i,2),angle(2,i)] = axisangle(u1(:,:,i));
    if i == m
        infid = infid + 1 - (abs(u0(1,1,i))^2+abs(u0(2,2,i))^2+abs(u1(1,2,i))^2+abs(u1(2,1,i))^2)/4;
    else
        infid = infid + ((1-axes(3,i,1)^2)*sin(angle(1,i)/2)^2+(1-axes(3,i,2)^2)*sin(angle(2,i)/2)^2)/2;
    end
end
infid = infid / n;
end

