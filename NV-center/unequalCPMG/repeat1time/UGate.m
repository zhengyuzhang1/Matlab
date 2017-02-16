function U = UGate(w, w1, t)
%w is the nuclear Larmor frequency; w1 is n*3 matrix with [sqrt((azz1+w)^2+azx1^2), azz1+w/sqrt(), azx1/sqrt();...]; 
%t is the list of free evolution time;
%U is evolution operator for electron and all nuclei.
t = t(:);
tn = length(t);
n = size(w1,1);
U0 = 1;
U1 = 1;
for i = 1:n
    u0 = eye(2);
    u1 = eye(2);
    for j = 1:(tn-1)/2
        t1 = t(2*j-1);
        t2 = t(2*j);
        u0 = (cos(w1(i,1)*t2/2)*eye(2)-1i*sin(w1(i,1)*t2/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) *...
            [exp(-0.5i*w*t1),0;0,exp(0.5i*w*t1)] * u0;
        u1 = [exp(-0.5i*w*t2),0;0,exp(0.5i*w*t2)] *...
            (cos(w1(i,1)*t1/2)*eye(2)-1i*sin(w1(i,1)*t1/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1;
    end
    u0 = [exp(-0.5i*w*t(end)),0;0,exp(0.5i*w*t(end))] * u0;
    u1 = (cos(w1(i,1)*t(end)/2)*eye(2)-1i*sin(w1(i,1)*t(end)/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1;
    U0 = kron(U0, u0);
    U1 = kron(U1, u1);
end
U = blkdiag(U0, U1);
end

