function infid = infidGate(w, w1, t, Ut, Ub)
%w is the nuclear Larmor frequency; w1 is n*3 matrix with [sqrt((azz1+w)^2+azx1^2), azz1+w/sqrt(), azx1/sqrt();...]; 
%t is the list of free evolution time; Ut is a series of target unitary
%operators; Ub are bases of unitary operators for the sake of computing
%efficiency.
%infid is the averaged infidelity

t = t(:);
tn = length(t);
n = size(w1,1);
fid = 0.25 * n;
for i = 1:n
    u0 = eye(2);
    u1 = eye(2);
    for j = 1:tn/2
        t1 = t(2*j-1);
        t2 = t(2*j);
        u0 = (cos(w1(i,1)*t2/2)*eye(2)-1i*sin(w1(i,1)*t2/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) *...
            [exp(-0.5i*w*t1),0;0,exp(0.5i*w*t1)] * u0;
        u1 = [exp(-0.5i*w*t2),0;0,exp(0.5i*w*t2)] *...
            (cos(w1(i,1)*t1/2)*eye(2)-1i*sin(w1(i,1)*t1/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1;
    end
    u0 = [exp(-0.5i*w*t(end)),0;0,exp(0.5i*w*t(end))] * u0;
    u1 = (cos(w1(i,1)*t(end)/2)*eye(2)-1i*sin(w1(i,1)*t(end)/2)*[w1(i,2),w1(i,3);w1(i,3),-w1(i,2)]) * u1;
    V = blkdiag(u0,u1);
    for j = 1:15
        fid = fid + 1/80 * trace(Ut(:,:,i)*Ub(:,:,j)'*Ut(:,:,i)'*V*Ub(:,:,j)*V');
    end      
end
infid = abs(1 - fid / n);

end

