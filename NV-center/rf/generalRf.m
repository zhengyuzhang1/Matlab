function [U0,U1] = generalRf(t, trf, phi, w, wrf, w1, rabi, b)

b = b/norm(b);
tn = length(t);
U0 = eye(2);
U1 = eye(2);
for i = 1:tn
    U0 = Urf(t(i),trf(i),phi(i),w,wrf,w1,rabi,b,mod(i+1,2)) * U0;
    U1 = Urf(t(i),trf(i),phi(i),w,wrf,w1,rabi,b,mod(i,2)) * U1;
end
        
end
